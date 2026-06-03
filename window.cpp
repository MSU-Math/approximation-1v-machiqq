#include <QPainter>
#include <QKeyEvent>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "window.h"
#include "functions.h"
#include "approx2.h"

#define DEFAULT_A -10
#define DEFAULT_B  10
#define DEFAULT_N  10

Window::Window(QWidget *parent) : QWidget(parent)
{
    a = DEFAULT_A;
    b = DEFAULT_B;
    n = DEFAULT_N;
    func_id = 0;
    view_mode = 0;
    scale_exp = 0;

    m_nodes  = nullptr;
    m_fvals  = nullptr;
    m_coeffs1 = nullptr;
    m_coeffs2 = nullptr;

    change_func();
}

Window::~Window()
{
    free(m_nodes);
    free(m_fvals);
    free(m_coeffs1);
    free(m_coeffs2);
}

QSize Window::minimumSizeHint() const { return QSize(100, 100); }
QSize Window::sizeHint()        const { return QSize(1000, 1000); }


int Window::parse_command_line(int argc, char *argv[])
{
    if (argc == 1)
        return 0;
    if (argc == 2)
        return -1;

    if (sscanf(argv[1], "%lf", &a) != 1 ||
        sscanf(argv[2], "%lf", &b) != 1 ||
        b - a < 1.e-6)
        return -2;

    if (argc > 3 && sscanf(argv[3], "%d", &n) != 1)
        return -2;
    if (n <= 0)
        return -2;

    if (argc > 4) {
        int k = 0;
        if (sscanf(argv[4], "%d", &k) != 1 || k < 0 || k > 6)
            return -2;
        func_id = k - 1;
    }

    change_func();
    return 0;
}

void Window::change_func()
{
    func_id = (func_id + 1) % 7;

    f_name = func_name(func_id);

    rebuild();
    update();
}

void Window::change_view()
{
    view_mode = (view_mode + 1) % 4;
    update();
}

void Window::zoom_out()
{
    scale_exp--;
    update();
}

void Window::zoom_in()
{
    scale_exp++;
    update();
}

void Window::increase_n()
{
    n *= 2;
    rebuild();
    update();
}

void Window::decrease_n()
{
    if (n <= 2) return;
    n /= 2;
    if (n < 2) n = 2;
    rebuild();
    update();
}

void Window::increase_perturb()
{
    perturb++;
    update();
}

void Window::decrease_perturb()
{
    perturb--;
    update();
}

void Window::rebuild()
{
    if (n < 4) return;

    free(m_nodes);
    free(m_fvals);
    free(m_coeffs1);
    free(m_coeffs2);

    m_nodes = nullptr;
    m_fvals = nullptr;
    m_coeffs1 = nullptr;
    m_coeffs2 = nullptr;

    m_nodes  = (double *)malloc(n * sizeof(double));
    m_fvals  = (double *)malloc(n * sizeof(double));
    m_coeffs1 = (double *)malloc(4 * (n - 1) * sizeof(double));
    m_coeffs2 = (double *)malloc(4 * (n - 1) * sizeof(double));

    double center = 0.5 * (a + b);
    double half = 0.5 * (b - a) / pow(2.0, scale_exp);
    double ca = center - half;
    double cb = center + half

    for (int i = 0; i < n; i++) {
        m_nodes[i] = ca + i * (cb - ca) / (n - 1);
        double fv = func(func_id, m_nodes[i]);
        if (perturb != 0 && i == n / 2) {
	    double maxF = 0.0;
	    int steps = 200;
	    for (int k = 0; k <= steps; k++) {
		double xx = ca + k * (cb - ca) / steps;
		double fabs_v = fabs(func(func_id, xx));
		if (fabs_v > maxF) maxF = fabs_v;
	    }
	    fv += perturb * 0.1 * maxF;
	}
        m_fvals[i] = fv;
    }
    if (n >= 4) {
        build_piecewise_cubic(n, m_nodes, m_fvals, m_coeffs1);
	build_piecewise_cubic(n, m_nodes, m_fvals, m_coeffs2);
    }
}

static void compute_range(int func_id, const double *nodes, int n, const double *c1,
                          const double *c2, int view_mode, double &min_y, double &max_y)
{
    min_y = 1e300;
    max_y = -1e300;

    int steps = 4 * n;
    double ca = nodes[0];
    double cd = nodes[n-1];

    for (int i = 0; i <= steps; i++) {
	double x = ca + i * (cb - ca) / steps;
	double yf = func(func_id, x);

        bool show_func = (view_mode == 0 || view_mode == 1 || view_mode == 2);
        bool show_approx = (view_mode == 0 || view_mode == 1 || view_mode == 2);
	bool show_err = (view_mode == 3);

	if (show_func) {
	    if (yf < min_y) min_y = yf;
	    if (yf > max_y) max_y = yf;

	}

        if (show_approx && n >= 4) {
	    if (view_mode == 0 || view_mode == 2) {
		double ya = eval_piecewise_cubic_m1(n, nodes, c1, x);
		if (ya < min_y) min_y = ya;
		if (ya > max_y) max_y = ya;
	    }
	    if (view_mode == 1 || view_mode == 2) {
		double ya = eval_piecewise_cubic(n, nodes, c2, x);
		if (ya < min_y) min_y = ya;
		if (ya > max_y) max_y = ya;
	    }
	}

	if (show_err && n >= 4) {
	    double e1 = 0;
	    double e2 = 0;
	    if (view_mode == 3) {
		e1 = fabs(yf - eval_piecewise_cubic_m1(n, nodes, c1, x));
		e2 = fabs(yf - eval_piecewise_cubic(n, nodes, c2, x));
	    }
	    if (e1 < min_y) min_y = e1;
	    if (e1 > max_y) max_y = e1;
	    if (e2 < min_y) min_y = e2;
	    if (e2 > max_y) max_y = e2;
	}
    }

    if (min_y > max_y) {
	min_y = -1.0;
	max_y = 1.0;
    }

    double dy = 0.01 * (max_y - min_y);
    double psi = 1e-12;

    if (dy < psi) dy = psi;
    min_y -=dy;
    max_y += dy;
}

void Window::keyPressEvent(QKeyEvent *e)
{
    switch (e -> key()) {
	case Qt::Key_0: change_func(); break;
	case Qt::Key_1: change_view(); break;
	case Qt::Key_2: zoom_in(); break;
	case Qt::Key_3: zoom_out(); break;
	case Qt::Key_4: increase_n(); break;
	case Qt::Key_5: decrease_n(); break;
	case Qt::Key_6: increase_perturb(); break;
	case Qt::Key_7: decrease_perturb(); break;
	default: QWidget::keyPressEvent(e); break;
    }
}

void Window::paintEvent(QPaintEvent *)
{
    QPainter painter(this);

    if (!m_nodes || n < 2) return;

    double center = 0.5 * (a + b);
    double half = 0.5 * (b - a) / pow(2.0, scale_exp);
    double ca = center - half;
    double cb = center + half;

    double max_y, min_y;
    compute_range(func_id, m_nodes, n, m_coeffs1, m_coeffs2, view_mode, min_y, max_y);
    double maxF = (fabs(max_y) > fabs(min_y)) ? fabs(max_y) : fabs(min_y);
    printf("max|F| = %.6e\n", maxF);
    fflush(stdout);

    painter.save();
    painter.translate(0.5 * width(), 0.5 * height());
    painter.scale(width() / (cb - ca), -height() / (max_y - min_y));
    painter.translate(-0.5 * (ca + cb), -0.5 * (min_y + max_y));

    int steps = 4 * n;
    double step = (cb - ca) / steps;

    if (view_mode != 3) {
	QPen pen_f(Qt::blue);
	pen_f.setWidth(0);
        painter.setPen(pen_f);
	double x1 = ca;
	double y1 = func(func_id, x1);
	for (int i = 1; i <= steps; i++) {
	    double x2 = ca + i * steps;
	    double y2 = func(func_id, x2);
	    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
	    x1 = x2;
	    y1 = y2;
	}
    }

    if (n >= 4) {
	if (view_mode == 0 || view_mode == 2) {
	    QPen pen_a(Qt::red);
	    pen_a.setWidth(0);
            painter.setPen(pen_a);
	    double x1 = ca;
	    double y1 = eval_piecewise_cubic_m1(n, m_nodes, m_coeffs1, x1);
	    for (int i = 1; i <= steps; i++) {
	        double x2 = ca + i * steps;
	        double y2 = eval_piecewise_cubic_m1(n, m_nodes, m_coeffs1, x2);
	        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
	        x1 = x2;
	        y1 = y2;
	    }
        }
	if (view_mode == 1 || view_mode == 2) {
	    QPen pen_b(Qt::darkGreen);
	    pen_b.setWidth(0);
            painter.setPen(pen_b);
	    double x1 = ca;
	    double y1 = eval_piecewise_cubic(n, m_nodes, m_coeffs2, x1);
	    for (int i = 1; i <= steps; i++) {
	        double x2 = ca + i * steps;
	        double y2 = eval_piecewise_cubic(n, m_nodes, m_coeffs2, x2);
	        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
	        x1 = x2;
	        y1 = y2;
	    }
        }

	if (view_mode == 3) {
	    QPen pen_e1(Qt::red);
	    pen_e1.setWidth(0);
            painter.setPen(pen_e1);
	    double x1 = ca;
	    double y1 = fabs(func(func_id, x1) - eval_piecewise_cubic_m1(n, m_nodes, m_coeffs1, x1));
	    for (int i = 1; i <= steps; i++) {
	        double x2 = ca + i * steps;
	        double y2 = fabs(func(func_id, x2) - eval_piecewise_cubic_m1(n, m_nodes, m_coeffs1, x2));
	        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
	        x1 = x2;
	        y1 = y2;
	    }

            QPen pen_e2(Qt::darkGreen);
	    pen_e2.setWidth(0);
            painter.setPen(pen_e2);
	    double x1 = ca;
	    double y1 = fabs(func(func_id, x1) - eval_piecewise_cubic(n, m_nodes, m_coeffs2, x1));
	    for (int i = 1; i <= steps; i++) {
	        double x2 = ca + i * steps;
	        double y2 = fabs(func(func_id, x2) - eval_piecewise_cubic(n, m_nodes, m_coeffs2, x2));
	        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
	        x1 = x2;
	        y1 = y2;
	    }
        }
    }
    QPen pen_axis(Qt::black);
    pen_axis.setWidht(0);
    painter.setPen(pen_axis);
    painter.drawLine(QPointF(ca,0), QPointF(cb, 0));
    painter.drawLine(QPointF(0, min_y), QPointF(0, max_y));
    painter.restore();

    double maxErr1 = 0.0;
    double maxErr2 = 0.0;
    if (n >= 4) {
	int errSteps = 1000;
	double ca2 = m_nodes[0];
	double cb2 = m_nodes[n-1];
	for (int i = 0; i <= errSteps; i++) {
	    double x = ca2 + i * (cb2 - ca2) / errSteps;
	    double fv = func(func_id, x);
	    double e1 = fabs(fv - eval_piecewise_cubic_m1(n, m_nodes, m_coeffs1, x));
	    double e2 = fabs(fv - eval_piecewise_cubic(n, m_nodes, m_coeffs2, x));
	    if (e1 > maxErr1) maxErr1 = e1;
	    if (e2 > maxErr2) maxErr2 = e2;
	}
    }

    painter.setPen("blue");
    painter.drawText(0, 20, QString("k = %1  %2").arg(func_id + 1).arg(f_name));
    painter.drawText(0, 40, QString("n = %1  s = %2  p = %3").arg(n).arg(scale_exp).arg(perturb));
    painter.drawText(0, 60, QString("max|F| = %1").arg(maxF, 0, 'e', 4));

    const char *modes[] = {"f + approx1","f + approx2", "f + approx1+2", "errors"};
    painter.drawText(0, 80, QString("view: %1").arg(modes[view_mode]));

    if ( n >= 4) {
	painter.drawText(0, 100, QString("err1 = %1  err2 = %2").arg(maxErr1, 0, 'e', 4).arg(maxErr2, 0, 'e', 4));
    }
}
