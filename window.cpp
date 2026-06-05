#include "window.h"
#include "functions.h"
#include "method33.h"
#include "method43.h"

#include <QPainter>
#include <QKeyEvent>
#include <QResizeEvent>
#include <QDebug>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define DRAW_STEPS 2000

Window::Window(QWidget *parent)
    : QWidget(parent),
      m_a(-1.0), m_b(1.0), m_n(10), m_k(0),
      m_view_mode(0), m_scale_s(0), m_perturb_p(0),
      m_x(NULL), m_f(NULL),
      m_a33(NULL), m_a43(NULL),
      m_work33(NULL), m_work43(NULL),
      m_alloc_n(0)
{
    setMinimumSize(400, 300);
    setFocusPolicy(Qt::StrongFocus);
}

Window::~Window()
{
    free_arrays();
}

void Window::free_arrays()
{
    free(m_x);     m_x     = NULL;
    free(m_f);     m_f     = NULL;
    free(m_a33);   m_a33   = NULL;
    free(m_a43);   m_a43   = NULL;
    free(m_work33); m_work33 = NULL;
    free(m_work43); m_work43 = NULL;
    m_alloc_n = 0;
}

void Window::alloc_arrays(int n)
{
    if (n == m_alloc_n) return;
    free_arrays();
    if (n < 2) {
        m_alloc_n = 0;
        return;
    }
    m_x     = (double*)malloc(n * sizeof(double));
    m_f     = (double*)malloc(n * sizeof(double));
    m_a33   = (double*)malloc(4 * (n - 1) * sizeof(double));
    m_a43   = (double*)malloc(4 * (n - 1) * sizeof(double));
    m_work33 = (double*)malloc((3 * n + 4) * sizeof(double));
    m_work43 = (double*)malloc(4 * n * sizeof(double));
    m_alloc_n = n;
}

int Window::parse_command_line(int argc, char *argv[])
{
    if (argc < 5) return 1;

    bool ok = true;
    char *end;
    m_a = strtod(argv[1], &end); if (*end != '\0') ok = false;
    m_b = strtod(argv[2], &end); if (*end != '\0') ok = false;
    m_n = (int)strtol(argv[3], &end, 10); if (*end != '\0') ok = false;
    m_k = (int)strtol(argv[4], &end, 10); if (*end != '\0') ok = false;

    if (!ok) return 1;
    if (m_a >= m_b) return 1;
    if (m_n < 2) return 1;
    if (m_k < 0 || m_k >= FUNC_COUNT) return 1;

    rebuild();
    return 0;
}

void Window::change_func()
{
    m_k = (m_k + 1) % FUNC_COUNT;
    m_perturb_p = 0;
    rebuild();
    update();
}

void Window::get_xrange(double &xa, double &xb) const
{
    double mid = (m_a + m_b) * 0.5;
    double half = (m_b - m_a) * 0.5;
    double factor = pow(2.0, (double)m_scale_s);
    xa = mid - half / factor;
    xb = mid + half / factor;
}

double Window::max_abs_f() const
{
    int steps = 500;
    double maxv = 0.0;
    for (int i = 0; i <= steps; i++) {
        double t = m_a + (m_b - m_a) * i / steps;
        double v = fabs(func(m_k, t));
        if (v > maxv) maxv = v;
    }
    return maxv;
}

double Window::fval(double x) const
{
    double val = func(m_k, x);
    if (m_perturb_p != 0 && m_alloc_n >= 2) {
        int idx = m_n / 2;
        if (fabs(x - m_x[idx]) < 1e-15 * fabs(m_x[idx] + 1.0)) {
            double mxf = max_abs_f();
            val += m_perturb_p * 0.1 * mxf;
        }
    }
    return val;
}

void Window::rebuild()
{
    int n = m_n;
    double xa;
    double xb;
    get_xrange(xa, xb);

    alloc_arrays(n);
    get_xrange(xa, xb);
    if (!m_x || !m_f) return;

    for (int i = 0; i < n; i++) {
        m_x[i] = xa + (xb - xa) * i / (double)(n - 1);
        m_f[i] = func(m_k, m_x[i]);
    }

    if (m_perturb_p != 0) {
        int idx = n / 2;
        double mxf = max_abs_f();
        m_f[idx] = func(m_k, m_x[idx]) + m_perturb_p * 0.1 * mxf;
    }

    method33_build(n, m_x, m_f, m_a33, m_work33);

    method43_build(n, m_x, m_f, m_a43, m_work43);
}


struct GraphPoint {
    double x, y;
};

void Window::draw_graph(QPainter &painter)
{
    int W = width();
    int H = height();

    double xa, xb;
    get_xrange(xa, xb);

    int steps = DRAW_STEPS;

    GraphPoint *pts_f   = new GraphPoint[steps + 1];
    GraphPoint *pts_m33 = new GraphPoint[steps + 1];
    GraphPoint *pts_m43 = new GraphPoint[steps + 1];
    GraphPoint *pts_e33 = new GraphPoint[steps + 1];
    GraphPoint *pts_e43 = new GraphPoint[steps + 1];

    {
        double mxf_scaled = max_abs_f();
        for (int i = 0; i <= steps; i++) {
            double t = xa + (xb - xa) * i / steps;
            pts_f[i].x = t;
            double v = func(m_k, t);
            if (m_perturb_p != 0 && m_alloc_n >= 2) {
                int idx = m_n / 2;
                if (fabs(t - m_x[idx]) < 1e-14) {
                    v += m_perturb_p * 0.1 * mxf_scaled;
                }
            }
            pts_f[i].y = v;
        }
    }

    if (m_alloc_n >= 2) {
        for (int i = 0; i <= steps; i++) {
            double t = xa + (xb - xa) * i / steps;
            pts_m33[i].x = t;
            pts_m33[i].y = method33_eval(t, m_a, m_b, m_n, m_x, m_a33);
        }
        for (int i = 0; i <= steps; i++) {
            double t = xa + (xb - xa) * i / steps;
            pts_m43[i].x = t;
            pts_m43[i].y = method43_eval(t, m_a, m_b, m_n, m_x, m_a43);
        }
        for (int i = 0; i <= steps; i++) {
            pts_e33[i].x = pts_f[i].x;
            pts_e33[i].y = pts_m33[i].y - pts_f[i].y;
            pts_e43[i].x = pts_f[i].x;
            pts_e43[i].y = pts_m43[i].y - pts_f[i].y;
        }
    } else {
        for (int i = 0; i <= steps; i++) {
            pts_m33[i] = pts_m43[i] = pts_f[i];
            pts_e33[i].x = pts_e43[i].x = pts_f[i].x;
            pts_e33[i].y = pts_e43[i].y = 0.0;
        }
    }

    GraphPoint **show_pts = NULL;
    int show_count = 0;
    GraphPoint *show_arr[4];

    switch (m_view_mode) {
    case 0:
        show_arr[0] = pts_f;
        show_arr[1] = pts_m33;
        show_count = 2;
        break;
    case 1:
        show_arr[0] = pts_f;
        show_arr[1] = pts_m43;
        show_count = 2;
        break;
    case 2:
        show_arr[0] = pts_f;
        show_arr[1] = pts_m33;
        show_arr[2] = pts_m43;
        show_count = 3;
        break;
    case 3:
        show_arr[0] = pts_e33;
        show_arr[1] = pts_e43;
        show_count = 2;
        break;
    }
    show_pts = show_arr;

    double fmin = 1e300, fmax = -1e300;
    for (int g = 0; g < show_count; g++) {
        for (int i = 0; i <= steps; i++) {
            double v = show_pts[g][i].y;
            if (v < fmin) fmin = v;
            if (v > fmax) fmax = v;
        }
    }

    if (fmax - fmin < 1e-15) {
        double mid = (fmax + fmin) * 0.5;
        fmin = mid - 1e-7;
        fmax = mid + 1e-7;
    }

    double max_abs_val = fabs(fmin) > fabs(fmax) ? fabs(fmin) : fabs(fmax);

    printf("view=%d  k=%d %s  n=%d  s=%d  p=%d  max|F|=%.6g\n",
           m_view_mode, m_k, func_name(m_k), m_n, m_scale_s, m_perturb_p, max_abs_val);
    fflush(stdout);

    auto mapX = [&](double x) -> double {
        return (x - xa) / (xb - xa) * W;
    };
    auto mapY = [&](double y) -> double {
        return H - (y - fmin) / (fmax - fmin) * H;
    };

    painter.fillRect(0, 0, W, H, QColor(255, 255, 255));

    painter.setPen(QPen(QColor(180, 180, 180), 1));
    if (fmin <= 0.0 && fmax >= 0.0) {
        int y0 = (int)mapY(0.0);
        painter.drawLine(0, y0, W, y0);
    }
    if (xa <= 0.0 && xb >= 0.0) {
        int x0 = (int)mapX(0.0);
        painter.drawLine(x0, 0, x0, H);
    }

    QColor col_f   (200, 140, 0);
    QColor col_m33 (0, 80, 200);
    QColor col_m43 (0, 160, 60);
    QColor col_e33 (200, 30, 30);
    QColor col_e43 (130, 0, 200);

    QColor colors[4];
    const char *labels[4];
    switch (m_view_mode) {
    case 0:
        colors[0] = col_f;   labels[0] = "f(x)";
        colors[1] = col_m33; labels[1] = "M33";
        break;
    case 1:
        colors[0] = col_f;   labels[0] = "f(x)";
        colors[1] = col_m43; labels[1] = "M43";
        break;
    case 2:
        colors[0] = col_f;   labels[0] = "f(x)";
        colors[1] = col_m33; labels[1] = "M33";
        colors[2] = col_m43; labels[2] = "M43";
        break;
    case 3:
        colors[0] = col_e33; labels[0] = "err33";
        colors[1] = col_e43; labels[1] = "err43";
        break;
    }

    for (int g = 0; g < show_count; g++) {
        painter.setPen(QPen(colors[g], 2));
        bool first = true;
        QPoint prev;
        for (int i = 0; i <= steps; i++) {
            int px = (int)mapX(show_pts[g][i].x);
            int py = (int)mapY(show_pts[g][i].y);
            if (!first) {
                painter.drawLine(prev, QPoint(px, py));
            }
            prev = QPoint(px, py);
            first = false;
        }
        painter.setPen(colors[g]);
        painter.drawText(10, 20 + g * 20, QString(labels[g]));
    }

    if (m_n <= 100 && m_alloc_n >= 2) {
        painter.setPen(QPen(QColor(0, 0, 0, 180), 1));
        for (int i = 0; i < m_n; i++) {
            if (m_x[i] < xa || m_x[i] > xb) continue;
            int px = (int)mapX(m_x[i]);
            int py = (int)mapY(m_f[i]);
            painter.drawEllipse(QPoint(px, py), 3, 3);
        }
    }

    painter.setPen(QColor(30, 30, 30));
    painter.setFont(QFont("Monospace", 10));

    QString info = QString("k=%1 %2  n=%3  s=%4  p=%5  max|F|=%6")
        .arg(m_k)
        .arg(QString(func_name(m_k)))
        .arg(m_n)
        .arg(m_scale_s)
        .arg(m_perturb_p)
        .arg(max_abs_val, 0, 'g', 5);

    painter.drawText(10, H - 10, info);
    
    if (m_alloc_n >= 2) {
    	int err_steps = 1000;
    	double err33 = 0.0;
    	double err43 = 0.0;
    	for (int i = 0; i <= err_steps; i++) {
    	    double t = xa + (xb - xa) * i / err_steps;
    	    double fv = func(m_k, t);
    	    double e33 = fabs(method33_eval(t, m_a, m_b, m_n, m_x, m_a33) - fv);
    	    double e43 = fabs(method43_eval(t, m_a, m_b, m_n, m_x, m_a43) - fv);
    	    if (e33 > err33) err33 = e33;
    	    if (e43 > err43) err43 = e43;
    	}
    	QString err_str = QString("max|err33|=%1   max|err43|=%2").arg(err33, 0, 'e', 4).arg(err43, 0, 'e', 4);
    	painter.drawText(10, H - 46, err_str);
    } 

    const char *view_names[] = {
        "Mode: f + Method33",
        "Mode: f + Method43",
        "Mode: f + Method33 + Method43",
        "Mode: Errors"
    };
    painter.drawText(10, H - 28, QString(view_names[m_view_mode]));

    delete[] pts_f;
    delete[] pts_m33;
    delete[] pts_m43;
    delete[] pts_e33;
    delete[] pts_e43;
}

void Window::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    draw_graph(painter);
}

void Window::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_0:
        m_k = (m_k + 1) % FUNC_COUNT;
        m_perturb_p = 0;
        rebuild();
        break;

    case Qt::Key_1:
        m_view_mode = (m_view_mode + 1) % 4;
        break;

    case Qt::Key_2:
        m_scale_s++;
	rebuild();
	break;

    case Qt::Key_3:
        m_scale_s--;
	rebuild();
        break;

    case Qt::Key_4:
        m_n *= 2;
        m_perturb_p = 0;
        rebuild();
        break;

    case Qt::Key_5:
        if (m_n >= 4) {
            m_n /= 2;
            m_perturb_p = 0;
            rebuild();
        }
        break;

    case Qt::Key_6:
        m_perturb_p++;
        rebuild();
        break;

    case Qt::Key_7:
        m_perturb_p--;
        rebuild();
        break;

    default:
        QWidget::keyPressEvent(event);
        return;
    }
    update();
}

void Window::resizeEvent(QResizeEvent *event)
{
    QWidget::resizeEvent(event);
    update();
}
