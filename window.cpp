#include <QPainter>
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

    m_nodes  = nullptr;
    m_fvals  = nullptr;
    m_coeffs = nullptr;

    change_func();
}

Window::~Window()
{
    free(m_nodes);
    free(m_fvals);
    free(m_coeffs);
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

    f = nullptr;

    rebuild_approx2();
    update();
}

void Window::rebuild_approx2()
{
    if (n < 4) return;

    free(m_nodes);
    free(m_fvals);
    free(m_coeffs);

    m_nodes  = (double *)malloc(n * sizeof(double));
    m_fvals  = (double *)malloc(n * sizeof(double));
    m_coeffs = (double *)malloc(4 * (n - 1) * sizeof(double));

    for (int i = 0; i < n; i++) {
        m_nodes[i] = a + i * (b - a) / (n - 1);
        m_fvals[i] = func(func_id, m_nodes[i]);
    }

    build_piecewise_cubic(n, m_nodes, m_fvals, m_coeffs);
}

void Window::paintEvent(QPaintEvent * /* event */)
{
    QPainter painter(this);

    double x1, x2, y1, y2;
    double max_y, min_y;
    double delta_x = (b - a) / n;

    max_y = min_y = 0;
    for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
        y1 = func(func_id, x1);
        if (y1 < min_y) min_y = y1;
        if (y1 > max_y) max_y = y1;

        if (m_coeffs) {
            double y_approx = eval_piecewise_cubic(n, m_nodes, m_coeffs, x1);
            if (y_approx < min_y) min_y = y_approx;
            if (y_approx > max_y) max_y = y_approx;
        }
    }
    double delta_y = 0.01 * (max_y - min_y);
    min_y -= delta_y;
    max_y += delta_y;

    double maxF = (fabs(max_y) > fabs(min_y)) ? fabs(max_y) : fabs(min_y);
    printf("max|F| = %.6e\n", maxF);
    fflush(stdout);

    painter.save();
    painter.translate(0.5 * width(), 0.5 * height());
    painter.scale(width() / (b - a), -height() / (max_y - min_y));
    painter.translate(-0.5 * (a + b), -0.5 * (min_y + max_y));

    QPen pen("black");
    pen.setWidth(0);
    painter.setPen(pen);

    QPen pen_f(Qt::blue);
    pen_f.setWidth(0);
    painter.setPen(pen_f);
    x1 = a; y1 = func(func_id, x1);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
        y2 = func(func_id, x2);
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2; y1 = y2;
    }
    x2 = b; y2 = func(func_id, x2);
    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

    if (m_coeffs) {
        QPen pen_a(Qt::red);
        pen_a.setWidth(0);
        painter.setPen(pen_a);
        int steps = 4 * n;
        double step = (b - a) / steps;
        x1 = a;
        y1 = eval_piecewise_cubic(n, m_nodes, m_coeffs, x1);
        for (int i = 1; i <= steps; i++) {
            x2 = a + i * step;
            y2 = eval_piecewise_cubic(n, m_nodes, m_coeffs, x2);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            x1 = x2; y1 = y2;
        }
    }

    pen.setWidth(0);
    pen.setColor("red");
    painter.setPen(pen);
    painter.drawLine(a, 0, b, 0);
    painter.drawLine(0, max_y, 0, min_y);

    painter.restore();

    painter.setPen("blue");
    painter.drawText(0, 20, f_name);
    painter.drawText(0, 40, QString("n=%1  max|F|=%2")
                     .arg(n).arg(maxF, 0, 'e', 4));
}



