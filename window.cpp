#include "window.h"

#include "functions.h"
#include "method33.h"
#include "method43.h"

#include <QApplication>
#include <QKeyEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QFont>
#include <QFontMetrics>
#include <QRect>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>


struct GraphScale {
    double a, b;
    double Fmin, Fmax;
    int px, py, pw, ph;
};


static GraphScale draw_function_graph(QPainter &painter, const QRect &area,
                                       double a, double b,
                                       double (*func)(double),
                                       const QColor &color,
                                       const char *console_label)
{
    const int NSAMPLE = 2000;

    double Fmin = 1e300, Fmax = -1e300;
    std::vector<double> samples(NSAMPLE + 1);
    for (int i = 0; i <= NSAMPLE; ++i) {
        double t = a + (b - a) * i / NSAMPLE;
        double v = func(t);
        samples[i] = v;
        if (v < Fmin) Fmin = v;
        if (v > Fmax) Fmax = v;
    }

    if (Fmax - Fmin < 1e-12) {
        double mid = (Fmax + Fmin) / 2.0;
        Fmin = mid - 1.0;
        Fmax = mid + 1.0;
    }

    double maxabs = std::max(std::fabs(Fmin), std::fabs(Fmax));

    std::printf("%s: max|F| on [%.6g,%.6g] = %.10g\n", console_label, a, b, maxabs);

    QPen pen(color);
    pen.setWidth(2);
    painter.setPen(pen);

    double sx = area.width() / (b - a);
    double sy = area.height() / (Fmax - Fmin);

    QPointF prev;
    for (int i = 0; i <= NSAMPLE; ++i) {
        double t = a + (b - a) * i / NSAMPLE;
        double v = samples[i];
        double px = area.left() + (t - a) * sx;
        double py = area.top() + area.height() - (v - Fmin) * sy;
        QPointF cur(px, py);
        if (i > 0) {
            painter.drawLine(prev, cur);
        }
        prev = cur;
    }

    GraphScale gs;
    gs.a = a; gs.b = b; gs.Fmin = Fmin; gs.Fmax = Fmax;
    gs.px = area.left(); gs.py = area.top();
    gs.pw = area.width(); gs.ph = area.height();
    return gs;
}

static double nice_step(double range, int desired_ticks)
{
    if (range <= 0.0 || desired_ticks <= 0) {
        return 1.0;
    }
    double raw = range / desired_ticks;
    double mag = std::pow(10.0, std::floor(std::log10(raw)));
    double norm = raw / mag;

    double step;
    if (norm < 1.5) step = 1.0;
    else if (norm < 3.0) step = 2.0;
    else if (norm < 7.0) step = 5.0;
    else step = 10.0;

    return step * mag;
}

static void draw_axes(QPainter &painter, const QRect &area,
                       double a, double b, double Fmin, double Fmax)
{
    double sx = area.width() / (b - a);
    double sy = area.height() / (Fmax - Fmin);

    double py_axis;
    if (Fmin <= 0.0 && Fmax >= 0.0) {
        py_axis = area.top() + area.height() - (0.0 - Fmin) * sy;
    } else {
        py_axis = area.top() + area.height();
    }

    double px_axis;
    if (a <= 0.0 && b >= 0.0) {
        px_axis = area.left() + (0.0 - a) * sx;
    } else {
        px_axis = area.left();
    }

    QPen axisPen(QColor(90, 90, 90));
    axisPen.setWidth(1);
    painter.setPen(axisPen);
    painter.drawLine(QPointF(area.left(), py_axis), QPointF(area.right(), py_axis));
    painter.drawLine(QPointF(px_axis, area.top()), QPointF(px_axis, area.bottom()));

    QFont font = painter.font();
    font.setPointSize(8);
    painter.setFont(font);
    QFontMetrics fm(font);

    const int TICK_LEN = 5;
    char buf[64];

    double stepX = nice_step(b - a, 8);
    double startX = std::ceil(a / stepX) * stepX;
    for (double t = startX; t <= b + stepX * 1e-6; t += stepX) {
        double tClean = t;
        if (std::fabs(tClean) < stepX * 1e-9) tClean = 0.0;
        double px = area.left() + (t - a) * sx;
        painter.drawLine(QPointF(px, py_axis - TICK_LEN), QPointF(px, py_axis + TICK_LEN));
        std::snprintf(buf, sizeof(buf), "%.4g", tClean);
        int textW = fm.horizontalAdvance(buf);
        double ty = (py_axis + TICK_LEN + fm.ascent() + 1 <= area.bottom())
                        ? py_axis + TICK_LEN + fm.ascent() + 1
                        : py_axis - TICK_LEN - 2;
        painter.drawText(QPointF(px - textW / 2.0, ty), buf);
    }

    double stepY = nice_step(Fmax - Fmin, 8);
    double startY = std::ceil(Fmin / stepY) * stepY;
    for (double v = startY; v <= Fmax + stepY * 1e-6; v += stepY) {
        double vClean = v;
        if (std::fabs(vClean) < stepY * 1e-9) vClean = 0.0;
        double py = area.top() + area.height() - (v - Fmin) * sy;
        painter.drawLine(QPointF(px_axis - TICK_LEN, py), QPointF(px_axis + TICK_LEN, py));
        std::snprintf(buf, sizeof(buf), "%.4g", vClean);
        int textW = fm.horizontalAdvance(buf);
        double tx = (px_axis - TICK_LEN - textW - 2 >= area.left())
                        ? px_axis - TICK_LEN - textW - 2
                        : px_axis + TICK_LEN + 2;
        painter.drawText(QPointF(tx, py + fm.ascent() / 2.0 - 1), buf);
    }
}

struct Curve {
    const char *name;
    QColor color;
    double (*func)(double);
};

static void draw_curve_set(QPainter &painter, const QRect &area,
                            double a, double b,
                            Curve *curves, int ncurves,
                            double *out_maxabs, double *out_fmin, double *out_fmax)
{
    const int NSAMPLE = 1500;

    double Fmin = 1e300, Fmax = -1e300;
    std::vector<std::vector<double>> samples(ncurves, std::vector<double>(NSAMPLE + 1));

    for (int c = 0; c < ncurves; ++c) {
        for (int i = 0; i <= NSAMPLE; ++i) {
            double t = a + (b - a) * i / NSAMPLE;
            double v = curves[c].func(t);
            samples[c][i] = v;
            if (v < Fmin) Fmin = v;
            if (v > Fmax) Fmax = v;
        }
    }

    if (Fmax - Fmin < 1e-12) {
        double mid = (Fmax + Fmin) / 2.0;
        Fmin = mid - 1.0;
        Fmax = mid + 1.0;
    }

    double maxabs = std::max(std::fabs(Fmin), std::fabs(Fmax));
    if (out_maxabs) *out_maxabs = maxabs;
    if (out_fmin) *out_fmin = Fmin;
    if (out_fmax) *out_fmax = Fmax;

    draw_axes(painter, area, a, b, Fmin, Fmax);

    double sx = area.width() / (b - a);
    double sy = area.height() / (Fmax - Fmin);

    for (int c = 0; c < ncurves; ++c) {
        QPen pen(curves[c].color);
        pen.setWidth(2);
        painter.setPen(pen);

        QPointF prev;
        for (int i = 0; i <= NSAMPLE; ++i) {
            double t = a + (b - a) * i / NSAMPLE;
            double v = samples[c][i];
            double px = area.left() + (t - a) * sx;
            double py = area.top() + area.height() - (v - Fmin) * sy;
            QPointF cur(px, py);
            if (i > 0) {
                painter.drawLine(prev, cur);
            }
            prev = cur;
        }
    }

    std::printf("max{|Fmin|,|Fmax|} = %.10g  (Fmin=%.6g, Fmax=%.6g) on [%.6g,%.6g]\n",
                maxabs, Fmin, Fmax, a, b);
}

static Window *g_active_window = nullptr;

static double adapter_f(double t);
static double adapter_m1(double t);
static double adapter_m2(double t);
static double adapter_err1(double t);
static double adapter_err2(double t);

static double compute_max_abs(double a, double b, double (*func)(double))
{
    const int NSAMPLE = 2000;
    double m = 0.0;
    for (int i = 0; i <= NSAMPLE; ++i) {
        double t = a + (b - a) * i / NSAMPLE;
        double v = std::fabs(func(t));
        if (v > m) m = v;
    }
    return m;
}

Window::Window(QWidget *parent)
    : QWidget(parent), a0(-1.0), b0(1.0), n0(10), k0(0),
      k(0), n(10), mode(3), s(0), p(0)
{
    setFocusPolicy(Qt::StrongFocus);
    setWindowTitle("Graph");
}

Window::~Window()
{
    if (g_active_window == this) {
        g_active_window = nullptr;
    }
}

int Window::parse_command_line(int argc, char *argv[])
{
    if (argc != 5) {
        return 1;
    }

    char *endptr = nullptr;

    double aa = std::strtod(argv[1], &endptr);
    if (endptr == argv[1] || *endptr != '\0') return 1;

    double bb = std::strtod(argv[2], &endptr);
    if (endptr == argv[2] || *endptr != '\0') return 1;

    long nn = std::strtol(argv[3], &endptr, 10);
    if (endptr == argv[3] || *endptr != '\0') return 1;

    long kk = std::strtol(argv[4], &endptr, 10);
    if (endptr == argv[4] || *endptr != '\0') return 1;

    if (!(aa < bb)) return 1;
    if (nn < 3) return 1;
    if (kk < 0 || kk >= f_count()) return 1;

    a0 = aa;
    b0 = bb;
    n0 = (int)nn;
    k0 = (int)kk;

    n = n0;
    k = k0;
    mode = 3;
    s = 0;
    p = 0;

    g_active_window = this;
    rebuild_all();

    return 0;
}

void Window::current_view_range(double *a_view, double *b_view) const
{
    double mid = (a0 + b0) / 2.0;
    double half = (b0 - a0) / 2.0;
    double factor = std::pow(2.0, s);
    double half_view = half / factor;
    *a_view = mid - half_view;
    *b_view = mid + half_view;
}

void Window::rebuild_all()
{
    if (n < 3) n = 3;

    x.assign((size_t)n, 0.0);
    fv.assign((size_t)n, 0.0);

    for (long long i = 0; i < n; ++i) {
        double t = a0 + (b0 - a0) * (double)i / (double)(n - 1);
        x[(size_t)i] = t;
        fv[(size_t)i] = f_value(k, t);
    }

    if (n >= 1) {
        double maxabs_f = 0.0;
        for (long long i = 0; i < n; ++i) {
            double av = std::fabs(fv[(size_t)i]);
            if (av > maxabs_f) maxabs_f = av;
        }
        long long mid_idx = n / 2;
        fv[(size_t)mid_idx] += (double)p * 0.1 * maxabs_f;
    }

    coef33.assign((size_t)(4 * (n - 1)), 0.0);
    work33.assign((size_t)n, 0.0);
    method33_build((int)n, x.data(), fv.data(), coef33.data(), work33.data());

    coef43.assign((size_t)(4 * n), 0.0);
    work43.assign((size_t)(5 * (n + 1)), 0.0);
    method43_build((int)n, x.data(), fv.data(), coef43.data(), work43.data());
}

double Window::eval_f(double t) const
{
    return f_value(k, t);
}

double Window::eval_method1(double t) const
{
    return method33_value(t, a0, b0, (int)n, x.data(), coef33.data());
}

double Window::eval_method2(double t) const
{
    return method43_value(t, a0, b0, (int)n, x.data(), coef43.data());
}

static double adapter_f(double t)
{
    return g_active_window ? g_active_window->eval_f(t) : 0.0;
}
static double adapter_m1(double t)
{
    return g_active_window ? g_active_window->eval_method1(t) : 0.0;
}
static double adapter_m2(double t)
{
    return g_active_window ? g_active_window->eval_method2(t) : 0.0;
}
static double adapter_err1(double t)
{
    if (!g_active_window) return 0.0;
    return g_active_window->eval_method1(t) - g_active_window->eval_f(t);
}
static double adapter_err2(double t)
{
    if (!g_active_window) return 0.0;
    return g_active_window->eval_method2(t) - g_active_window->eval_f(t);
}

void Window::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event);

    g_active_window = this;

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.fillRect(rect(), Qt::white);

    QFont font = painter.font();
    font.setPointSize(10);
    painter.setFont(font);
    QFontMetrics fm(font);
    int textH = fm.height();
    int headerLines = 8;
    int headerHeight = textH * headerLines + 10;

    QRect plotArea(10, headerHeight, width() - 20, height() - headerHeight - 10);
    if (plotArea.width() < 10) plotArea.setWidth(10);
    if (plotArea.height() < 10) plotArea.setHeight(10);

    double a_view, b_view;
    current_view_range(&a_view, &b_view);
    {
        int thumbW = std::min(180, width() / 4);
        int thumbH = std::min(110, headerHeight - 4);
        if (thumbW > 20 && thumbH > 20) {
            QRect thumbArea(width() - thumbW - 8, 2, thumbW, thumbH);
            painter.setPen(QColor(150, 150, 150));
            painter.drawRect(thumbArea.adjusted(-1, -1, 0, 0));
            draw_function_graph(painter, thumbArea, a_view, b_view,
                                 adapter_f, QColor(0, 0, 0),
                                 "f(x) thumbnail (own scale)");
        }
    }

    double maxabs = 0.0;
    double fmin = 0.0, fmax = 0.0;

    switch (mode) {
    case 1: {
        Curve curves[2] = {
            { "f(x)", QColor(0, 0, 0), adapter_f },
            { "P1(x) [method33]", QColor(220, 30, 30), adapter_m1 },
        };
        draw_curve_set(painter, plotArea, a_view, b_view, curves, 2, &maxabs, &fmin, &fmax);
        break;
    }
    case 2: {
        Curve curves[2] = {
            { "f(x)", QColor(0, 0, 0), adapter_f },
            { "P2(x) [method43]", QColor(30, 60, 220), adapter_m2 },
        };
        draw_curve_set(painter, plotArea, a_view, b_view, curves, 2, &maxabs, &fmin, &fmax);
        break;
    }
    case 3: {
        Curve curves[3] = {
            { "f(x)", QColor(0, 0, 0), adapter_f },
            { "P1(x) [method33]", QColor(220, 30, 30), adapter_m1 },
            { "P2(x) [method43]", QColor(30, 60, 220), adapter_m2 },
        };
        draw_curve_set(painter, plotArea, a_view, b_view, curves, 3, &maxabs, &fmin, &fmax);
        break;
    }
    case 4: {
        Curve curves[2] = {
            { "P1(x)-f(x)", QColor(220, 30, 30), adapter_err1 },
            { "P2(x)-f(x)", QColor(30, 60, 220), adapter_err2 },
        };
        draw_curve_set(painter, plotArea, a_view, b_view, curves, 2, &maxabs, &fmin, &fmax);
        break;
    }
    default:
        break;
    }
    double err1 = compute_max_abs(a_view, b_view, adapter_err1);
    double err2 = compute_max_abs(a_view, b_view, adapter_err2);
    std::printf("max|P1(x)-f(x)| = %.10g,  max|P2(x)-f(x)| = %.10g   on [%.6g,%.6g]\n",
                err1, err2, a_view, b_view);

    painter.setPen(QColor(0, 0, 0));

    static const char *mode_names[5] = {
        "", "f + P1 (method 1)", "f + P2 (method 2)",
        "f + P1 + P2", "errors P1-f, P2-f"
    };

    char line[256];
    int ly = fm.ascent() + 2;

    std::snprintf(line, sizeof(line), "k=%d  %s", k, f_name(k));
    painter.drawText(10, ly, line); ly += textH;

    std::snprintf(line, sizeof(line), "n=%lld", n);
    painter.drawText(10, ly, line); ly += textH;

    std::snprintf(line, sizeof(line), "mode=%d (%s)", mode, mode_names[mode]);
    painter.drawText(10, ly, line); ly += textH;

    std::snprintf(line, sizeof(line), "scale s=%d   view=[%.6g,%.6g]", s, a_view, b_view);
    painter.drawText(10, ly, line); ly += textH;

    std::snprintf(line, sizeof(line), "perturbation p=%d", p);
    painter.drawText(10, ly, line); ly += textH;

    std::snprintf(line, sizeof(line), "max{|Fmin|,|Fmax|}=%.6g", maxabs);
    painter.drawText(10, ly, line); ly += textH;

    std::snprintf(line, sizeof(line), "max|P1-f| (method 1 error) = %.6g", err1);
    painter.drawText(10, ly, line); ly += textH;

    std::snprintf(line, sizeof(line), "max|P2-f| (method 2 error) = %.6g", err2);
    painter.drawText(10, ly, line); ly += textH;
}

void Window::keyPressEvent(QKeyEvent *event)
{
    bool handled = true;

    switch (event->key()) {
    case Qt::Key_0:
        change_func();
        break;

    case Qt::Key_1:
        mode = (mode % 4) + 1;
        update();
        break;

    case Qt::Key_2:
        s += 1;
        update();
        break;

    case Qt::Key_3:
        if (s > 0) s -= 1;
        update();
        break;

    case Qt::Key_4:
        n *= 2;
        rebuild_all();
        update();
        break;

    case Qt::Key_5:
        if (n / 2 >= 3) {
            n /= 2;
        }
        rebuild_all();
        update();
        break;

    case Qt::Key_6:
        p += 1;
        rebuild_all();
        update();
        break;

    case Qt::Key_7:
        p -= 1;
        rebuild_all();
        update();
        break;

    default:
        handled = false;
        break;
    }

    if (!handled) {
        QWidget::keyPressEvent(event);
    }
}

void Window::change_func()
{
    k = (k + 1) % f_count();
    rebuild_all();
    update();
}
