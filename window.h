#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QString>

class Window : public QWidget
{
    Q_OBJECT

public:
    explicit Window(QWidget *parent = 0);
    ~Window();

    int parse_command_line(int argc, char *argv[]);

public slots:
    void change_func();

protected:
    void paintEvent(QPaintEvent *event) Q_DECL_OVERRIDE;
    void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
    void resizeEvent(QResizeEvent *event) Q_DECL_OVERRIDE;

private:
    double  m_a;
    double  m_b;
    int     m_n;
    int     m_k;

    int     m_view_mode;
    int     m_scale_s;
    int     m_perturb_p;

    double *m_x;
    double *m_f;
    double *m_a33;
    double *m_a43;
    double *m_work33;
    double *m_work43;

    int     m_alloc_n;

    void rebuild();
    void free_arrays();
    void alloc_arrays(int n);

    void draw_graph(QPainter &painter);

    void get_xrange(double &xa, double &xb) const;

    double fval(double x) const;

    double max_abs_f() const;

    QString status_string() const;
};

#endif
