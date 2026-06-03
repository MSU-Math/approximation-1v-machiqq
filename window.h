#ifndef WINDOW_H
#define WINDOW_H
#include <QWidget>
#include <QString>
#include "approx2.h"

class Window : public QWidget
{
    Q_OBJECT
  private:
    int func_id;
    const char *f_name;
    double a;
    double b;
    int n;
    void (*f)(void);
    int view_mode;
    int scale_exp;
    int perturb;
    QString f_name;

    double *m_nodes;
    double *m_fvals;
    double *m_coeffs1;
    double *m_coeffs2;
    void rebuild();

  public:
    Window(QWidget *parent = nullptr);
    ~Window();
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    int parse_command_line(int argc, char *argv[]);
  public slots:
    void change_func();
    void change_view();
    void zoom_in();
    void zoom_out();
    void increase_n();
    void decrease_n();
    void increase_perturb();
    void decrease_perturb();

  protected:
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *);


};
#endif
