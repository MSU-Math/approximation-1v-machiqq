#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <vector>


class Window : public QWidget
{
    Q_OBJECT

public:
    explicit Window(QWidget *parent = nullptr);
    ~Window();
    int parse_command_line(int argc, char *argv[]);

protected:
    void paintEvent(QPaintEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;

public slots:
    void change_func();

public:

    double eval_method1(double t) const;
    double eval_method2(double t) const;
    double eval_f(double t) const;

private:
    double a0, b0;
    int n0;
    int k0;

    int k;
    long long n;
    int mode;
    int s;
    int p;

    std::vector<double> x;
    std::vector<double> fv;

    std::vector<double> coef33;
    std::vector<double> work33;

    std::vector<double> coef43;
    std::vector<double> work43;

    void rebuild_all();

    void current_view_range(double *a_view, double *b_view) const;
};

#endif
