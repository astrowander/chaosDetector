#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private:
    Ui::MainWindow *ui;
    int numberOfCurrentThreads;
    std::vector<double> specialOrbits;

public slots:
    void calculate();
    void oneThreadFinished();
    void addSpecialOrbit();

signals:
    void threadsChanged(int);
};

#endif // MAINWINDOW_H
