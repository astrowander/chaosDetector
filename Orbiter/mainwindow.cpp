#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "restricted_problem.h"
#include "glutplot.h"
#include <QThread>
#include <QTime>
#include <QInputDialog>
#define maxThreads 8
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    numberOfCurrentThreads=0;
    connect(ui->calculateButton,SIGNAL(clicked()),this,SLOT(calculate()));
    connect(ui->addButton,SIGNAL(clicked()),this,SLOT(addSpecialOrbit()));
    connect(this, SIGNAL(threadsChanged(int)), ui->threadsNumber,SLOT(display(int)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

class MyThread : public QThread {
Q_OBJECT
private:
    double dt, tmax, tol, C0;
    int whichIndicator;
    bool evolution;
    bool finished;

public:
    vector<double> *xaxis, *yaxis;
    double CPU_time, indicator, x;

    MyThread(double _C0, double _x,double _dt, double _tmax, vector<double> &_xaxis, vector<double> &_yaxis, int _whichIndicator, double _tol, bool _evolution)
    {
        C0 = _C0;
        x =_x;
        dt=_dt;
        tmax=_tmax;
        xaxis=&_xaxis;
        yaxis=&_yaxis;
        whichIndicator=_whichIndicator;
        tol=_tol;
        evolution=_evolution;

        CPU_time=0;
        indicator=0;
        finished=false;
    }

    void run()
    {
        integrate(C0, x,dt,tmax,*xaxis,*yaxis, whichIndicator, tol, indicator, CPU_time, evolution);
        finished=true;
        exec();
        return;
    }
    bool isFinished()
    {
        return finished;
    }

signals:

    //vector<double> *xaxisSignal();
    //vector<double> *yaxisSignal();
    //double indicatorSignal();
    //double CPU_time_signal();
};

int searchPlaceForInsert(double d, vector<double> &x)
{
    int n = x.size();
    int first = 0;
    int last = n;
    int mid;

    if (n==0) return 0;
    else if (x[0]>d) return 0;
    else if (x[n-1]<d) return n;

    while (first<last)
    {
        mid = first +(last-first)/2;
        if (d<=x[mid])
        {
            last=mid;
        }
        else
        {
            first=mid+1;
        }
    }
    return last;
}

void MainWindow::calculate()
{
    vector<vector<double> > x_all, y_all;

    double tmax = ui->tmaxLineEdit->text().toDouble();
    double dt = ui->dtLineEdit->text().toDouble();
    double tol = ui->toleranceLineEdit->text().toDouble();
    double xmin, xmax, ymin, ymax;
    int whichIndicator=0;
    if (ui->radioButton->isChecked()) whichIndicator=1;
    if (ui->radioButton_2->isChecked()) whichIndicator=7;
    if (ui->radioButton_3->isChecked()) whichIndicator=2;
    if (ui->radioButton_4->isChecked()) whichIndicator=4;
    if (ui->radioButton_5->isChecked()) whichIndicator=6;

    bool autoBorder = ui->checkBox->isChecked();
    if (!autoBorder)
    {
        xmin = ui->xminLineEdit->text().toDouble();
        xmax = ui->xmaxLineEdit->text().toDouble();
        ymin = ui->yminLineEdit->text().toDouble();
        ymax = ui->ymaxLineEdit->text().toDouble();
    }

    int nOfSpecialOrbits = specialOrbits.size();

    if(ui->tabWidget->currentIndex()==0)
    {
        double generalTime;
        double C0 = ui->C0LineEdit->text().toDouble();
        double x = ui->x0LineEdit->text().toDouble();
        double xm = ui->xmaxLineEdit->text().toDouble();
        int numberOfSteps = ui->dxLineEdit->text().toInt();
        double dx = (xm-x)/(numberOfSteps-1);
        ui->newdxLineEdit->setText(QString::number(dx));
        vector<double> xaxis = {}, yaxis = {};
        ui->progressBar->setMaximum(numberOfSteps);

       /* for (int i=0; i<(numberOfSteps+1)/8; i++)
        {
            double indicator[8], CPU_time[8], maxTime=0;
            ui->progressBar->setValue(i);
            MyThread* thread[8];
            for(int j=0; j<8; j++)
            {
                thread[j] = new MyThread(x0+i*dx*8 + j*dx,dt,tmax,xaxis,yaxis, whichIndicator, tol, false);
                thread[j]->start();
                emit threadsChanged(ui->threadsNumber->value()+1);
            }
            bool allFinished;
            do
            {
                allFinished=true;
                for (int j=0; j<8; j++) allFinished=allFinished && thread[j]->isFinished();
            } while (!allFinished);


            for (int j=0; j<8; j++)
            {
                xaxis.push_back(x0+i*dx*8 + j*dx);
                yaxis.push_back(thread[j]->indicator);
                if (CPU_time[j]>maxTime) maxTime=thread[j]->CPU_time;
            }
            generalTime+=maxTime/8;
        }
*/
        vector<MyThread*> threads;
        QTime cputime;
        double numberOfFinishedThreads=0;
        cputime.start();
        int countOfProcessedOrbits;
        do {
            while (threads.size()<maxThreads && x<xm)
            {
                if (countOfProcessedOrbits<nOfSpecialOrbits)
                {
                    x=specialOrbits[countOfProcessedOrbits++];
                }
                else if (countOfProcessedOrbits==nOfSpecialOrbits)
                {
                    x = ui->x0LineEdit->text().toDouble();
                    countOfProcessedOrbits++;
                }
                MyThread* temp = new MyThread(C0, x,dt,tmax,xaxis,yaxis, whichIndicator, tol, false);
                threads.push_back(temp);
                threads.back()->start();
                if (countOfProcessedOrbits > nOfSpecialOrbits) x+=dx;
            }

            for (int j=0; j<threads.size(); j++)
            {
                if (threads[j]->isFinished())
                {
                    int place = searchPlaceForInsert(threads[j]->x,xaxis);
                    xaxis.insert(xaxis.begin() + place,threads[j]->x);
                    yaxis.insert(yaxis.begin() + place,threads[j]->indicator);
                    threads[j]->quit();
                    threads[j]->wait();
                    threads.erase(threads.begin() + j);
                    //emit threadsChanged(ui->threadsNumber->value() + 1);
                    //numberOfFinishedThreads++;
                    ui->progressBar->setValue(++numberOfFinishedThreads);
                }
            }
        } while (x<=xm || threads.size()>0);
        int nMilliseconds = cputime.elapsed();
        ui->cputimeLineEdit->setText(QString::number((nMilliseconds+0.0)/1000));
        findmax(xaxis);
        findmax(yaxis);
        x_all.push_back(xaxis);
        y_all.push_back(yaxis);
        cout<< "General computation time: " << (nMilliseconds+0.0)/1000 << " seconds" << endl;

        char** temp;
        plot(0, temp, x_all, y_all, "", "", true, false, true, xmin, xmax, ymin, ymax, autoBorder);
    }

}

void MainWindow::oneThreadFinished()
{
    numberOfCurrentThreads--;
    emit threadsChanged(ui->threadsNumber->value()-1);
}

void MainWindow::addSpecialOrbit()
{
     specialOrbits.push_back(QInputDialog::getDouble(this, "Добавить орбиту", "Введите начальное условие x = ",0,-1,1,8));
     ui->listWidget->addItem( QString::number(specialOrbits.back()));

}

#include "mainwindow.moc"
