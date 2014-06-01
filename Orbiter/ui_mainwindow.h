/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Tue May 21 16:42:25 2013
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QAction>
#include <QApplication>
#include <QButtonGroup>
#include <QCheckBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLCDNumber>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QProgressBar>
#include <QPushButton>
#include <QRadioButton>
#include <QSpacerItem>
#include <QStatusBar>
#include <QTabWidget>
#include <QToolBar>
#include <QVBoxLayout>
#include <QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QTabWidget *tabWidget;
    QWidget *tab;
    QWidget *horizontalLayoutWidget_3;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_3;
    QLineEdit *x0LineEdit;
    QWidget *horizontalLayoutWidget_4;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_4;
    QLineEdit *dxLineEdit;
    QWidget *horizontalLayoutWidget_5;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_5;
    QLineEdit *xmaxLineEdit;
    QWidget *horizontalLayoutWidget_11;
    QHBoxLayout *horizontalLayout_11;
    QLabel *label_13;
    QLineEdit *newdxLineEdit;
    QLabel *label_15;
    QPushButton *addButton;
    QListWidget *listWidget;
    QWidget *tab_2;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QLineEdit *toleranceLineEdit;
    QHBoxLayout *horizontalLayout_10;
    QLabel *label_11;
    QLineEdit *tmaxLineEdit;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QLineEdit *dtLineEdit;
    QHBoxLayout *horizontalLayout_13;
    QLabel *label_16;
    QLineEdit *C0LineEdit;
    QLabel *label_6;
    QWidget *verticalLayoutWidget_2;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *checkBox;
    QSpacerItem *verticalSpacer;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_7;
    QLineEdit *xminLineEdit;
    QHBoxLayout *horizontalLayout_7;
    QLabel *label_8;
    QLineEdit *xmaxLineEdit_2;
    QHBoxLayout *horizontalLayout_9;
    QLabel *label_10;
    QLineEdit *yminLineEdit;
    QHBoxLayout *horizontalLayout_8;
    QLabel *label_9;
    QLineEdit *ymaxLineEdit;
    QPushButton *calculateButton;
    QLabel *label_12;
    QRadioButton *radioButton;
    QRadioButton *radioButton_2;
    QRadioButton *radioButton_3;
    QRadioButton *radioButton_4;
    QRadioButton *radioButton_5;
    QProgressBar *progressBar;
    QLCDNumber *threadsNumber;
    QWidget *horizontalLayoutWidget_12;
    QHBoxLayout *horizontalLayout_12;
    QLabel *label_14;
    QLineEdit *cputimeLineEdit;
    QMenuBar *menuBar;
    QMenu *menuOrbiter_0_2;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(687, 546);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setGeometry(QRect(10, 10, 421, 241));
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        horizontalLayoutWidget_3 = new QWidget(tab);
        horizontalLayoutWidget_3->setObjectName(QString::fromUtf8("horizontalLayoutWidget_3"));
        horizontalLayoutWidget_3->setGeometry(QRect(10, 10, 160, 31));
        horizontalLayout_3 = new QHBoxLayout(horizontalLayoutWidget_3);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_3 = new QLabel(horizontalLayoutWidget_3);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        horizontalLayout_3->addWidget(label_3);

        x0LineEdit = new QLineEdit(horizontalLayoutWidget_3);
        x0LineEdit->setObjectName(QString::fromUtf8("x0LineEdit"));

        horizontalLayout_3->addWidget(x0LineEdit);

        horizontalLayoutWidget_4 = new QWidget(tab);
        horizontalLayoutWidget_4->setObjectName(QString::fromUtf8("horizontalLayoutWidget_4"));
        horizontalLayoutWidget_4->setGeometry(QRect(10, 50, 160, 31));
        horizontalLayout_4 = new QHBoxLayout(horizontalLayoutWidget_4);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, 0, 0);
        label_4 = new QLabel(horizontalLayoutWidget_4);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        horizontalLayout_4->addWidget(label_4);

        dxLineEdit = new QLineEdit(horizontalLayoutWidget_4);
        dxLineEdit->setObjectName(QString::fromUtf8("dxLineEdit"));

        horizontalLayout_4->addWidget(dxLineEdit);

        horizontalLayoutWidget_5 = new QWidget(tab);
        horizontalLayoutWidget_5->setObjectName(QString::fromUtf8("horizontalLayoutWidget_5"));
        horizontalLayoutWidget_5->setGeometry(QRect(10, 90, 160, 31));
        horizontalLayout_5 = new QHBoxLayout(horizontalLayoutWidget_5);
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(0, 0, 0, 0);
        label_5 = new QLabel(horizontalLayoutWidget_5);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        horizontalLayout_5->addWidget(label_5);

        xmaxLineEdit = new QLineEdit(horizontalLayoutWidget_5);
        xmaxLineEdit->setObjectName(QString::fromUtf8("xmaxLineEdit"));

        horizontalLayout_5->addWidget(xmaxLineEdit);

        horizontalLayoutWidget_11 = new QWidget(tab);
        horizontalLayoutWidget_11->setObjectName(QString::fromUtf8("horizontalLayoutWidget_11"));
        horizontalLayoutWidget_11->setGeometry(QRect(10, 130, 160, 31));
        horizontalLayout_11 = new QHBoxLayout(horizontalLayoutWidget_11);
        horizontalLayout_11->setSpacing(6);
        horizontalLayout_11->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_11->setObjectName(QString::fromUtf8("horizontalLayout_11"));
        horizontalLayout_11->setContentsMargins(0, 0, 0, 0);
        label_13 = new QLabel(horizontalLayoutWidget_11);
        label_13->setObjectName(QString::fromUtf8("label_13"));

        horizontalLayout_11->addWidget(label_13);

        newdxLineEdit = new QLineEdit(horizontalLayoutWidget_11);
        newdxLineEdit->setObjectName(QString::fromUtf8("newdxLineEdit"));
        newdxLineEdit->setEnabled(true);

        horizontalLayout_11->addWidget(newdxLineEdit);

        label_15 = new QLabel(tab);
        label_15->setObjectName(QString::fromUtf8("label_15"));
        label_15->setGeometry(QRect(250, 10, 131, 17));
        addButton = new QPushButton(tab);
        addButton->setObjectName(QString::fromUtf8("addButton"));
        addButton->setGeometry(QRect(190, 170, 221, 29));
        listWidget = new QListWidget(tab);
        listWidget->setObjectName(QString::fromUtf8("listWidget"));
        listWidget->setGeometry(QRect(190, 40, 221, 121));
        tabWidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        tabWidget->addTab(tab_2, QString());
        verticalLayoutWidget = new QWidget(centralWidget);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(10, 260, 222, 136));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label = new QLabel(verticalLayoutWidget);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);

        toleranceLineEdit = new QLineEdit(verticalLayoutWidget);
        toleranceLineEdit->setObjectName(QString::fromUtf8("toleranceLineEdit"));

        horizontalLayout->addWidget(toleranceLineEdit);


        verticalLayout->addLayout(horizontalLayout);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setSpacing(6);
        horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
        label_11 = new QLabel(verticalLayoutWidget);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        horizontalLayout_10->addWidget(label_11);

        tmaxLineEdit = new QLineEdit(verticalLayoutWidget);
        tmaxLineEdit->setObjectName(QString::fromUtf8("tmaxLineEdit"));

        horizontalLayout_10->addWidget(tmaxLineEdit);


        verticalLayout->addLayout(horizontalLayout_10);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        label_2 = new QLabel(verticalLayoutWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout_2->addWidget(label_2);

        dtLineEdit = new QLineEdit(verticalLayoutWidget);
        dtLineEdit->setObjectName(QString::fromUtf8("dtLineEdit"));

        horizontalLayout_2->addWidget(dtLineEdit);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_13 = new QHBoxLayout();
        horizontalLayout_13->setSpacing(6);
        horizontalLayout_13->setObjectName(QString::fromUtf8("horizontalLayout_13"));
        label_16 = new QLabel(verticalLayoutWidget);
        label_16->setObjectName(QString::fromUtf8("label_16"));

        horizontalLayout_13->addWidget(label_16);

        C0LineEdit = new QLineEdit(verticalLayoutWidget);
        C0LineEdit->setObjectName(QString::fromUtf8("C0LineEdit"));

        horizontalLayout_13->addWidget(C0LineEdit);


        verticalLayout->addLayout(horizontalLayout_13);

        label_6 = new QLabel(centralWidget);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(480, 10, 151, 17));
        verticalLayoutWidget_2 = new QWidget(centralWidget);
        verticalLayoutWidget_2->setObjectName(QString::fromUtf8("verticalLayoutWidget_2"));
        verticalLayoutWidget_2->setGeometry(QRect(460, 40, 214, 211));
        verticalLayout_2 = new QVBoxLayout(verticalLayoutWidget_2);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        checkBox = new QCheckBox(verticalLayoutWidget_2);
        checkBox->setObjectName(QString::fromUtf8("checkBox"));

        verticalLayout_2->addWidget(checkBox);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        label_7 = new QLabel(verticalLayoutWidget_2);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        horizontalLayout_6->addWidget(label_7);

        xminLineEdit = new QLineEdit(verticalLayoutWidget_2);
        xminLineEdit->setObjectName(QString::fromUtf8("xminLineEdit"));

        horizontalLayout_6->addWidget(xminLineEdit);


        verticalLayout_2->addLayout(horizontalLayout_6);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        label_8 = new QLabel(verticalLayoutWidget_2);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        horizontalLayout_7->addWidget(label_8);

        xmaxLineEdit_2 = new QLineEdit(verticalLayoutWidget_2);
        xmaxLineEdit_2->setObjectName(QString::fromUtf8("xmaxLineEdit_2"));

        horizontalLayout_7->addWidget(xmaxLineEdit_2);


        verticalLayout_2->addLayout(horizontalLayout_7);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
        label_10 = new QLabel(verticalLayoutWidget_2);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        horizontalLayout_9->addWidget(label_10);

        yminLineEdit = new QLineEdit(verticalLayoutWidget_2);
        yminLineEdit->setObjectName(QString::fromUtf8("yminLineEdit"));

        horizontalLayout_9->addWidget(yminLineEdit);


        verticalLayout_2->addLayout(horizontalLayout_9);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        label_9 = new QLabel(verticalLayoutWidget_2);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        horizontalLayout_8->addWidget(label_9);

        ymaxLineEdit = new QLineEdit(verticalLayoutWidget_2);
        ymaxLineEdit->setObjectName(QString::fromUtf8("ymaxLineEdit"));

        horizontalLayout_8->addWidget(ymaxLineEdit);


        verticalLayout_2->addLayout(horizontalLayout_8);

        calculateButton = new QPushButton(centralWidget);
        calculateButton->setObjectName(QString::fromUtf8("calculateButton"));
        calculateButton->setGeometry(QRect(240, 260, 191, 101));
        label_12 = new QLabel(centralWidget);
        label_12->setObjectName(QString::fromUtf8("label_12"));
        label_12->setGeometry(QRect(470, 270, 181, 20));
        radioButton = new QRadioButton(centralWidget);
        radioButton->setObjectName(QString::fromUtf8("radioButton"));
        radioButton->setGeometry(QRect(470, 300, 117, 22));
        radioButton->setChecked(true);
        radioButton_2 = new QRadioButton(centralWidget);
        radioButton_2->setObjectName(QString::fromUtf8("radioButton_2"));
        radioButton_2->setGeometry(QRect(470, 330, 117, 22));
        radioButton_3 = new QRadioButton(centralWidget);
        radioButton_3->setObjectName(QString::fromUtf8("radioButton_3"));
        radioButton_3->setGeometry(QRect(470, 360, 117, 22));
        radioButton_4 = new QRadioButton(centralWidget);
        radioButton_4->setObjectName(QString::fromUtf8("radioButton_4"));
        radioButton_4->setGeometry(QRect(470, 390, 117, 22));
        radioButton_5 = new QRadioButton(centralWidget);
        radioButton_5->setObjectName(QString::fromUtf8("radioButton_5"));
        radioButton_5->setGeometry(QRect(470, 420, 117, 22));
        progressBar = new QProgressBar(centralWidget);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(10, 400, 411, 23));
        progressBar->setValue(0);
        threadsNumber = new QLCDNumber(centralWidget);
        threadsNumber->setObjectName(QString::fromUtf8("threadsNumber"));
        threadsNumber->setGeometry(QRect(10, 430, 81, 41));
        horizontalLayoutWidget_12 = new QWidget(centralWidget);
        horizontalLayoutWidget_12->setObjectName(QString::fromUtf8("horizontalLayoutWidget_12"));
        horizontalLayoutWidget_12->setGeometry(QRect(270, 440, 160, 31));
        horizontalLayout_12 = new QHBoxLayout(horizontalLayoutWidget_12);
        horizontalLayout_12->setSpacing(6);
        horizontalLayout_12->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_12->setObjectName(QString::fromUtf8("horizontalLayout_12"));
        horizontalLayout_12->setContentsMargins(0, 0, 0, 0);
        label_14 = new QLabel(horizontalLayoutWidget_12);
        label_14->setObjectName(QString::fromUtf8("label_14"));

        horizontalLayout_12->addWidget(label_14);

        cputimeLineEdit = new QLineEdit(horizontalLayoutWidget_12);
        cputimeLineEdit->setObjectName(QString::fromUtf8("cputimeLineEdit"));
        cputimeLineEdit->setEnabled(true);

        horizontalLayout_12->addWidget(cputimeLineEdit);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 687, 25));
        menuOrbiter_0_2 = new QMenu(menuBar);
        menuOrbiter_0_2->setObjectName(QString::fromUtf8("menuOrbiter_0_2"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuOrbiter_0_2->menuAction());

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0));
        label_3->setText(QApplication::translate("MainWindow", "x0", 0));
        x0LineEdit->setText(QApplication::translate("MainWindow", "-0.49", 0));
        label_4->setText(QApplication::translate("MainWindow", "n", 0));
        dxLineEdit->setText(QApplication::translate("MainWindow", "64", 0));
        label_5->setText(QApplication::translate("MainWindow", "xmax", 0));
        xmaxLineEdit->setText(QApplication::translate("MainWindow", "-0.01", 0));
        label_13->setText(QApplication::translate("MainWindow", "dx", 0));
        newdxLineEdit->setText(QString());
        label_15->setText(QApplication::translate("MainWindow", "\320\236\321\201\320\276\320\261\321\213\320\265 \320\276\321\200\320\261\320\270\321\202\321\213", 0));
        addButton->setText(QApplication::translate("MainWindow", "\320\224\320\276\320\261\320\260\320\262\320\270\321\202\321\214 \320\276\321\200\320\261\320\270\321\202\321\203", 0));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("MainWindow", "Tab 1", 0));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("MainWindow", "Tab 2", 0));
        label->setText(QApplication::translate("MainWindow", "tolerance", 0));
        toleranceLineEdit->setText(QApplication::translate("MainWindow", "1e-11", 0));
        label_11->setText(QApplication::translate("MainWindow", "tmax", 0));
        tmaxLineEdit->setText(QApplication::translate("MainWindow", "10000", 0));
        label_2->setText(QApplication::translate("MainWindow", "dt", 0));
        dtLineEdit->setText(QApplication::translate("MainWindow", "0.01", 0));
        label_16->setText(QApplication::translate("MainWindow", "C0", 0));
        C0LineEdit->setText(QApplication::translate("MainWindow", "4", 0));
        label_6->setText(QApplication::translate("MainWindow", "\320\237\320\260\321\200\320\260\320\274\320\265\321\202\321\200\321\213 \320\263\321\200\320\260\321\204\320\265\321\200\320\260", 0));
        checkBox->setText(QApplication::translate("MainWindow", "\320\220\320\262\321\202\320\276\320\274\320\260\321\202\320\270\321\207\320\265\321\201\320\272\320\270\320\265 \320\263\321\200\320\260\320\275\320\270\321\206\321\213", 0));
        label_7->setText(QApplication::translate("MainWindow", "xmin", 0));
        xminLineEdit->setText(QApplication::translate("MainWindow", "-0.51", 0));
        label_8->setText(QApplication::translate("MainWindow", "xmax", 0));
        xmaxLineEdit_2->setText(QApplication::translate("MainWindow", "0.01", 0));
        label_10->setText(QApplication::translate("MainWindow", "ymin", 0));
        yminLineEdit->setText(QApplication::translate("MainWindow", "0.9", 0));
        label_9->setText(QApplication::translate("MainWindow", "ymax", 0));
        ymaxLineEdit->setText(QApplication::translate("MainWindow", "410", 0));
        calculateButton->setText(QApplication::translate("MainWindow", "\320\235\320\260\321\207\320\260\321\202\321\214 \321\200\320\260\321\201\321\207\320\265\321\202", 0));
        label_12->setText(QApplication::translate("MainWindow", "\320\222\321\213\320\261\320\276\321\200 \320\270\320\275\320\264\320\270\320\272\320\260\321\202\320\276\321\200\320\260 \321\205\320\260\320\276\321\201\320\260", 0));
        radioButton->setText(QApplication::translate("MainWindow", "FLI", 0));
        radioButton_2->setText(QApplication::translate("MainWindow", "OFLI", 0));
        radioButton_3->setText(QApplication::translate("MainWindow", "MEGNO", 0));
        radioButton_4->setText(QApplication::translate("MainWindow", "OMEGNO", 0));
        radioButton_5->setText(QApplication::translate("MainWindow", "SALI", 0));
        label_14->setText(QApplication::translate("MainWindow", "CPU=time: ", 0));
        cputimeLineEdit->setText(QString());
        menuOrbiter_0_2->setTitle(QApplication::translate("MainWindow", "Orbiter 0.2", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
