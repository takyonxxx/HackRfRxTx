#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle("HackRf");
    hackRfDevice = new HackRfDevice();
    ui->pushToogleHackrf->setStyleSheet("font-size: 24pt; font: bold; color: #ffffff; background-color: #097532;");
    ui->pushExit->setStyleSheet("font-size: 24pt; font: bold; color: #ffffff; background-color: #900C3F;");
}

MainWindow::~MainWindow()
{
    delete hackRfDevice;
    delete ui;
}


void MainWindow::on_pushToogleHackrf_clicked()
{
    if(ui->pushToogleHackrf->text() == "Start")
    {
        hackRfDevice->startHackrf();
        ui->pushToogleHackrf->setText("Stop");
    }
    else
    {
        hackRfDevice->stopHackrf();
        ui->pushToogleHackrf->setText("Start");
    }
}


void MainWindow::on_pushExit_clicked()
{
    exit(0);
}

