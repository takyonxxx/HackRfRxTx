#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "constants.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle("HackRf");
    hackRfDevice = new HackRfDevice();
    connect(hackRfDevice, &HackRfDevice::setNewFttData, this, &MainWindow::getNewFttData);

    ui->pushToogleHackrf->setStyleSheet("font-size: 24pt; font: bold; color: #ffffff; background-color: #097532;");
    ui->pushExit->setStyleSheet("font-size: 24pt; font: bold; color: #ffffff; background-color: #900C3F;");

    ui->plotter->setTooltipsEnabled(true);

    ui->plotter->setSampleRate(DEFAULT_SAMPLE_RATE);
    ui->plotter->setSpanFreq(static_cast<quint32>(DEFAULT_SAMPLE_RATE));
    ui->plotter->setCenterFreq(static_cast<quint64>(DEFAULT_FREQUENCY));

    ui->plotter->setFftRange(-140.0f, 20.0f);
    ui->plotter->setPandapterRange(-140.f, 20.f);
    auto m_LowCutFreq = -1* DEFAULT_CUT_OFF;
    auto m_HiCutFreq = DEFAULT_CUT_OFF;
    ui->plotter->setHiLowCutFrequencies(m_LowCutFreq, m_HiCutFreq);
    ui->plotter->setDemodRanges(m_LowCutFreq, -5000, 5000,m_HiCutFreq, true);

    ui->plotter->setFreqUnits(1000);
    ui->plotter->setPercent2DScreen(50);
    ui->plotter->setFilterBoxEnabled(true);
    ui->plotter->setCenterLineEnabled(true);
    ui->plotter->setClickResolution(1);

    ui->plotter->setFftPlotColor(QColor("#CEECF5"));
    ui->plotter->setFreqStep(50000);

    //ui->plotter->setPeakDetection(true ,2);
    ui->plotter->setFftFill(true);
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


void MainWindow::on_radioPtt_clicked(bool checked)
{
    if(checked)
        hackRfDevice->setPtt(true);
    else
        hackRfDevice->setPtt(false);
}

void MainWindow::getNewFttData(float *d_iirFftData, float *d_realFftData, int fftsize)
{
    ui->plotter->setNewFttData(d_iirFftData, d_realFftData, static_cast<int>(fftsize));
}

