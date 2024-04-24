#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <hackrfdevice.h>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushToogleHackrf_clicked();

    void on_pushExit_clicked();

    void on_radioPtt_clicked(bool checked);

private:
    HackRfDevice *hackRfDevice{};

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
