#ifndef HACKRFDEVICE_H
#define HACKRFDEVICE_H

#include <QObject>
#include <QDebug>
#include <libhackrf/hackrf.h>

#define hackRfLnaGain 40
#define hackRfVgaGain 40

enum HackRF_Format {
    HACKRF_FORMAT_FLOAT32	=0,
    HACKRF_FORMAT_INT16	=1,
    HACKRF_FORMAT_INT8	=2,
    HACKRF_FORMAT_FLOAT64 =3,
};

typedef enum {
    HACKRF_TRANSCEIVER_MODE_OFF = 0,
    HACKRF_TRANSCEIVER_MODE_RX = 1,
    HACKRF_TRANSCEIVER_MODE_TX = 2,
} HackRF_transceiver_mode_t;

class HackRfDevice: public QObject
{
    Q_OBJECT
public:
    explicit HackRfDevice(QObject *parent = nullptr);
    ~HackRfDevice();

    static std::vector<std::string> listDevices();
    bool startHackrf();
    bool stopHackrf();

private:
    hackrf_device* m_device;

};

#endif // HACKRFDEVICE_H
