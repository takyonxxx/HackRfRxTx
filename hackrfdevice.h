#ifndef HACKRFDEVICE_H
#define HACKRFDEVICE_H

#include <QObject>
#include <QDebug>
#include <libhackrf/hackrf.h>
#include "audiooutput.h"

#define hackRfLnaGain                   40
#define hackRfVgaGain                   40
#define AUDIO_SAMPLE_RATE               48*1000
#define hackrfCenterFrequency           100*1000*1000
#define hackrfSampleRate                2000000

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

    bool startHackrf();
    bool stopHackrf();
    bool force_sample_rate( double fs_hz );
    bool ptt() const;
    void setPtt(bool newPtt);

private:
    static int rx_callbackStream(hackrf_transfer* transfer);
    static int tx_callbackStream(hackrf_transfer* transfer);
    std::vector<std::string> listDevices();
    std::vector<std::string> device_serials;
    std::vector<int> device_board_ids;
    hackrf_device* m_device;
    AudioOutput *audioOutput{};
    bool m_ptt;
    int centerFrequency;
    int sampleRate;
};

#endif // HACKRFDEVICE_H
