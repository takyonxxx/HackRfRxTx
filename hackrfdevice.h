#ifndef HACKRFDEVICE_H
#define HACKRFDEVICE_H

#include <QObject>
#include <QDebug>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <complex>
#include <libhackrf/hackrf.h>
#include "audiooutput.h"

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

    void set_frequency( uint64_t freq );
    void set_sample_rate( uint64_t srate );
    void set_amp_enabled( bool enabled );
    void set_lna_gain( uint32_t gain );
    void set_vga_gain( uint32_t gain );
    bool is_streaming() const;

    bool force_sample_rate( double fs_hz );
    bool ptt() const;
    void setPtt(bool newPtt);

signals:
    void setNewFttData(float *, float *, int);

private:
    static int rx_callbackStream(hackrf_transfer* transfer);
    static int tx_callbackStream(hackrf_transfer* transfer);

    void process_fft(hackrf_transfer *transfer);    
    std::vector<std::string> listDevices();
    std::vector<std::string> device_serials;
    std::vector<int> device_board_ids;
    hackrf_device* m_device;
    AudioOutput *audioOutput{};

    bool m_ptt;
    int centerFrequency;
    int sampleRate;
    int audioSampleRate;
    int cutoff_freq;

    float *d_realFftData;
    float *d_iirFftData;
    float *d_pwrFftData;
    float  d_fftAvg;


};
#endif // HACKRFDEVICE_H
