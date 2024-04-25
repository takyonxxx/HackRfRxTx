#ifndef HACKRFDEVICE_H
#define HACKRFDEVICE_H

#include <QObject>
#include <QDebug>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <libhackrf/hackrf.h>
#include "audiooutput.h"

#define _GHZ(x) ((uint64_t)(x) * 1000000000)
#define _MHZ(x) ((x) * 1000000)
#define _KHZ(x) ((x) * 1000)
#define _HZ(x) ((x) * 1)

#define hackrfCenterFrequency           _MHZ(100)
#define hackrfSampleRate                _MHZ(20)
#define AUDIO_SAMPLE_RATE               _KHZ(48)
#define DEFAULT_CUT_OFF                 _KHZ(300)
#define HACKRF_TX_VGA_MAX_DB            47.0
#define HACKRF_RX_VGA_MAX_DB            40.0
#define HACKRF_RX_LNA_MAX_DB            40.0
#define HACKRF_AMP_MAX_DB               14.0
#define DEFAULT_FFT_SIZE                8192 * 4
#define DEFAULT_FFT_RATE                25 //Hz
#define DEFAULT_AUDIO_GAIN              50
#define MAX_FFT_SIZE                 DEFAULT_FFT_SIZE
#define RESET_FFT_FACTOR             -72

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

signals:
    void setNewFttData(float *, float *, int);

private:
    static int rx_callbackStream(hackrf_transfer* transfer);
    static int tx_callbackStream(hackrf_transfer* transfer);
    void fft(double* real, double* imag, int n);
    double* process_fft(hackrf_transfer *transfer, int& length);
    void fm_demodulation(const float* input, size_t input_len, float* output, float sample_rate, float center_freq);
    std::vector<float> create_lowpass_filter(float cutoff_freq, float sample_rate, int num_taps);
    void apply_fir_filter(const std::vector<float>& input, const std::vector<float>& taps, std::vector<float>& output);
    std::vector<std::string> listDevices();
    std::vector<std::string> device_serials;
    std::vector<int> device_board_ids;
    hackrf_device* m_device;
    AudioOutput *audioOutput{};
    bool m_ptt;
    int centerFrequency;
    int sampleRate;
    int cutoff_freq;

    float *d_realFftData;
    float *d_iirFftData;
    float *d_pwrFftData;
    float  d_fftAvg;
};

#endif // HACKRFDEVICE_H
