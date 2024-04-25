#include "hackrfdevice.h"
#include <complex>
const double PI = 3.14159265358979323846;

std::string removeZerosFromBegging(const std::string &string) {
    uint32_t i = 0;
    while (i < string.length() && string[i] == '0') {
        i++;
    }
    return string.substr(i, string.length() - i);
}

HackRfDevice::HackRfDevice(QObject *parent):
    m_ptt(false), cutoff_freq(DEFAULT_CUT_OFF)
{
    audioOutput = new AudioOutput(this, AUDIO_SAMPLE_RATE);

    d_realFftData = new float[MAX_FFT_SIZE];
    d_pwrFftData = new float[MAX_FFT_SIZE]();
    d_iirFftData = new float[MAX_FFT_SIZE];
    for (int i = 0; i < MAX_FFT_SIZE; i++)
        d_iirFftData[i] = RESET_FFT_FACTOR;  // dBFS

    if (hackrf_init() != HACKRF_SUCCESS) {
        throw std::runtime_error("can not init hackrf");
    }

    listDevices();
}

HackRfDevice::~HackRfDevice()
{
    stopHackrf();
}

std::vector<std::string> HackRfDevice::listDevices()
{    
    auto list = hackrf_device_list();    
    if (!list) {
        throw std::runtime_error("can not read hackrf devices list");
    }
    for (int i = 0; i < list->devicecount; ++i) {
        if (!list->serial_numbers[i]) {
            throw std::runtime_error("can not read hackrf serial");
        }
        device_serials.push_back(removeZerosFromBegging(list->serial_numbers[i]));
        device_board_ids.push_back(list->usb_board_ids[i]);
        qDebug() << "Found HackRf " << removeZerosFromBegging(list->serial_numbers[i]) << list->usb_board_ids[i];
    }
    hackrf_device_list_free(list);
    return device_serials;
}

bool HackRfDevice::ptt() const
{
    return m_ptt;
}

void HackRfDevice::setPtt(bool newPtt)
{
    m_ptt = newPtt;
}

bool HackRfDevice::startHackrf()
{
    auto serial = device_serials[0];
    // auto board_id = device_board_ids[0];

    if (hackrf_open(&m_device) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not open hackrf device");
    }

    if (hackrf_set_amp_enable(m_device, false) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp rx");
    }
    if (m_ptt && hackrf_set_amp_enable(m_device, true) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp tx");
    }
    if (m_ptt && hackrf_set_txvga_gain(m_device, HACKRF_TX_VGA_MAX_DB) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp tx");
    }
    if (hackrf_set_antenna_enable(m_device, 0) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set antenna");
    }  

    if (hackrf_set_lna_gain(m_device, HACKRF_RX_LNA_MAX_DB) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set lna gain");
    }
    if (hackrf_set_vga_gain(m_device, HACKRF_RX_LNA_MAX_DB) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set vga gain");
    }

    if (!force_sample_rate(hackrfSampleRate)) {
        throw std::runtime_error("can not set sample rate");
    }
    else
        sampleRate = hackrfSampleRate;

    if (hackrf_set_freq(m_device, hackrfCenterFrequency) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set frequency");
    }
    else
        centerFrequency = hackrfCenterFrequency;

    if (m_ptt && hackrf_start_tx(m_device, &HackRfDevice::tx_callbackStream, this) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not start tx stream");
    }
    if (!m_ptt && hackrf_start_rx(m_device, &HackRfDevice::rx_callbackStream, this) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not start rx stream");
    }

    if(m_ptt)
        qDebug() << "Started HackRf Tx Mode" << serial.c_str();
    else
        qDebug() << "Started HackRf Rx Mode" << serial.c_str();


    return true;
}

bool HackRfDevice::stopHackrf()
{
    if (m_ptt && hackrf_stop_tx(m_device) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not stop rx stream");
    }

    if (!m_ptt && hackrf_stop_rx(m_device) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not stop tx stream");
    }

    if (hackrf_close(m_device) != HACKRF_SUCCESS) {
        qDebug() << "Can not close hackrf";
        return false;
    }

    m_device = nullptr;
    qDebug() << "HackRf closed...";
    return true;
}

bool HackRfDevice::force_sample_rate( double fs_hz )
{
    auto baseband_filter_bw_hz = hackrf_compute_baseband_filter_bw_round_down_lt( uint32_t(fs_hz) );
    qDebug() << "HackRFDevice: Setting filter to " << baseband_filter_bw_hz;
    hackrf_set_baseband_filter_bandwidth( m_device, baseband_filter_bw_hz );

    qDebug() << "HackRFDevice: Setting sample rate to " << fs_hz;
    if (hackrf_set_sample_rate(m_device, hackrfSampleRate) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set sample rate");
    }
    else
    {
        return true;
    }
    return false;
}

void HackRfDevice::fm_demodulation(const float* input, size_t input_len, float* output, float sample_rate, float center_freq) {
    float prev_phase = 0.0;

    // Compute instantaneous frequency
    for (size_t i = 0; i < input_len; ++i) {
        float phase = std::atan2(input[i], input[i-1]);
        float deviation = phase - prev_phase;
        // Compute frequency deviation (change in phase over time)
        output[i] = deviation * sample_rate / (2 * M_PI);
        prev_phase = phase;
    }
}


std::vector<float> HackRfDevice::create_lowpass_filter(float cutoff_freq, float sample_rate, int num_taps)
{
    std::vector<float> taps(num_taps);
    float omega_c = 2 * M_PI * cutoff_freq / sample_rate;
    float window_sum = 0.0;
    float n_minus_center = num_taps / 2.0;
    for (int n = 0; n < num_taps; ++n) {
        if (n == n_minus_center) {
            taps[n] = 2 * cutoff_freq / sample_rate;
        } else {
            taps[n] = (std::sin(omega_c * (n - n_minus_center)) / (M_PI * (n - n_minus_center))) * (0.54 - 0.46 * std::cos(2 * M_PI * n / (num_taps - 1)));
        }
        window_sum += taps[n];
    }
    // Normalize taps
    for (int n = 0; n < num_taps; ++n) {
        taps[n] /= window_sum;
    }
    return taps;
}

void HackRfDevice::apply_fir_filter(const std::vector<float>& input, const std::vector<float>& taps, std::vector<float>& output)
{
    int num_taps = taps.size();
    int num_samples = input.size();
    int center_tap = num_taps / 2;

    // Apply FIR filter
    for (int i = 0; i < num_samples; ++i) {
        float acc = 0.0;
        int start_idx = std::max(0, i - center_tap);
        int end_idx = std::min(num_samples - 1, i + center_tap);
        for (int j = start_idx; j <= end_idx; ++j) {
            acc += input[j] * taps[center_tap + (j - i)];
        }
        output[i] = acc;
    }
}

// Function to perform an iterative Cooley-Tukey FFT
void HackRfDevice::fft(double* real, double* imag, int n)
{
    // Bit reversal permutation
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (i < j) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }
        int m = n / 2;
        while (j >= m && m > 0) {
            j -= m;
            m /= 2;
        }
        j += m;
    }

    // FFT computation
    for (int len = 2; len <= n; len *= 2) {
        double angle = -2.0 * PI / len;
        double wlen_cos = cos(angle);
        double wlen_sin = sin(angle);
        for (int i = 0; i < n; i += len) {
            double w_real = 1.0;
            double w_imag = 0.0;
            for (int j = 0; j < len / 2; j++) {
                double u_real = real[i + j];
                double u_imag = imag[i + j];
                double v_real = real[i + j + len / 2];
                double v_imag = imag[i + j + len / 2];

                double t_real = w_real * v_real - w_imag * v_imag;
                double t_imag = w_real * v_imag + w_imag * v_real;

                real[i + j] = u_real + t_real;
                imag[i + j] = u_imag + t_imag;
                real[i + j + len / 2] = u_real - t_real;
                imag[i + j + len / 2] = u_imag - t_imag;

                double w_temp = w_real * wlen_cos - w_imag * wlen_sin;
                w_imag = w_real * wlen_sin + w_imag * wlen_cos;
                w_real = w_temp;
            }
        }
    }
}

double* HackRfDevice::process_fft(hackrf_transfer *transfer, int& length)
{
    uint8_t *data = transfer->buffer;
    length = transfer->buffer_length / 2; // Each I/Q sample consists of 2 bytes
    double* log_amplitude = new double[length];
    d_fftAvg = static_cast<float>(1.0 - 1.0e-2 * 90);
    auto fftsize = static_cast<unsigned int>(length);

    if (fftsize > MAX_FFT_SIZE)
        fftsize = MAX_FFT_SIZE;

    if (fftsize == 0)
    {
        return 0;
    }

    auto pwr_scale = static_cast<float>(1.0 / fftsize);
    double fullScalePower = 1.0;

    // Create arrays for real and imaginary parts
    double *real = new double[fftsize];
    double *imag = new double[fftsize];
//    std::complex<float> pt;


    // Calculate the magnitude and logarithmic amplitude
    for (int i = 0; i < fftsize; i++) {
        if (i < fftsize / 2)
        {
            imag[i]  = data[fftsize / 2 + i];
        }
        else
        {
            imag[i] = data[i - fftsize / 2];
        }

//      fft(real, imag, fftsize);
        double magnitude = pwr_scale * (real[i] * real[i] + imag[i] * imag[i]);
//      double magnitude = pwr_scale * (pt.imag() * pt.imag() + pt.real() * pt.real());

        // Calculate the logarithmic amplitude
        log_amplitude[i] = 10 * log10(magnitude + 1.0e-20f) / fullScalePower;
        d_realFftData[i] = log_amplitude[i];
        d_iirFftData[i] += d_fftAvg * (d_realFftData[i] - d_iirFftData[i]);
    }

    emit setNewFttData(d_iirFftData, d_realFftData, static_cast<int>(fftsize));

    // Clean up the real and imaginary arrays
    delete[] real;
    delete[] imag;

    // Return the logarithmic amplitude array
    return log_amplitude;
}


int HackRfDevice::tx_callbackStream(hackrf_transfer *transfer) {
    qDebug() << transfer->valid_length;
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->tx_ctx);
    QByteArray data = QByteArray::fromRawData(reinterpret_cast<const char*>(transfer->buffer), transfer->valid_length);
}

int HackRfDevice::rx_callbackStream(hackrf_transfer *transfer)
{
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->rx_ctx);

    int length = 0;
    double* amplitude = device->process_fft(transfer, length);
    delete[] amplitude;

//    // Extract data from HackRF transfer buffer
//    const float* rf_data = reinterpret_cast<const float*>(transfer->buffer);
//    size_t rf_data_len = transfer->valid_length / sizeof(float);

//    // Downconvert to baseband by mixing with center frequency
//    std::vector<float> baseband_data(rf_data_len);
//    for (size_t i = 0; i < rf_data_len; ++i) {
//        float phase = 2 * M_PI * device->centerFrequency * i / device->sampleRate;
//        baseband_data[i] = rf_data[i] * std::cos(phase);
//    }

//    // Perform FM demodulation
//    std::vector<float> demodulated_data(rf_data_len);
//    device->fm_demodulation(baseband_data.data(), rf_data_len, demodulated_data.data(), device->sampleRate, device->centerFrequency);

//    // Create low-pass filter coefficients
//    float cutoff_freq = device->cutoff_freq; // Adjust cutoff frequency as needed
//    int num_taps = 101; // Adjust number of taps as needed
//    std::vector<float> filter_taps = device->create_lowpass_filter(cutoff_freq, device->sampleRate, num_taps);

//    // Apply low-pass filter to demodulated data
//    std::vector<float> filtered_data(rf_data_len);
//    device->apply_fir_filter(demodulated_data, filter_taps, filtered_data);


//    QByteArray data(reinterpret_cast<const char*>(filtered_data.data()), rf_data_len * sizeof(float));

//    // Write filtered data to audio output
//    if (!device->m_ptt && device->audioOutput && !data.isEmpty()) {
//        device->audioOutput->writeBuffer(data);
//    }
    return 0;
}
