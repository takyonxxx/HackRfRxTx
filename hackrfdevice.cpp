#include "hackrfdevice.h"
#include "constants.h"

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
    audioOutput = new AudioOutput(this, DEFAULT_AUDIO_SAMPLE_RATE);

    d_realFftData = new float[DEFAULT_FFT_SIZE];
    d_pwrFftData = new float[DEFAULT_FFT_SIZE]();
    d_iirFftData = new float[DEFAULT_FFT_SIZE];
    for (int i = 0; i < DEFAULT_FFT_SIZE; i++)
        d_iirFftData[i] = RESET_FFT_FACTOR;  // dBFS

    if (hackrf_init() != HACKRF_SUCCESS) {
        throw std::runtime_error("can not init hackrf");
    }

    sampleRate = DEFAULT_SAMPLE_RATE;
    audioSampleRate = DEFAULT_AUDIO_SAMPLE_RATE;

    listDevices();
}

HackRfDevice::~HackRfDevice()
{
    stopHackrf();
    hackrf_exit();
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

    if (hackrf_set_antenna_enable(m_device, 0) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set antenna");
    }

    set_amp_enabled(false);
    set_lna_gain(HACKRF_RX_LNA_MAX_DB);
    set_vga_gain(HACKRF_RX_VGA_MAX_DB);
    set_sample_rate(DEFAULT_SAMPLE_RATE);
    set_frequency(DEFAULT_FREQUENCY);

    if (m_ptt && hackrf_set_txvga_gain(m_device, HACKRF_TX_VGA_MAX_DB) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp tx");
    }

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

void HackRfDevice::set_frequency(uint64_t freq)
{
    if (hackrf_set_freq(m_device, freq) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set frequency");
    }
    else
        centerFrequency = freq;
}

void HackRfDevice::set_sample_rate(uint64_t srate)
{
    if (!force_sample_rate(srate)) {
        throw std::runtime_error("can not set sample rate");
    }
    else
        sampleRate = srate;
}

void HackRfDevice::set_amp_enabled(bool enabled)
{
    if (hackrf_set_amp_enable(m_device, enabled) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp");
    }
}

void HackRfDevice::set_lna_gain(uint32_t gain)
{
    if (hackrf_set_lna_gain(m_device, gain) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set lna gain");
    }
}

void HackRfDevice::set_vga_gain(uint32_t gain)
{
    if (hackrf_set_vga_gain(m_device, gain) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set vga gain");
    }
}

bool HackRfDevice::is_streaming() const
{
    return ( hackrf_is_streaming(m_device) == HACKRF_TRUE );
}

bool HackRfDevice::force_sample_rate( double fs_hz )
{
    auto baseband_filter_bw_hz = hackrf_compute_baseband_filter_bw_round_down_lt( uint32_t(fs_hz) );
    qDebug() << "HackRFDevice: Setting filter to " << baseband_filter_bw_hz;
    hackrf_set_baseband_filter_bandwidth( m_device, baseband_filter_bw_hz );

    qDebug() << "HackRFDevice: Setting sample rate to " << fs_hz;
    if (hackrf_set_sample_rate(m_device, sampleRate) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set sample rate");
    }
    else
    {
        return true;
    }
    return false;
}

void HackRfDevice::process_fft(hackrf_transfer *transfer)
{
    if (transfer->valid_length <= 0)
    {
        return;
    }

    auto *data = transfer->buffer;
    int fftsize = static_cast<int>(transfer->valid_length / 2.0);
    d_fftAvg = static_cast<float>(1.0 - 1.0e-2 * 75);

    if (fftsize > DEFAULT_FFT_SIZE)
        fftsize = DEFAULT_FFT_SIZE;

    auto pwr_scale = static_cast<float>(1.0 / fftsize);
    double fullScalePower = 1.0;
    std::vector<std::complex<float>> fftData(fftsize);

    // Calculate magnitude spectrum and logarithmic amplitude
    for (int i = 0; i < fftsize; i++) {
        if (i < fftsize / 2)
        {
            fftData[i] = data[fftsize / 2 + i];
        }
        else
        {
            fftData[i] = data[i - fftsize / 2];
        }

        // auto magnitude = pwr_scale * (pt.imag() * pt.imag() + pt.real() * pt.real());
        float magnitude = std::norm(fftData[i]) * pwr_scale; // Calculate power
        float level = 100 * log10(magnitude + 1e-20) / fullScalePower;
        d_realFftData[i] = level;
        d_iirFftData[i] += d_fftAvg * (d_realFftData[i] - d_iirFftData[i]);
    }

    emit setNewFttData(d_iirFftData, d_realFftData, fftsize);
}

int HackRfDevice::tx_callbackStream(hackrf_transfer *transfer) {
    qDebug() << transfer->valid_length;
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->tx_ctx);
    QByteArray data = QByteArray::fromRawData(reinterpret_cast<const char*>(transfer->buffer), transfer->valid_length);
}

// Function to resample IQ data
void resampleIQData(const std::vector<std::complex<float>>& iq_data,
                    std::vector<std::complex<float>>& resampled_iq_data,
                    double resample_ratio) {
    // Calculate the size of the resampled IQ data
    size_t resampled_size = static_cast<size_t>(iq_data.size() * resample_ratio);
    resampled_iq_data.resize(resampled_size);

    // Resample the IQ data using linear interpolation
    for (size_t i = 0; i < resampled_size; ++i) {
        // Calculate the original index for the resampled sample
        double original_index = i / resample_ratio;

        // Linear interpolation between lower and upper indices
        size_t lower_index = static_cast<size_t>(original_index);
        size_t upper_index = std::min(lower_index + 1, iq_data.size() - 1);
        double fraction = original_index - lower_index;

        std::complex<float> interpolated_value = iq_data[lower_index] * std::complex<float>(1 - fraction, 0.0f) +
                                                 iq_data[upper_index] * std::complex<float>(fraction, 0.0f);

        resampled_iq_data[i] = interpolated_value;
    }
}

// Function to perform WFM demodulation
void demodulateWFM(const std::vector<std::complex<float>>& resampled_iq_data,
                   std::vector<float>& demodulated_signal) {
    demodulated_signal.reserve(resampled_iq_data.size());
    for (size_t i = 1; i < resampled_iq_data.size(); ++i) {
        // Calculate the phase difference
        std::complex<float> previous_sample = resampled_iq_data[i - 1];
        std::complex<float> current_sample = resampled_iq_data[i];

        float phase_diff = std::arg(current_sample * std::conj(previous_sample));

        // Store the phase difference in the demodulated signal
        demodulated_signal.push_back(phase_diff);
    }
}

void applyLowPassFilter(const std::vector<float>& input_signal, std::vector<float>& output_signal,
                        int sample_rate, float cutoff_frequency, size_t num_taps = 9500) {
    // Calculate the filter coefficients using the Hamming window method
    std::vector<float> filter_coefficients(num_taps);
    float nyquist = 0.5f * sample_rate;
    float normalized_cutoff = cutoff_frequency / nyquist;

    // Calculate the filter coefficients using the sinc function and Hamming window
    size_t mid = num_taps / 2;
    for (size_t n = 0; n < num_taps; ++n) {
        float hamming_window = 0.54 - 0.46 * cos((2 * M_PI * n) / (num_taps - 1));
        float sinc_value = sin(2 * M_PI * normalized_cutoff * (n - mid)) /
                           (M_PI * (n - mid));
        if (n == mid) {
            sinc_value = 2 * normalized_cutoff;
        }
        filter_coefficients[n] = hamming_window * sinc_value;
    }

    // Normalize the filter coefficients
    float sum = std::accumulate(filter_coefficients.begin(), filter_coefficients.end(), 0.0f);
    for (float& coeff : filter_coefficients) {
        coeff /= sum;
    }

    // Apply the filter to the input signal using convolution
    output_signal.resize(input_signal.size());
    for (size_t i = 0; i < input_signal.size(); ++i) {
        float filtered_sample = 0.0f;
        for (size_t j = 0; j < num_taps; ++j) {
            int index = static_cast<int>(i) - static_cast<int>(j) + static_cast<int>(mid);
            if (index >= 0 && index < static_cast<int>(input_signal.size())) {
                filtered_sample += input_signal[index] * filter_coefficients[j];
            }
        }
        output_signal[i] = filtered_sample;
    }
}


int HackRfDevice::rx_callbackStream(hackrf_transfer *transfer) {
    // Retrieve the HackRfDevice context
    HackRfDevice *_this = reinterpret_cast<HackRfDevice *>(transfer->rx_ctx);

    // Process the FFT (if necessary)
    _this->process_fft(transfer);

    // Retrieve the received data as IQ data
    auto rf_data = reinterpret_cast<int8_t*>(transfer->buffer);
    auto rf_data_len = transfer->valid_length;

    // Check for invalid data length
    if (rf_data_len % 2 != 0) {
        return -1; // Invalid data length
    }

    // Convert received data to complex IQ data
    std::vector<std::complex<float>> iq_data;
    iq_data.reserve(rf_data_len / 2);
    for (size_t i = 0; i < rf_data_len; i += 2) {
        float i_sample = rf_data[i] / 128.0f; // Normalize from int8 to float
        float q_sample = rf_data[i + 1] / 128.0f; // Normalize from int8 to float
        iq_data.emplace_back(i_sample, q_sample);
    }

    // Step 1: Resample the IQ data
    std::vector<std::complex<float>> resampled_iq_data;
    double resample_ratio = static_cast<double>(_this->audioSampleRate) / _this->sampleRate;

    // Call resampleIQData function
    resampleIQData(iq_data, resampled_iq_data, resample_ratio);

    // Step 2: Perform WFM demodulation on the resampled IQ data
    std::vector<float> demodulated_signal;
    demodulateWFM(resampled_iq_data, demodulated_signal);

    // Step 3: Apply low-pass filtering to the demodulated signal
    std::vector<float> filtered_signal;
    applyLowPassFilter(demodulated_signal, filtered_signal, _this->audioSampleRate, DEFAULT_CUT_OFF);

    // Convert the filtered signal to int16_t audio data
    std::vector<int16_t> audio_data(filtered_signal.size());
    std::transform(filtered_signal.begin(), filtered_signal.end(), audio_data.begin(),
                   [](float sample) {
                       return static_cast<int16_t>(sample * 32767); // Convert float to int16
                   });

    // Write the processed audio data to the audio output
    if (!_this->m_ptt && _this->audioOutput && !audio_data.empty()) {
        QByteArray data(reinterpret_cast<char*>(audio_data.data()), audio_data.size() * sizeof(int16_t));
        _this->audioOutput->writeBuffer(data);
    }

    return 0;
}
