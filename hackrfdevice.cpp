#include "hackrfdevice.h"

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

    d_realFftData = new float[DEFAULT_FFT_SIZE];
    d_pwrFftData = new float[DEFAULT_FFT_SIZE]();
    d_iirFftData = new float[DEFAULT_FFT_SIZE];
    for (int i = 0; i < DEFAULT_FFT_SIZE; i++)
        d_iirFftData[i] = RESET_FFT_FACTOR;  // dBFS

    if (hackrf_init() != HACKRF_SUCCESS) {
        throw std::runtime_error("can not init hackrf");
    }

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
    set_sample_rate(hackrfSampleRate);
    set_frequency(hackrfCenterFrequency);

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

void HackRfDevice::performFFT(std::vector<std::complex<float>> &fftData)
{
    int n = fftData.size();
    if (n <= 1)
    {
        // No need to perform FFT on a single element or empty vector.
        return;
    }

    // Helper function to perform the bit reversal permutation
    auto bit_reversal_permutation = [&](int index) {
        int reversed = 0;
        int mask = n >> 1;
        while (mask > 0)
        {
            if (index & mask)
            {
                reversed |= (n >> (mask == n ? 0 : 1));
            }
            mask >>= 1;
        }
        return reversed;
    };

    // Bit reversal permutation
    for (int i = 0; i < n; i++)
    {
        int rev_i = bit_reversal_permutation(i);
        if (i < rev_i)
        {
            std::swap(fftData[i], fftData[rev_i]);
        }
    }

    // FFT computation
    for (int len = 2; len <= n; len *= 2)
    {
        float angle = -2.0f * M_PI / len;
        std::complex<float> wlen = std::polar(1.0f, angle);
        for (int i = 0; i < n; i += len)
        {
            std::complex<float> w = 1.0f;
            for (int j = 0; j < len / 2; j++)
            {
                std::complex<float> u = fftData[i + j];
                std::complex<float> t = w * fftData[i + j + len / 2];

                fftData[i + j] = u + t;
                fftData[i + j + len / 2] = u - t;

                w *= wlen;
            }
        }
    }
}

void HackRfDevice::process_fft(hackrf_transfer *transfer)
{
    if (transfer->valid_length <= 0)
    {
        return;
    }

    auto *data = transfer->buffer;
    int fftsize = static_cast<int>(transfer->buffer_length / 2.0);
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

int HackRfDevice::rx_callbackStream(hackrf_transfer *transfer)
{
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->rx_ctx);   
    device->process_fft(transfer);
    // return 0;

   // Extract data from HackRF transfer buffer
   const float* rf_data = reinterpret_cast<const float*>(transfer->buffer);
   size_t rf_data_len = static_cast<int>(transfer->buffer_length);

   // Downconvert to baseband by mixing with center frequency
   std::vector<float> baseband_data(rf_data_len);
   for (size_t i = 0; i < rf_data_len; ++i) {
       float phase = 2 * M_PI * device->centerFrequency * i / device->sampleRate;
       baseband_data[i] = rf_data[i] * std::cos(phase);
   }

   // Perform FM demodulation
   std::vector<float> demodulated_data(rf_data_len);
   device->fm_demodulation(baseband_data.data(), rf_data_len, demodulated_data.data(), device->sampleRate, device->centerFrequency);

   // Create low-pass filter coefficients
   float cutoff_freq = device->cutoff_freq;
   int num_taps = 5;
   std::vector<float> filter_taps = device->create_lowpass_filter(cutoff_freq, device->sampleRate, num_taps);

   // Apply low-pass filter to demodulated data
   std::vector<float> filtered_data(rf_data_len);
   device->apply_fir_filter(demodulated_data, filter_taps, filtered_data);

   QByteArray data(reinterpret_cast<const char*>(filtered_data.data()), rf_data_len * sizeof(float));

   // Write filtered data to audio output
   if (!device->m_ptt && device->audioOutput && !data.isEmpty()) {
       device->audioOutput->writeBuffer(data);
   }
   return 0;
}
