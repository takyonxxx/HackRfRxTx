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


int HackRfDevice::tx_callbackStream(hackrf_transfer *transfer) {
    qDebug() << transfer->valid_length;
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->tx_ctx);
    QByteArray data = QByteArray::fromRawData(reinterpret_cast<const char*>(transfer->buffer), transfer->valid_length);
}

int HackRfDevice::rx_callbackStream(hackrf_transfer *transfer)
{
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->rx_ctx);

    // Extract data from HackRF transfer buffer
    const float* rf_data = reinterpret_cast<const float*>(transfer->buffer);
    size_t rf_data_len = transfer->valid_length / sizeof(float);

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
    float cutoff_freq = device->cutoff_freq; // Adjust cutoff frequency as needed
    int num_taps = 101; // Adjust number of taps as needed
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
