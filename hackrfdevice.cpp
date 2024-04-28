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
    HackRfDevice *_this = reinterpret_cast<HackRfDevice *>(transfer->rx_ctx);
    _this->process_fft(transfer);

    const float* rf_data = reinterpret_cast<const float*>(transfer->buffer);
    size_t rf_data_len = static_cast<int>(transfer->buffer_length);

    QByteArray data(reinterpret_cast<const char*>(rf_data), rf_data_len * sizeof(float));

    // Write filtered data to audio output
    if (!_this->m_ptt && _this->audioOutput && !data.isEmpty()) {
        _this->audioOutput->writeBuffer(data);
    }
   return 0;
}
