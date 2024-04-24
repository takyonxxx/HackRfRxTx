#include "hackrfdevice.h"

std::string removeZerosFromBegging(const std::string &string) {
    uint32_t i = 0;
    while (i < string.length() && string[i] == '0') {
        i++;
    }
    return string.substr(i, string.length() - i);
}

HackRfDevice::HackRfDevice(QObject *parent):
    m_ptt(false)
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
    auto board_id = device_board_ids[0];

    if (hackrf_open(&m_device) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not open hackrf device");
    }

    if (hackrf_set_amp_enable(m_device, false) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp rx");
    }
    if (m_ptt && hackrf_set_amp_enable(m_device, true) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp tx");
    }
    if (m_ptt && hackrf_set_txvga_gain(m_device, hackRfVgaGain) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp tx");
    }
    if (hackrf_set_antenna_enable(m_device, 0) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set antenna");
    }  

    if (hackrf_set_lna_gain(m_device, hackRfLnaGain) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set lna gain");
    }
    if (hackrf_set_vga_gain(m_device, hackRfVgaGain) != HACKRF_SUCCESS) {
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


int HackRfDevice::tx_callbackStream(hackrf_transfer *transfer) {
    qDebug() << transfer->valid_length;
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->rx_ctx);
    QByteArray data = QByteArray::fromRawData(reinterpret_cast<const char*>(transfer->buffer), transfer->valid_length);
}

int HackRfDevice::rx_callbackStream(hackrf_transfer *transfer) {
    // qDebug() << transfer->valid_length;
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->rx_ctx);

    QByteArray data = QByteArray::fromRawData(reinterpret_cast<const char*>(transfer->buffer), transfer->valid_length);
    if (device->audioOutput && !data.isEmpty()) {
        device->audioOutput->writeBuffer(data);
    }
    return 0;
}
