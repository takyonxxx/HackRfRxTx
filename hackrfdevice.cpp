#include "hackrfdevice.h"

std::string removeZerosFromBegging(const std::string &string) {
    uint32_t i = 0;
    while (i < string.length() && string[i] == '0') {
        i++;
    }
    return string.substr(i, string.length() - i);
}

HackRfDevice::HackRfDevice(QObject *parent)
{
    audioOutput = new AudioOutput(this, DEFAULT_SAMPLE_RATE);

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

bool HackRfDevice::startHackrf()
{
    auto serial = device_serials[0];
    auto board_id = device_board_ids[0];
    qDebug() << "Starting HackRf " << serial.c_str() << board_id;

    if (hackrf_open(&m_device) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not open hackrf device");
    }

    if (hackrf_set_amp_enable(m_device, false) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set amp");
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

    if (hackrf_start_rx(m_device, &HackRfDevice::rx_callbackStream, this) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not start stream");
    }

    if (hackrf_set_freq(m_device, DEFAULT_FREQUENCY) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not set frequency");
    }

    return true;
}

bool HackRfDevice::stopHackrf()
{
    if (hackrf_stop_rx(m_device) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not stop rx stream");
    }

    if (hackrf_close(m_device) != HACKRF_SUCCESS) {
        qDebug() << "Can not close hackrf";
        return false;
    }

    m_device = nullptr;

    qDebug() << "HackRf closed...";
    return true;
}

int HackRfDevice::rx_callbackStream(hackrf_transfer *transfer) {
    qDebug() << transfer->valid_length;
    HackRfDevice *device = reinterpret_cast<HackRfDevice *>(transfer->rx_ctx);

    QByteArray data = QByteArray::fromRawData(reinterpret_cast<const char*>(transfer->buffer), transfer->valid_length);
    if (device->audioOutput && !data.isEmpty()) {
        device->audioOutput->writeBuffer(data);
    }
    return 0;
}
