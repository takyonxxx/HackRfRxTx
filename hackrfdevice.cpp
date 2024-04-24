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
    listDevices();
}

HackRfDevice::~HackRfDevice()
{
    stopHackrf();
}

std::vector<std::string> HackRfDevice::listDevices()
{
    std::vector<std::string> serials;
    auto list = hackrf_device_list();
    if (!list) {
        throw std::runtime_error("can not read hackrf devices list");
    }
    for (int i = 0; i < list->devicecount; ++i) {
        if (!list->serial_numbers[i]) {
            throw std::runtime_error("can not read hackrf serial");
        }
        serials.push_back(removeZerosFromBegging(list->serial_numbers[i]));
    }
    hackrf_device_list_free(list);
    return serials;
}

bool HackRfDevice::startHackrf()
{
    if (hackrf_open(&m_device) != HACKRF_SUCCESS) {
        throw std::runtime_error("can not open hackrf device");
    }
    if (hackrf_set_amp_enable(m_device, 0) != HACKRF_SUCCESS) {
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

    qDebug() << "HackRf started...";
    return true;
}

bool HackRfDevice::stopHackrf()
{
    if (hackrf_exit() != HACKRF_SUCCESS) {
        qDebug() << "Can not exit hackrf";
        return false;
    }

    qDebug() << "HackRf stopped...";
    return true;
}
