#ifndef AUDIOOUTPUT_H
#define AUDIOOUTPUT_H

#include <QObject>
#include <QThread>
#include <QMutex>
#include <QQueue>
#include <QBuffer>
#include <QAudioSink>
#include <QMediaDevices>
#include <QAudioFormat>
#include <QAudioOutput>
#include <QIODevice>


class AudioOutput: public QObject
{
    Q_OBJECT
public:
    explicit AudioOutput(QObject *parent, int sampleFormat);
    ~AudioOutput();
    void stop();
    void writeBuffer(const QByteArray &buffer);

private slots:
    void handleAudioOutputStateChanged(QAudio::State newState);
    void generateSineWave(short *buffer, int bufferSize, int sampleRate, int frequency);

private:

    QMutex *mutex{};
    bool m_abort {false};

    QAudioFormat m_format;
    QScopedPointer<QAudioSink> m_audioOutput;
    QIODevice *audioDevice;
};

#endif // AUDIOOUTPUT_H
