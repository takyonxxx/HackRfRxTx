#ifndef CONSTANTS_H
#define CONSTANTS_H

#define _GHZ(x) ((uint64_t)(x) * 1000000000)
#define _MHZ(x) ((x) * 1000000)
#define _KHZ(x) ((x) * 1000)
#define _HZ(x) ((x) * 1)

#define DEFAULT_FREQUENCY              _MHZ(100)
#define DEFAULT_SAMPLE_RATE            _MHZ(20)
#define DEFAULT_AUDIO_SAMPLE_RATE      _KHZ(48)
#define DEFAULT_CUT_OFF                _KHZ(300)
#define HACKRF_TX_VGA_MAX_DB            47.0
#define HACKRF_RX_VGA_MAX_DB            40.0
#define HACKRF_RX_LNA_MAX_DB            40.0
#define HACKRF_AMP_MAX_DB               14.0
#define DEFAULT_FFT_SIZE                DEFAULT_SAMPLE_RATE
#define RESET_FFT_FACTOR                0

#define DB_M_PI     3.14159265358979323846
#define FL_M_PI     3.1415926535f

#define DB_M_SQRT2  1.4142135623730951
#define FL_M_SQRT2  1.4142135623f

#endif // CONSTANTS_H
