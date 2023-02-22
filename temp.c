#include <windows.h>


HWAVEIN   hWaveIn;
WAVEHDR waveBlock;
WAVEFORMATEX wfmt;


#define FFT_LEN 2048
#define FFT_STP 256

static short abuf[FFT_STP];
static fp32t aud[FFT_STP * 2];
static int aud_p;


int main(int ac, char **av)
{
    FastDFT FFT;
    STFT_Tuner Tuner;
    fftInit(&FFT, FFT_LEN);
    tunrInit(&Tuner, &FFT, FFT_STP, 96000);
    tunrZero(&Tuner);

    wfmt.wFormatTag = WAVE_FORMAT_PCM;
    wfmt.nSamplesPerSec = 48000;
    wfmt.wBitsPerSample = 16;
    wfmt.nChannels   = 1;
    wfmt.nBlockAlign = 2;
    wfmt.nAvgBytesPerSec = wfmt.nSamplesPerSec * wfmt.nBlockAlign;
    wfmt.cbSize = 0;

    waveBlock.lpData = (LPSTR)abuf;
    waveBlock.dwBufferLength = FFT_STP * 2;
    waveBlock.dwBytesRecorded = 0;
    waveBlock.dwUser = 0;
    waveBlock.dwFlags = WHDR_BEGINLOOP | WHDR_ENDLOOP;
    waveBlock.dwLoops = 1;

    HANDLE w_event = CreateEvent(NULL, 0, 0, NULL);
    waveInOpen(&hWaveIn, WAVE_MAPPER, &wfmt, (DWORD_PTR)w_event, 0L, CALLBACK_EVENT);
    waveInStart(hWaveIn);

    while (1)
    {
        WaitForSingleObject(w_event, 1000);

        if (waveBlock.dwBytesRecorded > 0)
        {
            int n = waveBlock.dwBytesRecorded;
            short *p = abuf;

            while (n--)
                aud[aud_p++] = (float)*p++ / 32767.0;

            if (aud_p >= FFT_STP)
            {
                if (tunrExec(&Tuner, aud, tune_func, aud, NULL))
                {
                    puts("ablock");
                }

                aud_p -= FFT_STP;

                for (int i = 0; i < aud_p; ++i)
                    aud[i] = aud[i + FFT_STP];
            }
        }

        waveInPrepareHeader(hWaveIn, &waveBlock, sizeof(WAVEHDR));
        waveInAddBuffer(hWaveIn, &waveBlock, sizeof(WAVEHDR));
    }

    waveInReset(hWaveIn);
    waveInClose(hWaveIn);
    fftFree(&FFT), tunrFree(&Tuner);
    return 0;
}
