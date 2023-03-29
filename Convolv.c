#include <complex.h>
#include <stdlib.h>
#include "Convolv.h"

#define mconj conjf


/* 卷积滤波器 : 创建 */
void cnvInit(Convolv *cx,
                 FastDFT *pfft)
{
    int N = pfft->fftLen / 2;

    mfloat *m = malloc(N * 8
                * sizeof(*m));

    cx->obuf = m;
    cx->pfft = pfft;
    cx->convLen = N;
    cx->fir = (void *)(m + N);
    cx->arr = cx->fir + N + 1;
}

/* 卷积滤波器 : 重置 */
void cnvZero(Convolv *cx)
{
    int i, clen = cx->convLen;

    for (i = 0; i < clen; ++i)
        cx->obuf[i] = .0;
}

/* 卷积滤波器 : 销毁 */
void cnvFree(Convolv *cx)
{
    free(cx->obuf);
}

/* 卷积滤波器 : 设置IR */
void cnvSetIR(Convolv *cx,
                 mfloat *imprs)
{
    int i, clen = cx->convLen;
    mComplex *fir = cx->fir;

    for (i = 0; i < clen; ++i)
    {
        fir[i] = imprs[i];
        fir[i + clen] = 0;
    }

    fftExec(cx->pfft, fir, 0);
}

/* 卷积滤波器 : 执行 */
void cnvExec(Convolv *cx,
                 mfloat *audio)
{
    int i, clen = cx->convLen;
    int fn = cx->pfft->fftLen;
    mComplex *fir = cx->fir;
    mComplex *arr = cx->arr;

    for (i = 0; i < clen; ++i)
    {
        arr[i] = audio[i];
        arr[i + clen] = 0;
    }

    fftExec(cx->pfft, arr, 0);

    arr[0x00] *= fir[0x00];
    arr[clen] *= fir[clen];

    for (i = 1; i < clen; ++i)
    {
        arr[i] *= fir[i];

        arr[fn - i]
              = mconj(arr[i]);
    }

    fftExec(cx->pfft, arr, 1);

    mfloat *obuf = cx->obuf;

    for (i = 0; i < clen; ++i)
    {
        audio[i]
           = obuf[i] + arr[i];

        obuf[i] = arr[i+clen];
    }
}
