#include <complex.h>
#include <stdlib.h>
#include "Convolv.h"

#define mconj conjf


/* 卷积滤波器 : 创建 */
void convInit(Convolv *cx,
                     FastDFT *pfft)
{
    int cnlen = pfft->fftLen / 2;

    mfloat *m = malloc(cnlen * 8
                * sizeof(mfloat));

    cx->convLen = cnlen;
    cx->pfft = pfft, cx->obuf = m;
    cx->fir = (void *)(m + cnlen);
    cx->arr = cx->fir + cnlen * 3;
}

/* 卷积滤波器 : 重置 */
void convZero(Convolv *cx)
{
    int clen = cx->convLen;

    for (int i = 0; i < clen; ++i)
        cx->obuf[i] = .0;
}

/* 卷积滤波器 : 销毁 */
void convFree(Convolv *cx)
{
    free(cx->obuf), cx->convLen = 0;
}

/* 卷积滤波器 : 设置IR */
void convSetIR(Convolv *cx,
                     mfloat *imprs)
{
    int clen = cx->convLen;

    for (int i = 0; i < clen; ++i)
    {
        cx->fir[i] = imprs[i];
        cx->fir[i + clen] = 0;
    }

    fftExec(cx->pfft, cx->fir, 0);
}

/* 卷积滤波器 : 执行 */
void convExec(Convolv *cx,
                     mfloat *audio)
{
    int fftlen = cx->pfft->fftLen;
    mfloat *obuf = cx->obuf;
    int clen = cx->convLen;

    for (int i = 0; i < clen; ++i)
    {
        cx->arr[i] = audio[i];
        cx->arr[i + clen] = 0;
    }

    fftExec(cx->pfft, cx->arr, 0);

    cx->arr[0] *= cx->fir[0];

    for (int i = 1; i < clen; ++i)
    {
        cx->arr[i] *= cx->fir[i];

        cx->arr[fftlen - i]
              = mconj(cx->arr[i]);
    }

    cx->arr[clen]*=cx->fir[clen];

    fftExec(cx->pfft, cx->arr, 1);

    for (int i = 0; i < clen; ++i)
    {
        audio[i]
           = obuf[i] + cx->arr[i];

        obuf[i] = cx->arr[i+clen];
    }
}
