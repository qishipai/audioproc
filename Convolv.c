#include <complex.h>
#include <stdlib.h>
#include "Convolv.h"


/* FFT 卷积器 : 创建 */
void convInit(Convolv *cx,
                     FastDFT *pfft)
{
    int csize = pfft->fftLen;
}

/* FFT 卷积器 : 重置 */
void convZero(Convolv *cx)
{
}

/* FFT 卷积器 : 销毁 */
void convFree(Convolv *cx)
{
}

/* FFT 卷积器 : 设置IR */
void convSetIR(Convolv *cx,
                     mfloat *imprs)
{
}

/* FFT 卷积器 : 执行 */
void convExec(Convolv *cx,
                     mfloat *audio)
{
}
