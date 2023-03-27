/* ===== [ FFT 卷积器 ] ===== */

#ifndef __Convolv_H__
#define __Convolv_H__

#include "FastDFT.h"

typedef struct {
  mfloat *obuf;   /* 累加缓冲 */
  FastDFT *pfft;
  mComplex *fir;  /* 频域响应 */
  mComplex *arr;  /* 频域存储 */
} Convolv;

/* FFT 卷积器 : 创建 */
void convInit(Convolv *cx,
                 FastDFT *pfft);

/* FFT 卷积器 : 重置 */
void convZero(Convolv *cx);

/* FFT 卷积器 : 销毁 */
void convFree(Convolv *cx);

/* FFT 卷积器 : 设置IR */
void convSetIR(Convolv *cx,
                 mfloat *imprs);

/* FFT 卷积器 : 执行 */
void convExec(Convolv *cx,
                 mfloat *audio);

#endif  /* { 云中龙++ 2023 } */