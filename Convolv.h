/* ===== [ 卷积滤波器 ] ===== */

#ifndef __Convolv_H__
#define __Convolv_H__

#include "FastDFT.h"

typedef struct {
  int convLen;    /* 卷积长度 */
  mfloat *obuf;   /* 累加缓冲 */
  FastDFT *pfft;
  mComplex *fir;  /* 频域响应 */
  mComplex *arr;  /* 频域存储 */
} Convolv;

/* 卷积滤波器 : 创建 */
void cnvInit(Convolv *cx,
                 FastDFT *pfft);

/* 卷积滤波器 : 重置 */
void cnvZero(Convolv *cx);

/* 卷积滤波器 : 销毁 */
void cnvFree(Convolv *cx);

/* 卷积滤波器 : 设置IR */
void cnvSetIR(Convolv *cx,
                 mfloat *imprs);

/* 卷积滤波器 : 执行 */
void cnvExec(Convolv *cx,
                 mfloat *audio);

#endif  /* { 云中龙++ 2023 } */
