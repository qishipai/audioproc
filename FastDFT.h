/* ====== [ Fast DFT ] ====== */

#ifndef __FastDFT_H__
#define __FastDFT_H__

#ifndef M_PI
#define M_PI 3.1415926535897932L
#endif

typedef float _Complex mComplex;
typedef float mfloat;

typedef struct {
  int fftLen;     /* 变换点数 */
  int (*swp)[2];  /* 蝶形索引 */
  mComplex *ur0;  /* 顺单位根 */
  mComplex *ur1;  /* 逆单位根 */
} FastDFT;

/* Fast DFT : 创建 */
void fftInit(FastDFT *cx, int N);

/* Fast DFT : 销毁 */
void fftFree(FastDFT *cx);

/* Fast DFT : 执行 */
void fftExec(FastDFT *cx,
        mComplex *arr, char inv);

#endif  /* ={ 云中龙++ 2023 }= */
