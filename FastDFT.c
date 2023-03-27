#include <complex.h>
#include <stdlib.h>
#include "FastDFT.h"

#define mcexp cexpf
#define mconj conjf


/* Fast DFT : 创建 */
void fftInit(FastDFT *cx, int N)
{
    int ri = 0, i = 0, (*s)[2];
    mComplex *u, unit;

    u = malloc(sizeof(*u) * N);
    s = malloc(sizeof(ri) * N);

    cx->fftLen = N;
    cx->swp = s, cx->ur1 = u;

    while (++i < N - 1)
    {
        ri ^= N + N / (i ^-i);

        if (i < ri)
        {
            s[0][0] = i;
            s[0][1] = ri, ++s;
        }
    }

    s[0][0] = 0;
    N >>= 1, cx->ur0 = u + N;
    unit = _Complex_I * M_PI/N;

    for (i = 0; i < N; ++i)
    {
        u[i] = mcexp(unit * i);
        u[i + N] = mconj(u[i]);
    }
}


/* Fast DFT : 销毁 */
void fftFree(FastDFT *cx)
{
    free(cx->swp), free(cx->ur1);
}


/* Fast DFT : 执行 */
void fftExec(FastDFT *cx,
             mComplex *arr, char inv)
{
    mComplex cp, *ur; int st, s, p;
    int N = cx->fftLen, (*w)[2];
    ur = inv? cx->ur1 : cx->ur0;

    for (w = cx->swp; w[0][0]; ++w)
    {
        cp = arr[w[0][0]];
        arr[w[0][0]] = arr[w[0][1]];
        arr[w[0][1]] = cp;
    }

    for (st = 2; st <= N; st *= 2)
    {
        int h = st / 2, d = N / st;

        for (s = 0; s < N; s += st)
        {
          for (p = 0; p < h; ++p)
          {
              int u = s+p, v = u+h;
              cp = arr[v] * ur[p*d];
              arr[v] = arr[u] - cp;
              arr[u] = arr[u] + cp;
          }
        }
    }

    if (inv)
    {
        mfloat adj = (mfloat)1 / N;

        for (st = 0; st < N; ++st)
            arr[st] = arr[st] * adj;
    }
}
