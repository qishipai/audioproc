#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846264338328L
#endif 

typedef float fp32t;



/* ========== [ 复数运算 ] ========== */

typedef struct {
  fp32t re, im;  /* [实部], [虚部] */
} fComplex;

#define fCPLX(m, n)  ((fComplex){m, n})

#define fcp_I  ((fComplex){0.0, 1.0})


fComplex fcpfMul(fComplex x, fp32t y)
{
    return fCPLX(x.re * y, x.im * y);
}

fComplex fcpUnit(fp32t _a)
{
    return fCPLX(cosf(_a), sinf(_a));
}

fComplex fcpConj(fComplex __cp)
{
    return fCPLX(__cp.re, - __cp.im);
}

fComplex fcpAdd(fComplex x, fComplex y)
{
    return fCPLX(
             x.re + y.re, x.im + y.im);
}

fComplex fcpSub(fComplex x, fComplex y)
{
    return fCPLX(
             x.re - y.re, x.im - y.im);
}

fComplex fcpMul(fComplex x, fComplex y)
{
    return fCPLX(
            x.re * y.re - x.im * y.im,
            x.re * y.im + x.im * y.re);
}



/* ========== [ Fast DFT ] ========== */

typedef struct {
  int fftLen;     /* 单次变换点数 */
  int (*swp)[2];  /* 交换位置列表 */
  fComplex *urt;  /* 单位根存储区 */
} FastDFT;


void fftInit(FastDFT *e, int samplN)
{
    int N = samplN, r = 0, i = 0;
    int (*sw)[2]; fp32t ang;
    fComplex *ur, *ur_i;

    ur = malloc(sizeof(*ur) * N);
    sw = malloc(sizeof(int) * N);

    e->fftLen = samplN;
    e->swp = sw, e->urt = ur;

    while (++i < N - 1)
    {
        r ^= N + N / (i ^ - i);

        if (i < r)
        {
            (*sw)[0] = i;
            (*sw)[1] = r, ++sw;
        }
    }

    (*sw)[0] = 0, N = samplN / 2;
    ang = -M_PI/N, ur_i = ur + N;

    for (i = 0; i < N; ++i)
    {
        ur[i] = fcpUnit(ang * i);
        ur_i[i] = fcpConj(ur[i]);
    }
}

void fftFree(FastDFT *e)
{
    free(e->urt), free(e->swp);
}

void fftExec(FastDFT *e, fComplex *arr, char m)
{
    int N = e->fftLen, (*sptr)[2] = e->swp;
    fComplex *ur = e->urt + (m? N / 2 : 0);
    fComplex cx, cy;

    while ((*sptr)[0])
    {
        fComplex __stmp = arr[(*sptr)[0]];
        arr[(*sptr)[0]] = arr[(*sptr)[1]];
        arr[(*sptr)[1]] = __stmp, ++sptr;
    }

    for (int c = 2; c <= N; c *= 2)
    {
        int h = c / 2, d = N / c;

        for (int s = 0; s < N; s += c)
        {
            for (int p = 0; p < h; ++p)
            {
                int u = s | p, v = u | h;
                cx = arr[u], cy = arr[v];
                cy = fcpMul(cy, ur[p * d]);
                arr[u] = fcpAdd(cx, cy);
                arr[v] = fcpSub(cx, cy);
            }
        }
    }

    if (m)
    {
        fp32t invC = (fp32t) 1 / N;

        for (int i = 0; i < N; ++i)
            arr[i] = fcpfMul(arr[i], invC);
    }
}




void println(int i, int t)
{
    printf("%03d> ", i);

    if (t > 60)
        t = 60;
    else
    if (t < 0)
        t = 0;

    while (t-- > 0)
        putchar(')');

    putchar('\n');
}

int main()
{
    FastDFT FFT;
    fftInit(&FFT, 256);

    fComplex arr[256];

    for (int i = 0; i < 256; ++i)
    {
        char s = abs(i - 40) < 16;
        arr[i].re = s? i % 2 : 0;
        arr[i].im = 0;
    }

    fftExec(&FFT, arr, 0);

    /*for (int i = 0; i < 128; ++i)
    {
        // fComplex t   = arr[i + 128];
        arr[i + 128] = arr[i];
        // arr[i] = t;
    }*/

    fftExec(&FFT, arr, 1);

    for (int i = 0; i < 256; ++i)
    {
        fp32t v = sqrtf(arr[i].re * arr[i].re + arr[i].im * arr[i].im);
        println(i, v * 20 + 0.5f);
    }

    fftFree(&FFT);
    return 0;
}
