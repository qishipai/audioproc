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



/* ========== [ 多项式插值 ] ========== */

typedef fp32t (*itp_efn_t)(
  fp32t x, fp32t const *po, int N);

typedef struct {
  FastDFT  FFT;
  fComplex *t0;
  fp32t *po, *t1;
} PolyITP;


void itpInit(PolyITP *e, int N)
{
    fp32t *m;

    m = malloc(sizeof(*m) * N * 5);

    fftInit(&e->FFT, N);
    e->po = m, e->t1 = m + N;
    e->t0 = (fComplex*)(m + N * 3);
}

void itpFree(PolyITP *e)
{
    fftFree(&e->FFT), free(e->po);
}

void itpZero(PolyITP *e)
{
    int N = e->FFT.fftLen * 5;

    for (int i = 0; i < N; ++i)
        e->po[i] = 0;

    e->t0[0].re = e->t1[0] = 1;
    fftExec(&e->FFT, e->t0, 0);
}

void itpAddSampl(PolyITP *e,
 itp_efn_t eval, fp32t x, fp32t y)
{
    fComplex *t1 = (void *)e->t1;
    int N = e->FFT.fftLen;

    fp32t y0 = eval(x, e->po, N);
    fp32t y1 = eval(x, e->t1, N);
    fp32t K  = (y - y0) / y1;

    for (int i = 0; i < N; ++i)
    {
        e->po[i] += K * e->t1[i];

        e->t1[i] =
                e->t1[i + N] = 0;
    }

    t1[0].re = - x, t1[1].re = 1;
    fftExec(&e->FFT, t1, 0);

    for (int i = 0; i < N; ++i)
        e->t0[i] = t1[i] =
         fcpMul(t1[i], e->t0[i]);

    fftExec(&e->FFT, t1, 1);

    for (int i = 1; i < N; ++i)
        e->t1[i] = e->t1[i << 1];
}

fp32t itpCalcPly(PolyITP *e,
        itp_efn_t calc, fp32t _x)
{
    return calc(_x,
            e->po, e->FFT.fftLen);
}


fp32t polyCalc(fp32t x,
               fp32t const *po, int N)
{
    fp32t resu = po[0], px = 1;

    for (int i = 1; i < N; ++i)
        resu += po[i] * (px *= x);

    return resu;
}

fp32t polyDeri(fp32t x,
               fp32t const *po, int N)
{
    fp32t resu = po[1], px = 1;

    for (int i = 2; i < N; ++i)
        resu += po[i] * (px *= x) * i;

    return resu;
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
    PolyITP ITP;
    itpInit(&ITP, 16);
    itpZero(&ITP);

    itpAddSampl(&ITP, polyCalc, 0.1, 0.1);
    itpAddSampl(&ITP, polyCalc, 0.3, 0.5);
    itpAddSampl(&ITP, polyCalc, 0.7, 0.5);
    itpAddSampl(&ITP, polyCalc, 0.9, 0.9);

    itpAddSampl(&ITP, polyDeri, 0.1, 0);
    itpAddSampl(&ITP, polyDeri, 0.9, 0);

    for (int i = 0; i < 16; ++i)
    {
        printf("%12f\n", ITP.t1[i]);
    }

    for (int i = 0; i <= 100; ++i)
    {
        int t = itpCalcPly(&ITP, polyCalc, (fp32t)i / 100) * 60 + 1;
        println(i, t);
    }

    itpFree(&ITP);
    return 0;
}
