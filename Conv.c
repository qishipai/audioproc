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



/* ========== [ FFT 卷积器 ] ========== */

typedef struct {
  int o_size;     /* 卷积段落长度 */
  FastDFT *pfft;
  fComplex *arr;  /* 复数点存储区 */
  fComplex *irf;  /* 频域传递函数 */
  fp32t *ol_buf;  /* 重叠段存储区 */
} FFT_Conv;


void convInit(FFT_Conv *e, FastDFT *fft)
{
    e->o_size = fft->fftLen / 2;
    e->pfft = fft;

    e->ol_buf = malloc(sizeof(fp32t)
                     * e->o_size * 9);

    e->irf = (fComplex *)(
           e->ol_buf + e->o_size * 1);

    e->arr = (fComplex *)(
           e->ol_buf + e->o_size * 5);
}

void convZero(FFT_Conv *e)
{
    int buf_size = e->o_size;

    for (int i = 0; i < buf_size; ++i)
        e->ol_buf[i] = 0.;
}

void convFree(FFT_Conv *e)
{
    free(e->ol_buf), e->ol_buf = NULL;
}

void convSetFIR(FFT_Conv *e, fp32t *irs)
{
    int n = e->pfft->fftLen;
    fComplex *irf = e->irf;

    for (int i = 0; i < e->o_size; ++i)
        irf[i] = fCPLX(irs[i], 0.);

    for (int i = e->o_size; i < n; ++i)
        irf[i].re = irf[i].im = 0.;

    fftExec(e->pfft, irf, 0);
}

void convExec(FFT_Conv *e, fp32t *audio)
{
    int n = e->pfft->fftLen;
    fp32t *olb = e->ol_buf;
    fComplex *arr = e->arr;

    for (int i = 0; i < e->o_size; ++i)
        arr[i] = fCPLX(audio[i], 0.);

    for (int i = e->o_size; i < n; ++i)
        arr[i].re = arr[i].im = 0.;

    fftExec(e->pfft, arr, 0);

    for (int i = 0; i < n; ++i)
        arr[i] = fcpMul(arr[i],
                            e->irf[i]);

    fftExec(e->pfft, arr, 1);

    for (int i = 0; i < e->o_size; ++i)
    {
        audio[i] = olb[i] + arr[i].re;
        olb[i] = arr[i + e->o_size].re;
    }
}









int pos_idx[1551];

struct {
  fp32t pos[1550][3], dat[1550][256][2];
} *HRIR;

static void swp_pos(int a, int b)
{
    int t = pos_idx[a];
    pos_idx[a] = pos_idx[b];
    pos_idx[b] = t;
}

static char cmp_pos(int a, int b)
{
    return (HRIR->pos[pos_idx[a]][1]
         == HRIR->pos[pos_idx[b]][1]
         && HRIR->pos[pos_idx[a]][0]
          < HRIR->pos[pos_idx[b]][0])
         || HRIR->pos[pos_idx[a]][1]
          < HRIR->pos[pos_idx[b]][1];
}

static void sort_pos(int n)
{
    int p = 0, u, v;

    while (++p <= n)
        for (u = p; (v = u / 2) > 0; u = v)
            if (cmp_pos(u, v))
                swp_pos(u, v);

    while (--p > 1)
    {
        swp_pos(1, p);

        for (u = 1; (v = u * 2) < p; u = v)
        {
            if (v + 1 < p && cmp_pos(v + 1, v))
                v += 1;

            if (cmp_pos(v, u))
                swp_pos(u, v);
        }
    }
}




#define CO_LEN 256
#define HR_STP 20

int main()
{
    FILE *hr = fopen("HRTF\\HRIR.floats", "rb");
    HRIR = malloc(sizeof(*HRIR));
    fread(HRIR, sizeof(*HRIR), 1, hr);
    fclose(hr);

    for (int i = 1; i <= 1550; ++i)
        pos_idx[i] = i - 1;

    sort_pos(1550);




    FILE *au_rd = popen("ffmpeg -v error -i \"20230228_153045.m4a\" -ac 1 -ar 48000 -f f32le -", "rb");
    FILE *au_wr = popen("ffmpeg -y -f f32le -ac 2 -ar 48000 -i - 3d.aac", "wb");


    FastDFT FFT;
    fftInit(&FFT, CO_LEN * 2);

    FFT_Conv convL, convR;
    convInit(&convL, &FFT);
    convInit(&convR, &FFT);
    convZero(&convL);
    convZero(&convR);

    static fp32t buf[CO_LEN][2];
    static fp32t aud[2][CO_LEN];

    int t = 0, idx = 1, idx1 = 1;

    while (1)
    {
        fp32t a = (fp32t)(++t % HR_STP) / HR_STP, b = 1 - a;

        if (t % HR_STP == 0)
            idx1 = idx, idx = pos_idx[1 + (t / HR_STP) % 1550];

        for (int i = 0; i < CO_LEN; ++i)
        {
            aud[0][i] = (HRIR->dat[idx][i][0] * a + HRIR->dat[idx1][i][0] * b);
            aud[1][i] = (HRIR->dat[idx][i][1] * a + HRIR->dat[idx1][i][1] * b);
        }

        convSetFIR(&convL, aud[0]);
        convSetFIR(&convR, aud[1]);

        if (fread(aud[0], 4 * CO_LEN, 1, au_rd) != 1)
            break;

        for (int i = 0; i < CO_LEN; ++i)
        {
            aud[0][i] *= 9;
            aud[1][i] = aud[0][i];
        }

        convExec(&convL, aud[0]);
        convExec(&convR, aud[1]);

        for (int i = 0; i < CO_LEN; ++i)
        {
            buf[i][0] = aud[0][i] * 2;
            buf[i][1] = aud[1][i] * 2;
        }

        fwrite(buf, 8 * CO_LEN, 1, au_wr);
    }

    free(HRIR);
    fftFree(&FFT);
    convFree(&convL);
    convFree(&convR);
    return (pclose(au_rd), pclose(au_wr), 0);
}
