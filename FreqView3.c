#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif



/* ========== [ 复数运算 ] ========== */

typedef struct {
  float re, im;
} fComplex;

#define fCPLX(m,n) ((fComplex){m, n})

#define fcp_I fCPLX(0.0f, 1.0f)


fComplex fcpfMul(fComplex x, float y)
{
    return fCPLX(x.re * y, x.im * y);
}

fComplex fcpConj(fComplex __cx)
{
    return fCPLX(__cx.re, - __cx.im);
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

fComplex fcpUnit(float __a)
{
    return fCPLX(cosf(__a), sinf(__a));
}



/* ========== [ Fast DFT ] ========== */

typedef struct {
  int fftLen;
  int (*swp)[2];
  fComplex *urt;
} FastDFTeng;


void fftInit(FastDFTeng *e, int N)
{
    int r = 0, (*sw)[2];
    fComplex *ur, *ur_i;
    float ang = -M_PI;

    ur = malloc(sizeof(*ur) * N);
    sw = malloc(sizeof(int) * N);

    e->fftLen = N;
    e->swp = sw, e->urt = ur;

    for (int i = 1; i < N - 1; ++i)
    {
        r ^= N + N / (i ^ -i);

        if (i < r)
        {
            (*sw)[0] = i;
            (*sw)[1] = r, ++sw;
        }
    }

    (*sw)[0] = 0, r = N / 2;
    ang /= r, ur_i = ur + r;

    for (int i = 0; i < r; ++i)
    {
        ur[i] = fcpUnit(ang * i);
        ur_i[i] = fcpConj(ur[i]);
    }
}

void fftFree(FastDFTeng *e)
{
    free(e->urt), free(e->swp);
}

void fftExec(FastDFTeng *e, fComplex *arr, char m)
{
    int N = e->fftLen, (*sp)[2] = e->swp;
    fComplex *ur = e->urt + (m? N / 2 : 0);
    fComplex cx, cy;

    while ((*sp)[0])
    {
        fComplex _tmp = arr[(*sp)[0]];
        arr[(*sp)[0]] = arr[(*sp)[1]];
        arr[(*sp)[1]] = _tmp, ++sp;
    }

    for (int c = 2; c <= N; c *= 2)
    {
        int h = c / 2, d = N / c;

        for (int s = 0; s < N; s += c)
            for (int p = 0; p < h; ++p)
            {
                int u = s | p, v = u | h;
                cx = arr[u], cy = arr[v];
                cy = fcpMul(cy, ur[p * d]);
                arr[u] = fcpAdd(cx, cy);
                arr[v] = fcpSub(cx, cy);
            }
    }

    if (m)
    {
        float invC = (float) 1 / N;

        for (int i = 0; i < N; ++i)
            arr[i] = fcpfMul(arr[i], invC);
    }
}



/* ========== [ 多项式插值 ] ========== */

typedef float (*itp_efn_t)(float x,
               float const *po, int N);

typedef struct {
  fComplex *t0;
  float *po, *t1;
  FastDFTeng FFT;
} PolyITPeng;


float polyCalc(float x,
               float const *po, int N)
{
    float resu = po[0], px = 1;

    for (int i = 1; i < N; ++i)
        resu += po[i] * (px *= x);

    return resu;
}

float polyDeri(float x,
               float const *po, int N)
{
    float resu = po[1], px = 1;

    for (int i = 2; i < N; ++i)
        resu += po[i] * (px *= x) * i;

    return resu;
}


void itpInit(PolyITPeng *e, int N)
{
    float *m;

    m = malloc(sizeof(*m) * N * 5);

    fftInit(&e->FFT, N);
    e->po = m, e->t1 = m + N;
    e->t0 = (fComplex*)(m + N * 3);
}

void itpFree(PolyITPeng *e)
{
    fftFree(&e->FFT), free(e->po);
}

void itpSetClean(PolyITPeng *e)
{
    int N = e->FFT.fftLen * 5;

    for (int i = 0; i < N; ++i)
        e->po[i] = 0;

    e->t0[0].re = e->t1[0] = 1;
    fftExec(&e->FFT, e->t0, 0);
}

void itpAddSampl(PolyITPeng *e,
 float x, float y, itp_efn_t eval)
{
    fComplex *t1 = (void *)e->t1;
    int N = e->FFT.fftLen;

    float y0 = eval(x, e->po, N);
    float y1 = eval(x, e->t1, N);
    float K  = (y - y0) / y1;

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

float itpCalcPly(PolyITPeng *e,
                          float x)
{
    return polyCalc(x,
            e->po, e->FFT.fftLen);
}




#define FFT_LEN 8192
#define SCR_W   1024
#define SCR_H   192
#define COL_N   120

typedef struct {
  unsigned char r, g, b;
} PIX_t;

PIX_t SCR[SCR_W * SCR_H];


static void clearScr()
{
    for (int i = 0; i < SCR_W * SCR_H; ++i)
        SCR[i] = (PIX_t){10, 10, 10};
}

static PIX_t genPIX(int key)
{
    PIX_t p = {0, 0, 0};

    switch((key %= 1530) / 255)
    {
    case 0: p.r = 255, p.g = key;        break;
    case 1: p.g = 255, p.r = 510 - key;  break;
    case 2: p.g = 255, p.b = key - 510;  break;
    case 3: p.b = 255, p.g = 1020 - key; break;
    case 4: p.b = 255, p.r = key - 1020; break;
    case 5: p.r = 255, p.b = 1530 - key; break;
    }

    return p;
}

static void DrawCir(int px, int py, int r, PIX_t cl)
{
    int sqr = r * r, p = py * SCR_W + px;

    for (int x = 0; x <= r; ++x)
        for (int y = 0; y <= r; ++y)
            if (x + y <= r || x * x + y * y <= sqr)
            {
                int d = y * SCR_W;
                SCR[p + d + x] = SCR[p + d - x] = cl;
                SCR[p - d + x] = SCR[p - d - x] = cl;
            }
}

static void Draw(int x0, int y0, PIX_t col0,
                         int x1, int y1, PIX_t col1)
{
    if (y0 < 0)
        y0 = 0;

    if (y1 < 0)
        y1 = 0;

    y1 = SCR_H - y1 - 5, y0 = SCR_H - y0 - 5;

    for (int y = SCR_H - 1; y > y1; --y)
    {
        SCR[y * SCR_W + x1 - 1] = col1;
        SCR[y * SCR_W + x1 + 1] = col1;
    }

    float dy = y1 - y0, dxi = 1.0f / (x1 - x0);

    for (int x = x0; x < x1; ++x)
    {
        float p = dxi * (x - x0);
        int y = y0 + p * dy;

        PIX_t c = (PIX_t){
            col1.r * p + col0.r * (1 - p),
            col1.g * p + col0.g * (1 - p),
            col1.b * p + col0.b * (1 - p),
        };

        SCR[y * SCR_W + x] = c;
        SCR[(y - 1) * SCR_W + x] = c;
    }

    DrawCir(x1, y1, 4, col1);
}

static inline float winFunc(int i)
{
    return 1.0f - cosf(M_PI * i * 2 / FFT_LEN);
}

static inline float getFreqVal(fComplex c)
{
    float e = log10f(c.re * c.re + c.im * c.im + 1e-5f);
    return e > 0? e : 0;
}


static const char* strf(const char *fmt, ...)
{
    static char buf[1024];
    va_list va;

    va_start(va, fmt);
    vsnprintf(buf, 1024, fmt, va);
    va_end(va);

    return buf;
}

static const char* basefn(const char *fn)
{
    static char buf[1024];
    int i = 0;

    for ( ; i < 1024 && fn[i]; ++i)
        buf[i] = fn[i];

    while (--i > 0 && fn[i] != '.'){ }

    return (buf[i] = 0, buf);
}

int main(int ac, char **av)
{
    if (ac != 2)
    {
        printf("使用方法：%s [音频文件]\n", *av);
        sleep(2);
        return 1;
    }

    FILE *vid = popen(strf("ffmpeg  -y -hide_banner "
      "-f rawvideo -pix_fmt rgb24 -r 60 -s 1024x192 "
      "-i - -i \"%s\" -pix_fmt yuv420p \"%s.mp4\"",
                      av[1], basefn(av[1])), "wb");

    FILE *aud = popen(strf("ffmpeg -v error -i \"%s\" "
      "-ac 1 -ar 48000 -f f32le - ", av[1]), "rb");

    fComplex *arr = malloc(sizeof(fComplex)  * FFT_LEN);

    float *in_buf = malloc(sizeof(float) *
                              (FFT_LEN + COL_N * 4));

    int *map_freq = malloc(sizeof(int) *
                         (FFT_LEN / 2 + COL_N + 20));

    int *map_cnt= map_freq + FFT_LEN / 2, itpcnt = 0;
    int (*itprng)[2] = (void *)(map_cnt  + COL_N);
    float  (*flt)[4] = (void *)(in_buf + FFT_LEN);

    FastDFTeng FFT;
    PolyITPeng ITP;
    fftInit(&FFT, FFT_LEN), itpInit(&ITP, 16);


    for (int i = 0; i < COL_N; ++i)
        map_cnt[i] = 0;

    for (int i = 3; i < FFT_LEN / 2; ++i)
    {
        map_freq[i] = logf((float)i / 3) * COL_N * 0.138f;
        // map_freq[i] = (float)(i - 3) * COL_N / 4094;
        ++map_cnt[map_freq[i]];
    }

    for (int i = 0; i < COL_N - 1; ++i)
    {
        if (map_cnt[i] && !map_cnt[i + 1])
            itprng[itpcnt][0] = i;

        if (!map_cnt[i] && map_cnt[i + 1])
            itprng[itpcnt++][1] = i + 1;
    }


    int ret;

    do
    {
        for (int i = 0; i < FFT_LEN - 800; ++i)
            in_buf[i] = in_buf[i + 800];

        ret = fread(in_buf + FFT_LEN - 800, 4, 800, aud);

        for (int i = 0; i < FFT_LEN; ++i)
        {
            arr[i].re = in_buf[i] * winFunc(i);
            arr[i].im = 0;
        }

        fftExec(&FFT, arr, 0);

        for (int i = 0; i < COL_N; ++i)
        {
            flt[i][0] = flt[i][1];
            flt[i][1] = flt[i][2];
            flt[i][2] = flt[i][3];
            flt[i][3] = 0;
        }

        for (int i = 3; i < FFT_LEN / 2; ++i)
            flt[map_freq[i]][3] += getFreqVal(arr[i]);

        for (int i = 0; i < COL_N; ++i)
            if (map_cnt[i] > 1)
                    flt[i][3] /= map_cnt[i];

        for (int i = 0; i < itpcnt; ++i)
        {
            int sta = itprng[i][0], d0 = itprng[i][1] - sta;
            int b0 = flt[itprng[i][0]] < flt[itprng[i][1]];
            
            itpSetClean(&ITP);
            itpAddSampl(&ITP, 0, flt[itprng[i][0]][3], polyCalc);
            itpAddSampl(&ITP, 0, 0, polyDeri);

            if (i < itpcnt - 1)
            {
                int b1 = flt[itprng[i + 1][0]] < flt[itprng[i + 1][1]];

                if (b0 == b1 && ++i)
                {
                    int d1 = d0 + itprng[i][1] - itprng[i][0];
                    itpAddSampl(&ITP, (float)d0 / d1, flt[itprng[i][0]][3], polyCalc);
                    d0 = d1;
                }
            }

            itpAddSampl(&ITP, 1, flt[itprng[i][1]][3], polyCalc);
            itpAddSampl(&ITP, 1, 0, polyDeri);

            for (int p = 1; p <= d0; ++p)
                flt[sta + p][3] = itpCalcPly(&ITP, (float)p / d0);
        }

        clearScr();

        int _x = 0, _y = 0;
        PIX_t _col = (PIX_t){0, 0, 0};

        for (int i = 0; i < COL_N; ++i)
        {
            int y = (flt[i][0] + flt[i][1]
                  + flt[i][2] + flt[i][3] - 0.5f) * 3;

            int x = i * 8.3f + 10;
            PIX_t col = genPIX(y > 0? y << 3 : 0);

            Draw(_x, _y, _col, x, y, col);
            _x = x, _y = y, _col = col;
        }

        fwrite(SCR, 3 * SCR_W * SCR_H, 1, vid);

    } while (ret == 800);

    fftFree(&FFT), itpFree(&ITP);
    free(arr), free(in_buf), free(map_freq);
    return (pclose(vid), pclose(aud), 0);
}
