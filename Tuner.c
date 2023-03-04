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



/* ========== [ STFT 变换 ] ========== */

typedef void (*trans_func_t)(
  fComplex *arr, int fsize, void *usr);

typedef struct {
  int w_step;     /* 单次步进点数 */
  int pbw, pbo;   /* 环形队列指针 */
  fp32t au_amp;   /* 幅度补偿系数 */
  FastDFT *pfft;
  fComplex *arr;  /* 复数点存储区 */
  fp32t *wn_buf;  /* 滑动窗口队列 */
  fp32t *ol_buf;  /* 重叠窗口队列 */
  fp32t *win_fn;  /* 窗函数值列表 */
} STFT_Trans;


void stftInit(STFT_Trans *e,
              FastDFT *fft, int w_step)
{
    int size = fft->fftLen;

    fp32t *mem = malloc(sizeof(fp32t)
                          * size * 5);

    e->arr = (fComplex *) mem;
    e->wn_buf = mem + size * 2;
    e->ol_buf = mem + size * 3;
    e->win_fn = mem + size * 4;

    for (int i = 0; i < size; ++i)
        e->win_fn[i] = 1.0 - cosf(
           M_PI * 2 * i / (size - 1));

    e->w_step = w_step, e->pfft = fft;
    e->au_amp = (fp32t) w_step / size;
}

void stftZero(STFT_Trans *e)
{
    int sz = e->pfft->fftLen;
    e->pbw = 0, e->pbo = 0;

    for (int i = 0; i < sz; ++i)
    {
        e->wn_buf[i] = 0.;
        e->ol_buf[i] = 0.;
    }
}

void stftFree(STFT_Trans *e)
{
    free(e->arr), e->arr = NULL;
}

void stftExec(STFT_Trans *e, fp32t *U,
          trans_func_t tfunc, void *td)
{
    int fft_size = e->pfft->fftLen;
    int fftsize2 = fft_size / 2;
    int pmask = fft_size - 1;
    int step_len = e->w_step;
    fp32t *wn_buf = e->wn_buf;
    fp32t *ol_buf = e->ol_buf;
    fp32t a_amp = e->au_amp;

    for (int i = 0; i < step_len; ++i)
    {
        wn_buf[e->pbw] = U[i];
        e->pbw = (e->pbw + 1) & pmask;
    }

    for (int i = 0; i < fft_size; ++i)
    {
        int pb = (e->pbw + i) & pmask;

        e->arr[i] = fCPLX(wn_buf[pb]
                   * e->win_fn[i], 0);
    }

    fftExec(e->pfft, e->arr, 0);
    tfunc(e->arr, fftsize2, td);

    for (int i = 1; i < fftsize2; ++i)
    {
        e->arr[fft_size - i]
                 = fcpConj(e->arr[i]);
    }

    fftExec(e->pfft, e->arr, 1);

    for (int i = 0; i < fft_size; ++i)
    {
        int pb = (e->pbo + i) & pmask;

        ol_buf[pb] += e->win_fn[i]
                       * e->arr[i].re;
    }

    for (int i = 0; i < step_len; ++i)
    {
        U[i] = ol_buf[e->pbo] * a_amp;

        ol_buf[e->pbo] = 0.;
        e->pbo = (e->pbo + 1) & pmask;
    }
}



/* ========== [ STFT 移调机 ] ========== */

typedef struct {
  int step;       /* 单次处理点数 */
  int sampr;      /* 源音频采样率 */
  int fifo_cnt;   /* 缓冲区计数器 */
  void *tu_func;  /* 频比计算函数 */
  void *tu_data;  /* 频比用户数据 */
  STFT_Trans ft;
  fp32t *pcache;  /* 相位值暂存区 */
  fp32t *pcount;  /* 相位值累加区 */
  fp32t *m_freq;  /* 频率点测量值 */
  fp32t *m_magn;  /* 频强度测量值 */
  fp32t *s_freq;  /* 频率点合成值 */
  fp32t *s_magn;  /* 频强度合成值 */
} STFT_Tuner;


static fComplex fcpCvt2Cp(fp32t m, fp32t p)
{
    return fCPLX(m * cosf(p), m * sinf(p));
}

static void fcpCvt(fComplex x,
                     fp32t *mg, fp32t *pha)
{
    *mg = sqrtf(x.re * x.re + x.im * x.im);
    *pha = *mg > 0? atan2f(x.im, x.re) : 0;
}

static fp32t phase_mod(fp32t ph)
{
    int np = ph / M_PI;

    if (np >= 0)
        return ph - M_PI * (np + (np & 1));
    else
        return ph - M_PI * (np - (np & 1));
}

static void _tunr_trans_func(fComplex *arr,
                       int fsize, void *usr)
{
    STFT_Tuner *e = usr;
    fp32t freq_per_b, ph_step;
    int fftsize = e->ft.pfft->fftLen;

    fp32t (*tuner_fn)(STFT_Tuner *o,
                   void *usr) = e->tu_func;

    freq_per_b = (fp32t)e->sampr / fftsize;
    ph_step = M_PI * 2 * e->step / fftsize;

    fp32t *mea_freq = e->m_freq;
    fp32t *mea_magn = e->m_magn;
    fp32t *syn_freq = e->s_freq;
    fp32t *syn_magn = e->s_magn;

    for (int i = 0; i <= fsize; ++i)
    {
        fp32t pha, u;
        fcpCvt(arr[i], mea_magn + i, &pha);

        u = phase_mod(pha - e->pcache[i]
                  - i * ph_step) / ph_step;

        e->pcache[i] = pha;
        syn_freq[i] = 0, syn_magn[i] = 0;
        mea_freq[i] = (i + u) * freq_per_b;
    }

    fp32t TR = tuner_fn(e, e->tu_data);

    for (int i = 0, p; i <= fsize; ++i)
    {
        if ((p = i * TR + .5) <= fsize)
        {
            syn_magn[p] += mea_magn[i];
            syn_freq[p] = mea_freq[i] * TR;
        }
    }

    for (int i = 0; i <= fsize; ++i)
    {
        fp32t u = syn_freq[i]
                    / freq_per_b * ph_step;

        e->pcount[i] = u
             = phase_mod(u + e->pcount[i]);

        arr[i] = fcpCvt2Cp(syn_magn[i], u);
    }
}

void tunrInit(STFT_Tuner *e,
     FastDFT *fft, int stepL, int sampR)
{
    int fsize = fft->fftLen / 2 + 1;

    fp32t *mem = malloc(sizeof(fp32t)
                         * fsize * 6);

    e->pcache = mem + fsize * 0;
    e->pcount = mem + fsize * 1;
    e->m_freq = mem + fsize * 2;
    e->m_magn = mem + fsize * 3;
    e->s_freq = mem + fsize * 4;
    e->s_magn = mem + fsize * 5;
    stftInit(&e->ft, fft, stepL);
    e->step = stepL, e->sampr = sampR;
}

void tunrZero(STFT_Tuner *e)
{
    int size = e->ft.pfft->fftLen / 2;

    for (int i = 0; i <= size; ++i)
    {
        e->pcache[i] = 0.;
        e->pcount[i] = 0.;
    }

    e->fifo_cnt = 0, stftZero(&e->ft);
}

void tunrFree(STFT_Tuner *e)
{
    stftFree(&e->ft), free(e->pcache);
}

char tunrExec(STFT_Tuner *e, fp32t *audio,
   fp32t (*tuner)(STFT_Tuner *o, void *u),
         fp32t *audioT, void *tune_usr_data)
{
    int fftsize = e->ft.pfft->fftLen;
    int sft_step = e->step;

    e->tu_func = tuner;
    e->tu_data = tune_usr_data;

    if (!audio)
    {
        if (e->fifo_cnt < 1) { return 0; }

        for (int i = 0; i < sft_step; ++i)
            audioT[i] = 0.;

        stftExec(&e->ft, audioT,
                     _tunr_trans_func, e);
    }
    else
    {
        for (int i = 0; i < sft_step; ++i)
            audioT[i] = audio[i];

        stftExec(&e->ft, audioT,
                     _tunr_trans_func, e);

        if ((e->fifo_cnt += sft_step)
                  < fftsize) { return 0; }
    }

    return (e->fifo_cnt -= sft_step, 1);
}



/* ======== [ STFT 移调机 (立体声) ] ======== */

typedef struct {
  fp32t *sL, *sR;
  STFT_Tuner L, R;
} STFT_Tuner2;


void tunr2Init(STFT_Tuner2 *e,
           FastDFT *fft, int step, int samprate)
{
    e->sL = malloc(sizeof(fp32t) * step);
    e->sR = malloc(sizeof(fp32t) * step);
    tunrInit(&e->L, fft, step, samprate);
    tunrInit(&e->R, fft, step, samprate);
}

void tunr2Zero(STFT_Tuner2 *tu2)
{
    tunrZero(&tu2->L), tunrZero(&tu2->R);
}

void tunr2Free(STFT_Tuner2 *tu2)
{
    free(tu2->sL), free(tu2->sR);
    tunrFree(&tu2->L), tunrFree(&tu2->R);
}

void arr2chn(fp32t const (*au)[2],
                fp32t *Lchn, fp32t *Rchn, int N)
{
    for (int i = 0; i < N; ++i)
        Lchn[i] = au[i][0], Rchn[i] = au[i][1];
}

void chn2arr(fp32t const *Lchn, 
       fp32t const *Rchn, fp32t (*au)[2], int N)
{
    for (int i = 0; i < N; ++i)
        au[i][0] = Lchn[i], au[i][1] = Rchn[i];
}

char tunr2Exec(STFT_Tuner2 *e, fp32t (*S)[2],
        fp32t (*tuner)(STFT_Tuner *o, void *u),
             fp32t (*T)[2], void *tune_usr_data)
{
    int N = e->L.step, stat;

    if (!S)
    {
        stat = tunrExec(&e->L, NULL, tuner,
                         e->sL, tune_usr_data)
            && tunrExec(&e->R, NULL, tuner,
                         e->sR, tune_usr_data);
    }
    else
    {
        arr2chn(S, e->sL, e->sR, N);

        stat = tunrExec(&e->L, e->sL, tuner,
                         e->sL, tune_usr_data)
            && tunrExec(&e->R, e->sR, tuner,
                         e->sR, tune_usr_data);
    }

    return (chn2arr(e->sL, e->sR, T, N), stat);
}






static fp32t tfn(STFT_Tuner *o, void *td)
{
    return *(fp32t *)td;
}



#define FFT_LEN 2048
#define FFT_STP 128


int main()
{
    FILE *au_rd = popen("ffmpeg -v error -i StarSky.aac -ac 2 -ar 96000 -f f32le -", "rb");
    FILE *au_wr = popen("ffmpeg -y -f f32le -ac 2 -ar 96000 -i - tune.flac", "wb");

    FastDFT FFT;
    STFT_Tuner2 Tuner;
    fftInit(&FFT, FFT_LEN);
    tunr2Init(&Tuner, &FFT, FFT_STP, 96000);
    tunr2Zero(&Tuner);

    static fp32t aud[FFT_STP][2];
    static fp32t tun[FFT_STP][2];
    char run = 1;

    fp32t tr = powf(2, 0.5 / 12);

    while (run)
    {
        fp32t (*sp)[2] = NULL;
        run = 0;

        if (fread(aud, 8*FFT_STP, 1, au_rd) == 1)
            run |= 1, sp = aud;

        if (tunr2Exec(&Tuner, sp, tfn, tun, &tr))
        {
            for (int i = 0; i < FFT_STP; ++i)
            {
                aud[i][0] = tun[i][0] * 0.9;
                aud[i][1] = tun[i][1] * 0.9;
            }

            run |= 1;
            fwrite(aud, 8 * FFT_STP, 1, au_wr);
        }
    }

    fftFree(&FFT), tunr2Free(&Tuner);
    return (pclose(au_rd), pclose(au_wr), 0);
}
