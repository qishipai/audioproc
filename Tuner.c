#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
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



/* ========== [ STFT 移调机 ] ========== */

typedef struct {
  int step;       /* 单次处理点数 */
  int sampr;      /* 源音频采样率 */
  int fifo_cnt;   /* 队列缓冲点数 */
  FastDFT *pfft;
  fComplex *arr;  /* 变换值存储区 */
  fp32t *rd_buf;  /* 输入值缓冲区 */
  fp32t *wr_buf;  /* 输出值缓冲区 */
  fp32t *win_fn;  /* 窗函数值列表 */
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

static void au_fifo_roL(fp32t *aufifo,
      fp32t *audio_samp, int fL, int N)
{
    fp32t *fifo_prv_sp = aufifo + N;

    for (int i = N; i < fL; ++i)
        *aufifo++ = *fifo_prv_sp++; 

    if (!audio_samp)
        while (aufifo != fifo_prv_sp)
            *aufifo++ = 0.0;
    else
        while (aufifo != fifo_prv_sp)
            *aufifo++ = *audio_samp++;
}

void tunrInit(STFT_Tuner *e,
    FastDFT *fft, int stepL, int sampR)
{
    int size = fft->fftLen;

    fp32t  *m = malloc(sizeof(fp32t)
                    * (size * 7 + 6));

    e->win_fn = malloc(sizeof(fp32t)
                       * fft->fftLen);

    e->rd_buf = m + size * 2;
    e->wr_buf = m + size * 3;
    e->pcache = m + size * 4;

    for (int i = 0; i < size; ++i)
        e->win_fn[i] = 1.0 - cosf(
           M_PI * 2 * i / (size - 1));

    e->step = stepL, e->sampr = sampR;
    e->arr = (void *)m, e->pfft = fft;
    e->pcount = m + size *  9 / 2 + 1;
    e->m_freq = m + size * 10 / 2 + 2;
    e->m_magn = m + size * 11 / 2 + 3;
    e->s_freq = m + size * 12 / 2 + 4;
    e->s_magn = m + size * 13 / 2 + 5;
}

void tunrFree(STFT_Tuner *e)
{
    free (e->arr), free(e->win_fn);
}

void tunrZero(STFT_Tuner *e)
{
    int sz = e->pfft->fftLen*7 + 6;
    fp32t *memptr = (void *)e->arr;
    e->fifo_cnt = 0;

    for (int i = 0; i < sz; i += 2)
        memptr[i] = memptr[i + 1] = 0;
}

char tunrExec(STFT_Tuner *e, fp32t *audio,
   fp32t (*tuner)(STFT_Tuner *o, void *u),
         fp32t *audioT, void *tune_usr_data)
{
    int fftsize = e->pfft->fftLen;
    int size2 = fftsize / 2;
    fComplex *arr  = e->arr;
    fp32t *rd_buf = e->rd_buf;
    fp32t *wr_buf = e->wr_buf;
    fp32t *win_fn = e->win_fn;

    au_fifo_roL(e->rd_buf, audio,
                         fftsize, e->step);

    if (!audio)
    {
        if (e->fifo_cnt < 1) { return 0; }
    }
    else
    {
        if ((e->fifo_cnt += e->step)
                  < fftsize) { return 0; }
    }

    fp32t freq_per_b, ph_step, audios_amp;
    audios_amp = (fp32t)e->step  / fftsize;
    freq_per_b = (fp32t)e->sampr / fftsize;
    ph_step = M_PI * 2 * e->step / fftsize;

    fp32t *mea_freq = e->m_freq;
    fp32t *mea_magn = e->m_magn;
    fp32t *syn_freq = e->s_freq;
    fp32t *syn_magn = e->s_magn;

    for (int i = 0; i < fftsize; ++i)
    {
        arr[i].re = rd_buf[i] * win_fn[i];
        arr[i].im = 0;
    }

    fftExec(e->pfft, arr, 0);

    for (int i = 0; i <= size2; ++i)
    {
        fp32t pha, u;
        fcpCvt(arr[i], mea_magn + i, &pha);

        u = phase_mod(pha - e->pcache[i]
                  - i * ph_step) / ph_step;

        e->pcache[i] = pha;
        syn_freq[i] = 0, syn_magn[i] = 0;
        mea_freq[i] = (i + u) * freq_per_b;
    }

    fp32t TR = tuner(e, tune_usr_data);

    for (int i = 0, p; i <= size2; ++i)
    {
        if ((p = i * TR + .5) <= size2)
        {
            syn_magn[p] += mea_magn[i];
            syn_freq[p] = mea_freq[i] * TR;
        }
    }

    for (int i = 0; i <= size2; ++i)
    {
        fp32t u = syn_freq[i]
                    / freq_per_b * ph_step;

        e->pcount[i] = u
             = phase_mod(u + e->pcount[i]);

        arr[i] = fcpCvt2Cp(syn_magn[i], u);
    }

    for (int i = 1; i < size2; ++i)
        arr[fftsize - i] = fcpConj(arr[i]);

    fftExec(e->pfft, arr, 1);

    au_fifo_roL(e->wr_buf, NULL,
                         fftsize, e->step);

    for (int i = 0; i < fftsize; ++i)
        wr_buf[i] += arr[i].re * win_fn[i];

    for (int i = 0; i < e->step; ++i)
        audioT[i] = wr_buf[i] * audios_amp;

    return (e->fifo_cnt -= e->step, 1);
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

void tunr2Free(STFT_Tuner2 *tu2)
{
    free(tu2->sL), free(tu2->sR);
    tunrFree(&tu2->L), tunrFree(&tu2->R);
}

void tunr2Zero(STFT_Tuner2 *tu2)
{
    tunrZero(&tu2->L), tunrZero(&tu2->R);
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
    return *(fp32t*)td;
}



#define FFT_LEN 1024
#define FFT_STP 128


int main()
{
    FILE *au_rd = popen("ffmpeg -v error -i OaO.flac -ac 2 -ar 48000 -f f32le -", "rb");
    FILE *au_wr = popen("ffmpeg -y -f f32le -ac 2 -ar 48000 -i - out.flac", "wb");

    FastDFT FFT;
    STFT_Tuner2 Tuner;
    fftInit(&FFT, FFT_LEN);
    tunr2Init(&Tuner, &FFT, FFT_STP, 96000);
    tunr2Zero(&Tuner);

    static fp32t aud[FFT_STP][2];
    char run = 1;

    fp32t tr = powf(2, 2.0 / 12);

    while (run)
    {
        fp32t (*sp)[2] = NULL;
        run = 0;

        if (fread(aud, 8*FFT_STP, 1, au_rd) == 1)
            run |= 1, sp = aud;

        if (tunr2Exec(&Tuner, sp, tfn, aud, &tr))
        {
            for (int i = 0; i < FFT_STP; ++i)
                aud[i][0] *= 0.9, aud[i][1] *= 0.9;

            run |= 1;
            fwrite(aud, 8 * FFT_STP, 1, au_wr);
        }
    }

    fftFree(&FFT), tunr2Free(&Tuner);
    return (pclose(au_rd), pclose(au_wr), 0);
}
