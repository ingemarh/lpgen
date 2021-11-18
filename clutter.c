/****************************************************************************
clutter.c -- decodump routine for removing clutter via fft
*****************************************************************************/
#include <stdlib.h>
#include <string.h>		/* Included as we are using the bzero function */
#include <strings.h>
#include <unistd.h>
#include <math.h>

double get_tid (void);
void print_tid_clutter (double diff_time, int ncth, int nlth, char *decoder);

typedef struct
{
  float re;
  float im;
} DSPCMPLX;

typedef struct
{
  int16_t re;
  int16_t im;
} DSPCMPLXSHORT;

#include <fftw3.h>
#define CPU_SPEED 1.e9
#include <sys/time.h>
//##ifndef gethrtime /*gcc case*/
#ifndef NANOSEC			/*gcc case */
#include <time.h>
#define NANOSEC 1000000000
typedef unsigned long long int hrtime_t;
unsigned long
gethrtime (void)
{
  struct timespec ts;
  if (clock_gettime (CLOCK_MONOTONIC_RAW, &ts) != 0)
    return (-1);
  return ((ts.tv_sec * NANOSEC) + ts.tv_nsec);
}
#endif
#include <pthread.h>
struct pth
{
  ulong th;
  ulong k1;
  ulong k2;
  ulong nr_win;
  ulong nr_rep;
  int notch;
  DSPCMPLXSHORT *in;
  DSPCMPLX *out;
  ulong nr_pp;
  ulong nr_lag;
  ulong nr_clutter;
  int fft_clutter;
  ulong debug;
  fftwf_complex *inf;
  fftwf_plan pc;
  fftwf_plan pcb;
};

void *plwin_clutter (struct pth *);
void *plwin_lp (struct pth *);

int
clutter_6 (ulong nbits, int *par, ulong fbits, float *fpar,
	   DSPCMPLXSHORT * in, DSPCMPLX * out, double *upar)
{

/*printf("%ld %ld %ld %ld %ld %ld %ld\n",nbits,par,fbits,fpar,in,out,upar);*/
  ulong i, j, k, k2, clutthreads, float_notch = 0, ncth, nlth, nout;
  int nth;
  struct pth pth;
  //hrtime_t      start[]={0,0,0,0,0};
  double start[] = { 0, 0, 0, 0, 0 };
  //start[0]=gethrtime();
  start[0] = get_tid ();

  pth.nr_win = par[0];
  pth.nr_rep = par[1];
  pth.nr_pp = par[2];
  pth.nr_clutter = par[3];
  if (par[4] > -1)
    pth.notch = par[4];
  else
    pth.notch = (int) upar[0];
  clutthreads = par[5];
  pth.debug = par[6];
  pth.nr_lag = par[7];

  if (pth.notch < 0)
    float_notch = 1;
  pth.fft_clutter = pth.nr_rep / pth.nr_win;

  pthread_attr_t sched_glob;
  void *retval;
  pthread_t *thread = NULL;
  struct pth *ptth = NULL;
  nth = sysconf (_SC_NPROCESSORS_ONLN);
  ncth = clutthreads;
  if (ncth > pth.nr_clutter)
    ncth = pth.nr_clutter;
  nlth = pth.nr_lag;
  if (nlth > nth)
    nlth = nth;
  nth = (nlth > ncth) ? nlth : ncth;	/*max threads */
  if (nth > 1)
    {
      pthread_attr_init (&sched_glob);
      pthread_attr_setscope (&sched_glob, PTHREAD_SCOPE_SYSTEM);
      thread = (pthread_t *) malloc (nth * sizeof (pthread_t));
      ptth = (struct pth *) malloc (nth * sizeof (struct pth));
    }

  nout = pth.nr_pp * pth.nr_lag - pth.nr_lag * (pth.nr_lag - 1) / 2;
  if (nout == 0)
    nout = 1;
  bzero (out, nout * sizeof (DSPCMPLX));
  pth.in = in;
  pth.out = out;
  if (ncth > 0)
    {
      pth.inf =
	fftwf_malloc (ncth * 2 * pth.fft_clutter * sizeof (fftwf_complex));
      pth.pc =
	fftwf_plan_dft_1d (pth.fft_clutter, pth.inf,
			   pth.inf + pth.fft_clutter, FFTW_FORWARD,
			   FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
      pth.pcb =
	fftwf_plan_dft_1d (pth.fft_clutter, pth.inf + pth.fft_clutter,
			   pth.inf, FFTW_BACKWARD,
			   FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    }

  //start[1]=gethrtime();
  start[1] = get_tid ();
  for (k = 0, i = 0; i < ncth; i++, k = k2)
    {
      k2 = pth.nr_clutter * (i + 1) / ncth;
      pth.th = i;
      pth.k1 = k;
      pth.k2 = k2;
      if (ncth > 1)
	{
	  memcpy (ptth + i, &pth, sizeof (struct pth));
	  pthread_create (thread + i, &sched_glob,
			  (void *(*)(void *)) plwin_clutter,
			  (void *) (ptth + i));
	}
      else
	plwin_clutter (&pth);
    }
  if (ncth > 1)
    for (i = 0; i < ncth; i++)
      {
	pthread_join (thread[i], &retval);
	if (i == 0)
	  pth.notch = ptth[i].notch;
	else
	  pth.notch += ptth[i].notch;
	if (pth.debug > 0)
	  printf ("%d", ptth[i].notch);
      }
  if (nlth > 0)
    {
      for (k = 0, i = 0; i < nlth; i++, k = k2)
	{
	  k2 = pth.nr_lag * (i + 1) / nlth;
	  pth.th = i;
	  pth.k1 = k;
	  pth.k2 = k2;
	  if (nlth > 1)
	    {
	      memcpy (ptth + i, &pth, sizeof (struct pth));
	      pthread_create (thread + i, &sched_glob,
			      (void *(*)(void *)) plwin_lp,
			      (void *) (ptth + i));
	    }
	  else
	    plwin_lp (&pth);
	}
      if (nlth > 1)
	for (i = 0; i < nlth; i++)
	  pthread_join (thread[i], &retval);
    }

  if (ncth > 0)
    {
      out[0].im = (float) pth.notch / (float) ncth;
      fftwf_destroy_plan (pth.pc);
      fftwf_destroy_plan (pth.pcb);
      fftwf_free (pth.inf);
    }
  if (nth > 1)
    {
      free (thread);
      free (ptth);
    }
  if (pth.debug > 0)
    {
      //start[4]=gethrtime();
      start[4] = get_tid ();
      //printf("clutter: [%ld %ld] threads: %.2g s",ncth,nlth,(start[4]-start[0])/CPU_SPEED);
      print_tid_clutter (start[4] - start[0], ncth, nlth, "clutter");
      if (ncth > 0 && float_notch)
	printf (" Notch: %.2g", out[0].im);
      if (pth.debug > 1)
	for (j = 0; j < 4; j++)
	  printf (" %.2f", (start[j + 1] - start[j]) / CPU_SPEED);
      printf ("\n");
    }
  return 0;
}				/* clutter_6 */

void *
plwin_clutter (struct pth *pth)
{
  ulong k, j, i, win;
  int notch;
  DSPCMPLXSHORT *in;
  DSPCMPLX *clx, *cly;
  float clf, clm, *clt;

  clx = (DSPCMPLX *) pth->inf + pth->th * pth->fft_clutter * 2;
  cly = clx + pth->fft_clutter;
  i = pth->nr_rep / pth->nr_win;
  notch = pth->notch;
  if (notch < 0)
    {
      clt = (float *) calloc (pth->fft_clutter, sizeof (float));
      for (j = pth->k1; j < pth->k2; j++)
	for (win = 0; win < pth->nr_win; win++)
	  {
	    for (k = win; k < pth->nr_rep; k += pth->nr_win)
	      {
		in = pth->in + k * pth->nr_pp + j;
		clx[k / pth->nr_win].re = in->re;
		clx[k / pth->nr_win].im = in->im;
	      }
	    bzero (clx + i, (pth->fft_clutter - i) * sizeof (DSPCMPLX));
	    fftwf_execute_dft (pth->pc, (fftwf_complex *) clx,
			       (fftwf_complex *) cly);
	    for (k = 0; k < pth->fft_clutter; k++)
	      clt[k] += (cly[k].re * cly[k].re + cly[k].im * cly[k].im);
	  }
      notch = 0;
      clf = 2. * clt[0];
      clm = 0;
      j = pth->fft_clutter - 1;
      for (k = 1; k < pth->fft_clutter; k++)
	clm += clt[k];
      while (notch < -pth->notch + 1 && clf / 2. > clm / j)
	{
	  notch++;
	  clf = clt[notch] + clt[pth->fft_clutter - notch];
	  clm -= clf;
	  j -= 2;
	}
      notch--;
      free (clt);
      pth->notch = notch;
    }
  if (notch >= 0)
    {
      clf =
	1. / sqrt (1. -
		   (float) (2 * notch +
			    1) / (float) i) / (float) pth->fft_clutter;
      for (j = pth->k1; j < pth->k2; j++)
	for (win = 0; win < pth->nr_win; win++)
	  {
	    for (k = win; k < pth->nr_rep; k += pth->nr_win)
	      {
		in = pth->in + k * pth->nr_pp + j;
		clx[k / pth->nr_win].re = in->re;
		clx[k / pth->nr_win].im = in->im;
	      }
	    bzero (clx + i, (pth->fft_clutter - i) * sizeof (DSPCMPLX));
	    fftwf_execute_dft (pth->pc, (fftwf_complex *) clx,
			       (fftwf_complex *) cly);
	    bzero (cly, (notch + 1) * sizeof (DSPCMPLX));
	    bzero (cly + pth->fft_clutter - notch, notch * sizeof (DSPCMPLX));
	    fftwf_execute_dft (pth->pcb, (fftwf_complex *) cly,
			       (fftwf_complex *) clx);
	    for (k = win; k < pth->nr_rep; k += pth->nr_win)
	      {
		in = pth->in + k * pth->nr_pp + j;
		in->re = rint (clx[k / pth->nr_win].re * clf);
		in->im = rint (clx[k / pth->nr_win].im * clf);
	      }
	  }
    }

  return NULL;
}

void *
plwin_lp (struct pth *pth)
{
  ulong i, k, j;
  DSPCMPLXSHORT *in;
  DSPCMPLX *out, s0, s1;

  for (i = 0; i < pth->nr_rep; i++)
    {
      in = pth->in + i * pth->nr_pp;
      out = pth->out + pth->nr_pp * pth->k1 - pth->k1 * (pth->k1 - 1) / 2;
      for (j = pth->k1; j < pth->k2; j++)
	for (k = 0; k < pth->nr_pp - j; k++)
	  {
	    s0.re = in[k].re;
	    s0.im = in[k].im;
	    s1.re = in[k + j].re;
	    s1.im = in[k + j].im;
	    out->re += (s0.re * s1.re + s0.im * s1.im);
	    /*out->im += (s0.im * s1.re - s0.re * s1.im);*/
	    out->im += (s0.re * s1.im - s1.re * s0.im);
	    out++;
	  }
    }
  return NULL;
}

void
matface (int *par, int *nin, double *in_r, double *in_i,
	 int *nout, double *out_r, double *out_i, double *upar)
{
  ulong nb = 8;
  float *fpar = NULL;
  int i;
  DSPCMPLXSHORT *in;
  DSPCMPLX *out;
  in = (DSPCMPLXSHORT *) malloc (*nin * sizeof (DSPCMPLXSHORT));
  for (i = 0; i < *nin; i++)
    {
      in[i].re = in_r[i];
      in[i].im = in_i[i];
    }
  out = (DSPCMPLX *) malloc (*nout * sizeof (DSPCMPLX));
  clutter_6 (nb, par, 0, fpar, in, out, upar);
  for (i = 0; i < *nout; i++)
    {
      out_r[i] = out[i].re;
      out_i[i] = out[i].im;
    }
  for (i = 0; i < *nin; i++)
    {
      in_r[i] = in[i].re;
      in_i[i] = in[i].im;
    }
  free (in);
  free (out);
}

typedef struct
{
  double re;
  double im;
} DSPCMPLXDBL;

void
pyface (int *par, int nin, DSPCMPLXDBL * in_r, DSPCMPLX * out, double *upar)
{
  ulong nb = 8;
  float *fpar = NULL;
  int i;
  DSPCMPLXSHORT *in;
  in = (DSPCMPLXSHORT *) malloc (nin * sizeof (DSPCMPLXSHORT));
  for (i = 0; i < nin; i++)
    {
      in[i].re = in_r[i].re;
      in[i].im = in_r[i].im;
    }
  decoder_6 (nb, par, 0, fpar, in, out, upar);
  for (i = 0; i < nin; i++)
    {
      in_r[i].re = in[i].re;
      in_r[i].im = in[i].im;
    }

  free (in);
}
