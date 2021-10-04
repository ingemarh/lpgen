/****************************************************************************
pllow.c -- decodump routine for straight fft
*****************************************************************************/

#include <stdlib.h>
#include <string.h>		/* Included as we are using the bzero function */
#include <strings.h>
#include <unistd.h>
#include <math.h>

double get_tid (void);
void print_tid (double diff_time, int nth, char *decoder);

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
#ifdef __linux__
__inline__ unsigned long long int
gethrtime (void)
{
  unsigned long long int x = 0;
  __asm__ volatile (".byte 0x0f, 0x31":"=A" (x));
  return x;
}

#define CPU_SPEED 800.e6
#else
#define CPU_SPEED 1.e9
#endif
#include <pthread.h>
struct pth
{
  int nth;
  ulong th;
  ulong k1;
  ulong k2;
  ulong nr_gates;
  DSPCMPLXSHORT *in;
  DSPCMPLX *out;
  ulong frac;
  int fft_len;
  ulong nr_points;
  ulong n_deco;
  float *out_t;
  ulong nr_samp;
  ulong debug;
  fftwf_complex *inf;
  fftwf_plan p;
  fftwf_plan pb;
};

void *pllow_loop (struct pth *);

/* nbits is the number of integer values supplied by the user */
/* in is a pointer to the complex input data vector */
/* out is a pointer to complex output data vector */
int
pllow (ulong nbits, int *par, ulong fbits, float *fpar, DSPCMPLXSHORT * in,
       DSPCMPLX * out, double *upar)
/* pllow FFT version 3 fftw */
{
/* par: npar fft_len nr_gates nr_rep nr_samp frac*/

  float *out_tf;
  ulong i, j, k, k2, npar, maxthreads, nr_rep;
  struct pth pth;
  DSPCMPLX *ia;
#ifdef __linux__
  //unsigned long long int        start[]={0,0,0,0,0};
  double start[] = { 0, 0, 0, 0, 0 };

#else
  hrtime_t start[] = { 0, 0, 0, 0, 0 };
#endif
  //start[0]=gethrtime();
  start[0] = get_tid ();

  npar = par[0];
  pth.fft_len = par[1];
  pth.nr_gates = par[2];
  nr_rep = par[3];
  pth.nr_samp = par[4];
  pth.frac = par[5];
  pth.debug = par[6];
  if (pth.debug > 1)
    printf ("%ld %d %ld %f %d %lf\n", nbits, par[0], fbits, fpar[0], in[0].re,
	    out[0].re);
  maxthreads = par[7];
  in += par[8];
  pth.nr_points = pth.nr_samp / pth.nr_gates;
  if (pth.debug > 1)
    printf ("%ld %d %ld\n", pth.frac, pth.fft_len, pth.nr_points);

  pthread_attr_t sched_glob;
  void *retval;
  pthread_t *thread;
  struct pth *ptth;
  pth.nth = sysconf (_SC_NPROCESSORS_ONLN);
  if (pth.nth > maxthreads)
    pth.nth = maxthreads;
  if (pth.nth > 1)
    {
      pthread_attr_init (&sched_glob);
      pthread_attr_setscope (&sched_glob, PTHREAD_SCOPE_SYSTEM);
      thread = (pthread_t *) malloc (pth.nth * sizeof (pthread_t));
      ptth = (struct pth *) malloc (pth.nth * sizeof (struct pth));
    }

  pth.n_deco = pth.nr_gates * pth.fft_len;
  if (pth.debug > 1)
    printf ("%d %ld\n", pth.nth, pth.n_deco);
  pth.out_t = (float *) calloc (pth.nth * pth.n_deco, sizeof (float));

  pth.in = in;
  pth.out = out;
  pth.inf = fftwf_malloc (pth.nth * pth.fft_len * 2 * sizeof (fftwf_complex));
  pth.p =
    fftwf_plan_dft_1d (pth.fft_len, pth.inf, pth.inf + pth.fft_len,
		       FFTW_FORWARD, FFTW_ESTIMATE);
  pth.pb =
    fftwf_plan_dft_1d (pth.fft_len, pth.inf, pth.inf + pth.fft_len,
		       FFTW_BACKWARD, FFTW_ESTIMATE);

  //start[1]=gethrtime();
  start[1] = get_tid ();
  for (k = 0, i = 0; i < pth.nth; i++, k = k2)
    {
      k2 = nr_rep * (i + 1) / pth.nth;
      pth.th = i;
      pth.k1 = k;
      pth.k2 = k2;
      if (pth.debug > 1)
	printf ("%ld %ld %ld\n", i, k, k2);
      if (pth.nth > 1)
	{
	  memcpy (ptth + i, &pth, sizeof (struct pth));
	  pthread_create (thread + i, &sched_glob,
			  (void *(*)(void *)) pllow_loop,
			  (void *) (ptth + i));
	}
      else
	pllow_loop (&pth);
    }

  for (i = 0; i < pth.nth; i++)
    {
      if (pth.nth > 1)
	pthread_join (thread[i], &retval);
      if (i > 0)
	{
	  out_tf = pth.out_t + i * pth.n_deco;
	  for (j = 0; j < pth.n_deco; j++)
	    pth.out_t[j] += out_tf[j];
	}
    }
  //start[2]=gethrtime();
  start[2] = get_tid ();

  bzero (pth.inf, pth.fft_len * sizeof (fftwf_complex));
  out_tf = pth.out_t;
  for (i = 0; i < pth.nr_gates; i++)
    {
      ia = (DSPCMPLX *) pth.inf;
      for (j = 0; j < pth.fft_len; j++, ia++, out_tf++)
	ia->re = *out_tf;
      fftwf_execute_dft (pth.pb, pth.inf, pth.inf + pth.fft_len);
      ia = (DSPCMPLX *) (pth.inf + pth.fft_len);
      for (j = 0; j < pth.frac; j++, ia++)
	{
	  out[i + j * pth.nr_gates].re = ia->re / pth.fft_len;
	  out[i + j * pth.nr_gates].im = ia->im / pth.fft_len;
	}
    }

  if (pth.debug > 1)
    printf ("%ld\n", pth.frac);
  fftwf_destroy_plan (pth.p);
  fftwf_destroy_plan (pth.pb);
  fftwf_free (pth.inf);
  free (pth.out_t);
  if (pth.nth > 1)
    {
      free (thread);
      free (ptth);
    }
  if (pth.debug > 0)
    {
      //start[3]=gethrtime();
      start[3] = get_tid ();
      //printf("pllow fd_alt: %d threads: %.2g s",pth.nth,(start[3]-start[0])/CPU_SPEED);
      print_tid (start[3] - start[0], pth.nth, "pllow fd_alt");

      if (pth.debug > 1)
	for (j = 0; j < 3; j++)
	  printf (" %.2f", (start[j + 1] - start[j]) / CPU_SPEED);
      printf ("\n");
    }
  return 0;
}				/* pllow */


void *
pllow_loop (struct pth *pth)
{
  ulong k, j, i, th;
  DSPCMPLXSHORT *in;
  DSPCMPLX *ia, *ib;
  fftwf_complex *inf, *outf;
  float *out_t, *out_tf;

  th = pth->th;
  inf = pth->inf + 2 * th * pth->fft_len;
  outf = inf + pth->fft_len;
  ia = (DSPCMPLX *) inf;
  ib = (DSPCMPLX *) outf;
  bzero (inf + pth->nr_points,
	 (pth->fft_len - pth->nr_points) * sizeof (fftwf_complex));
  out_t = pth->out_t + th * pth->n_deco;

  for (k = pth->k1; k < pth->k2; k++)
    {
      in = pth->in + k * pth->nr_samp;
      out_tf = out_t;
      for (i = 0; i < pth->nr_gates; i++)
	{
	  for (j = 0; j < pth->nr_points; j++, in++)
	    {
	      ia[j].re = in->re;
	      ia[j].im = in->im;
	    }
	  fftwf_execute_dft (pth->p, inf, outf);
	  for (j = 0; j < pth->fft_len; j++, out_tf++)
	    *out_tf += (ib[j].re * ib[j].re + ib[j].im * ib[j].im);
	}
    }
  return NULL;
}
