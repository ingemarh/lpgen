#include <sys/types.h>		/*Needed for int16_t */
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <assert.h>
#include <pthread.h>

double get_tid (void);
void print_tid (double diff_time, int nth, char *decode);

#define         DSP_NODECO      0
#define         DSP_ACDECO      1
#define         DSP_TAILDECO    2
#define         MAX_NR_THREADS  4 
typedef struct			/* Make a new internal data type */
{
  float re;
  float im;
}
DSPCMPLX;

typedef struct			/* Make a new internal data type */
{
  int16_t re;
  int16_t im;
} DSPCMPLXSHORT;

typedef struct
{
  ulong rd_us;
  ulong rd_nbytes;
  ulong deco_us;
  ulong deco_flops;
} __USER__DSPRATE;

typedef struct
{
  __USER__DSPRATE rate;
  int tracelevel;
  int recording;
  char expid[64];
  char expname[64];
  char msg[128];
} DACCINFO;

typedef struct
{
  ulong rawsa;
  ulong fstlag;
  ulong maxlag;
  ulong nfract;
  ulong decomode;
  ulong zprflen1;
  ulong rmsa;
  ulong ndeco;
  ulong accum;
  ulong nbits;
  ulong vec_len;
  ulong nr_ipp;
  ulong nr_code_sets;
  ulong sub_int;
  ulong code_nr;
  ulong nr_threads;
  int thread_nr;
  int ipp;
  DSPCMPLX *lag_res1;
  DSPCMPLX *lag_res2;
  DSPCMPLXSHORT *in;
  DSPCMPLX *out_res;
  DSPCMPLX *out;
  DACCINFO *info;

  int *dpar;
  int *code;
  float Decofir[2][128];
} DECO;



int __user__decolen (int nsamples, int maxlag, int codelen, int nfract,
		     int taildeco);

void *__user__decoder_int_thread (DECO * calc);

void *__user__decoder_thread (DECO * deco_send);

void
__user__lag_profiling (DSPCMPLXSHORT * data, DSPCMPLX * res, ulong msize,
		       ulong max_lag, ulong lag_incr, int zer_away,
		       int loops);
long __user__decodoer (DECO * deco, DSPCMPLX * rawSA, DSPCMPLX * dataSA,
		       DACCINFO * info);

static int __user__i_deco (DECO * deco,
			   DSPCMPLX * rawSA,
			   DSPCMPLX * dataSA, ulong * flops, char *msg);

static int __user__ft_deco (DECO * deco,
			    DSPCMPLX * rawSA,
			    DSPCMPLX * dataSA, ulong * flops, char *msg);

static int __user__f_deco (DECO * deco,
			   DSPCMPLX * rawSA,
			   DSPCMPLX * dataSA, ulong * flops, char *msg);

static void __user__getfilter (float *h, int *nh, int *code, int lag,
			       int nbits);
static void __user__filter (DSPCMPLX * yy0, float *h0, DSPCMPLX * x0,
			    long *ny, long nh, long nx, int accum, int decim);

int
alter_dec (unsigned long dbits,
	   int *dpar, unsigned long fbits, float *fpar, DSPCMPLXSHORT * in,
	   DSPCMPLX * out, double *upar)
{
  int nr_ipp, nlags, k, sub_int, vec_len, nsamples,
    nr_code_sets, ret_code, code_nr, j, jjj, nr_slices = 1;

  DECO *deco_send[MAX_NR_THREADS];
  DSPCMPLX *lag_res1, *lag_res2, *out_res;
  pthread_t *threads;
  pthread_attr_t thread_sched;
  DECO *deco;
  DACCINFO *info;
  float tid = 0.0;
  double start = 0;
  deco = malloc (sizeof (DECO));
  info = malloc (sizeof (DACCINFO));

  /*for (k = 0; k < 7; k++)
    printf (" %d", dpar[k]);
  printf ("\n");*/
  deco->nbits = dpar[0];
  nr_ipp = dpar[1];
  vec_len = dpar[2];
  deco->nfract = dpar[3];
  deco->maxlag = (dpar[4] - deco->nfract + 1) / deco->nfract;
  sub_int = dpar[5];
  nr_code_sets = dpar[6];
  if (nr_code_sets < 0)
    {
      nr_code_sets *= -1;
      nr_slices = nr_ipp / nr_code_sets;
      nr_ipp = nr_code_sets;
    }
  deco->decomode = DSP_TAILDECO;
  deco->zprflen1 = (deco->nbits - 1) * deco->nfract + vec_len - 1;
  deco->rmsa = 0;
  deco->fstlag = 1;
  nlags = deco->maxlag * deco->nfract + deco->nfract - 1;
  deco->ndeco = __user__decolen (vec_len, nlags, deco->nbits, deco->nfract,
				 deco->decomode);
  deco->vec_len = vec_len;
  deco->nr_ipp = nr_ipp;
  deco->nr_code_sets = nr_code_sets;
  deco->sub_int = sub_int;
  deco->dpar = dpar;
  /*printf ("%ld %d %d %ld %ld %d %d %d\n", deco->nbits, nr_ipp, vec_len,
	  deco->nfract, deco->maxlag, sub_int, nr_code_sets, nr_slices);*/
  start = get_tid ();
  nsamples = deco->zprflen1 - (deco->nbits - 1) * deco->nfract + 1;
  if ((threads =
       (pthread_t *) malloc (MAX_NR_THREADS * sizeof (pthread_t))) == NULL)
    {
      perror ("Malloc threads");
      exit (EXIT_FAILURE);
    }


  if (nr_code_sets > MAX_NR_THREADS)
    {
      out_res = malloc (MAX_NR_THREADS * deco->ndeco * sizeof (DSPCMPLX));
      lag_res1 =
	malloc (MAX_NR_THREADS * (nlags + 1) * vec_len * sizeof (DSPCMPLX));
      /* Will hold the zero added lagprofile which will decoded */
      lag_res2 =
	malloc (MAX_NR_THREADS * (nlags + 1) * deco->zprflen1 *
		sizeof (DSPCMPLX));
      for (k = 0; k < MAX_NR_THREADS; k++)
	{
	  deco_send[k] = (DECO *) malloc (sizeof (DECO));
	}
    }
  else
    {
      out_res = malloc (nr_code_sets * deco->ndeco * sizeof (DSPCMPLX));
      lag_res1 =
	malloc (nr_code_sets * (nlags + 1) * vec_len * sizeof (DSPCMPLX));
      /* Will hold the zero added lagprofile which will decoded */
      lag_res2 =
	malloc (nr_code_sets * (nlags + 1) * deco->zprflen1 *
		sizeof (DSPCMPLX));
      for (k = 0; k < nr_code_sets; k++)
	{
	  deco_send[k] = (DECO *) malloc (sizeof (DECO));
	}
    }
  ret_code = pthread_attr_init (&thread_sched);
  if (ret_code != 0)
    {
      printf ("pthread_attr_init error=%d", ret_code);
      exit (EXIT_FAILURE);
    }
  if (pthread_attr_setinheritsched (&thread_sched, PTHREAD_INHERIT_SCHED) !=
     0)
     {
     perror ("pthread_attr_setinheritsched thread_sched");
     }

     if (pthread_attr_setschedpolicy (&thread_sched, SCHED_FIFO) != 0)
     {
     perror ("pthread_attr_setchedpolicy thread_sched");
     } 
  ret_code = pthread_attr_setscope (&thread_sched, PTHREAD_SCOPE_SYSTEM);
  if (ret_code != 0)
    {
      printf ("pthread_attr_setscope error=%d", ret_code);
      exit (EXIT_FAILURE);
    }


  for (jjj = 0; jjj < nr_slices; jjj++)
    {
      bzero (out, sizeof (DSPCMPLX) * deco->ndeco);
      deco->out = out;
      code_nr = 0;
      deco->ipp = -1;
      deco->nr_threads = MAX_NR_THREADS;
      for (j = 0; j < (int) (nr_code_sets / MAX_NR_THREADS); j++)
	{
	  for (k = 0; k < MAX_NR_THREADS; k++)
	    {
	      deco->ipp++;
	      memcpy (deco_send[k], deco, sizeof (DECO));
	      deco_send[k]->lag_res1 = lag_res1 + (nlags + 1) * vec_len * k;
	      deco_send[k]->lag_res2 =
		lag_res2 + (nlags + 1) * deco->zprflen1 * k;
	      deco_send[k]->in = in;
	      deco_send[k]->out_res = out_res + deco->ndeco * k;
	      deco_send[k]->info = info;
	      deco_send[k]->code_nr = code_nr++;
	      ret_code =
		pthread_create (&(threads[k]), &thread_sched,
				(void *(*)(void *)) __user__decoder_thread,
				(void *) deco_send[k]);
	      if (ret_code != 0)
		{
		  printf ("pthread create 1 error= %d\n", ret_code);
		  exit (EXIT_FAILURE);
		}
	    }
	  for (k = 0; k < MAX_NR_THREADS; k++)
	    {
	      /*
	       * Wait for all channels to be ready if running
	       * threaded
	       */
	      ret_code = pthread_join (threads[k], NULL);
	      if (ret_code != 0)
		{
		  printf ("pthread_join 1 error=%d\n", ret_code);
		  exit (EXIT_FAILURE);
		}
	    }
	  for (k = 0; k < MAX_NR_THREADS; k++)
	    {
	      deco->thread_nr = k;
	      memcpy (deco_send[k], deco, sizeof (DECO));
	      deco_send[k]->out_res = out_res;
	      deco_send[k]->out = out;
	      deco_send[k]->info = info;
	      ret_code =
		pthread_create (&(threads[k]), &thread_sched,
				(void *(*)(void *))
				__user__decoder_int_thread,
				(void *) deco_send[k]);
	      if (ret_code != 0)
		{
		  printf ("pthread create 1 error= %d\n", ret_code);
		  exit (EXIT_FAILURE);
		}

	    }
	  for (k = 0; k < MAX_NR_THREADS; k++)
	    {
	      /*
	       * Wait for all channels to be ready if running
	       * threaded
	       */
	      ret_code = pthread_join (threads[k], NULL);
	      if (ret_code != 0)
		{
		  printf ("pthread_join 1 error=%d\n", ret_code);
		  exit (EXIT_FAILURE);
		}
	    }
	}
      deco->nr_threads =
	nr_code_sets -
	((int) (nr_code_sets / MAX_NR_THREADS)) * MAX_NR_THREADS;
      for (k = 0; k < deco->nr_threads; k++)
	{
	  deco->ipp++;
	  memcpy (deco_send[k], deco, sizeof (DECO));
	  deco_send[k]->lag_res1 = lag_res1 + (nlags + 1) * vec_len * k;
	  deco_send[k]->lag_res2 =
	    lag_res2 + (nlags + 1) * deco->zprflen1 * k;
	  deco_send[k]->in = in;
	  deco_send[k]->out_res = out_res + deco->ndeco * k;
	  deco_send[k]->info = info;
	  deco_send[k]->code_nr = code_nr++;
	  ret_code =
	    pthread_create (&(threads[k]), &thread_sched,
			    (void *(*)(void *)) __user__decoder_thread,
			    (void *) deco_send[k]);
	  if (ret_code != 0)
	    {
	      printf ("pthread create 1 error= %d\n", ret_code);
	      exit (EXIT_FAILURE);
	    }

	}
      for (k = 0; k < deco->nr_threads; k++)
	{
	  /*
	   * Wait for all channels to be ready if running
	   * threaded
	   */
	  ret_code = pthread_join (threads[k], NULL);
	  if (ret_code != 0)
	    {
	      printf ("pthread_join 1 error=%d\n", ret_code);
	      exit (EXIT_FAILURE);
	    }
	}
      for (k = 0; k < deco->nr_threads; k++)
	{
	  deco->thread_nr = k;
	  memcpy (deco_send[k], deco, sizeof (DECO));
	  deco_send[k]->out_res = out_res;
	  deco_send[k]->out = out;
	  deco_send[k]->info = info;
	  ret_code =
	    pthread_create (&(threads[k]), &thread_sched,
			    (void *(*)(void *)) __user__decoder_int_thread,
			    (void *) deco_send[k]);
	  if (ret_code != 0)
	    {
	      printf ("pthread create 1 error= %d\n", ret_code);
	      exit (EXIT_FAILURE);
	    }
	}
      for (k = 0; k < deco->nr_threads; k++)
	{
	  /*
	   * Wait for all channels to be ready if running
	   * threaded
	   */
	  ret_code = pthread_join (threads[k], NULL);
	  if (ret_code != 0)
	    {
	      printf ("pthread_join 1 error=%d\n", ret_code);
	      exit (EXIT_FAILURE);
	    }
	}
      out += deco->ndeco;
      in += (vec_len * nr_code_sets);
      /*print_tid (get_tid () - start, MAX_NR_THREADS, "alt_decoder");
      printf ("\n");*/
    }

  free (deco);
  free (out_res);
  free (info);
  free (lag_res1);
  free (threads);
  free (lag_res2);
  if (nr_code_sets > MAX_NR_THREADS)
    {
      for (k = 0; k < MAX_NR_THREADS; k++)
	{
	  free (deco_send[k]);
	}
    }
  else
    {
      for (k = 0; k < nr_code_sets; k++)
	{
	  free (deco_send[k]);
	}
    }
  print_tid (get_tid() - start, deco->nr_threads, "alt_decoder");
  printf ("\n");
  return 0;
}

/****************************************************************************
///////////// decodoer //////////////////////////////////////////////////
*****************************************************************************
** A DECOlist based decoder. This is dispatcher routine, the actual
** work is done in the xx_deco routines, see below.
**
** 22-Aug-2000 Jm
** 22-Aug-2000 AW changed gettimeofday to gethrtime
****************************************************************************/
long
__user__decodoer (DECO * deco, DSPCMPLX * rawSA, DSPCMPLX * dataSA,
		  DACCINFO * info)
{
  ulong flops;
  if (deco->accum == 0)
    {
      info->rate.deco_us = 0;
      info->rate.deco_flops = 0;
      info->rate.rd_us = 0;
      info->rate.rd_nbytes = 0;
    }
  flops = 0;

  if (deco == NULL)
    {
      strcpy (info->msg, "decodoer(): no operation");
      return 0;
    }

  if (deco->decomode == DSP_NODECO)
    {
      strcpy (info->msg, "decodoer(): Decomode NODECO not allowed");
      return -1;
    }

  if (deco->nfract < 2)
    {				/*integer-only lags */
      if (__user__i_deco (deco, rawSA, dataSA, &flops, info->msg))
	{
	  return -1;
	}
    }
  else				/* fractional lags */
    {
      if (deco->decomode == DSP_TAILDECO)
	{
	  if (__user__ft_deco (deco, rawSA, dataSA, &flops, info->msg))
	    {
	      return -1;
	    }
	}
      else if (deco->decomode == DSP_ACDECO)
	{
	  if (__user__f_deco (deco, rawSA, dataSA, &flops, info->msg))
	    {
	      return -1;
	    }
	}
      else
	{
	  sprintf (info->msg, "decodoer():Unsupported decomode (%ld)",
		   deco->decomode);
	  return -1;
	}
    }

  info->rate.deco_flops += flops;

  /* WHEN decoding is ready -- ALL codes handled -- return number of points
   * written to the result memory. This happens when accum mode == 2.
   * Return 0 before that.
   */

  if (deco->accum == 2)
    {
      return deco->ndeco;
    }
  else
    {
      return 0;
    }

}				/*decodoer */

/****************************************************************************
///////////// i_deco ////////////////////////////////////////////////////////
*****************************************************************************
** A DECOlist based decoder for INTEGER LAGS. Also assumes lagstep = 1;
**
** 1 Apr, 1998 Jm
** 7-Jun-1998 Jm DECO version 1.1 taken into use ( profstep removed etc).
** 22-Aug-2000 Jm Simpified for KST.
****************************************************************************/
int
__user__i_deco (DECO * deco, DSPCMPLX * rawSA,	/* this points to the beginning of the intermediate  raw data */
		DSPCMPLX * dataSA,	/* this points to the beginning of result memory */
		ulong * flops, char *msg)
{
  int ntaps;
  int rawproflen;
  int profstep;
  long decoproflen;
  int lag;
  DSPCMPLX *raw_profile;
  DSPCMPLX *out_profile;
  float *decofir;

  raw_profile = rawSA;		/* + deco->rawsa; */
  rawproflen = deco->zprflen1;
  out_profile = dataSA + deco->rmsa;

  if (deco->decomode == DSP_TAILDECO)
    profstep = 2;
  else if (deco->decomode == DSP_ACDECO)
    profstep = 1;
  else
    {
      strcpy (msg, "i_deco(): illegal decomode");
      return 1;
    }

  decofir = &(deco->Decofir[0][0]);

  for (lag = deco->fstlag; lag <= deco->maxlag; lag++)
    {

      __user__getfilter (decofir, &ntaps, deco->code, lag, deco->nbits);

      __user__filter (out_profile,
		      decofir,
		      raw_profile, &decoproflen, ntaps, rawproflen,
		      deco->accum, 1);

      *flops += decoproflen * 2 * (ntaps + (deco->accum > 0));

      raw_profile += rawproflen;
      rawproflen -= profstep;
      out_profile += decoproflen;

    }				/*for */

  return 0;

}				/*i_deco */


/****************************************************************************
///////////// ft_deco ///////////////////////////////////////////////////////
*****************************************************************************
** 
** A DECOlist based decoder for FRACTIONAL LAGS, with TAIL DECODING.
** Also it is assumed that lagincrement from profile to profile is 1.
** OUTPUT DATA ORDER is as is assumed by GUISDAP experiment design package.
** This routine filters a single code, all lags, accumulating the output
** on top of previous codes.
** The GUISDAP output profile order e.g. in the case of "fractionality" 3 is
**
**  1/3[1] 2/3[1] 1 4/3[1] 4/3[2] 5/3[1] 5/3[2] 2 ... N N+1/3[N] N+2/3[N].
**
** Here the number in brackets [m] indicates which integer lag is used for
** decoding. 
**
** 15-Aug-1999 Jm
** 22-Aug-2000 Jm Simplified for KST.
****************************************************************************/
int
__user__ft_deco (DECO * deco, DSPCMPLX * rawSA,	/* this points to the beginning of the CURRENT rawdata */
		 DSPCMPLX * dataSA,	/* this points to the beginning of WHOLE result memory */
		 ulong * flops, char *msg)
{
  int m;
  int fir1, incount = 0;
  int ntaps0;
  int ntaps1;
  int rawproflen;
  long decoproflen;
  int lag;
  int nfract;
  DSPCMPLX *raw_profile;
  DSPCMPLX *out_profile;
  long outcount;
  float *decofir0;
  float *decofir1;
  int n_profiles;

  nfract = deco->nfract;

  /* Pointer to the start of input profile (including zeros infront) */

  raw_profile = rawSA;		/* + deco->rawsa; */
  rawproflen = deco->zprflen1;

  /*
   *  Pointer to the result memory
   */

  out_profile = dataSA + deco->rmsa;
  outcount = 0;

  /*
   *  Initialize the two decoding filters that we need for fractional lags: 
   *  decofir0 for decoding "from below" and decofir1 for decoding 
   *  "from above". At each new integer lag we then only need to compute
   *  one new filter.
   */

  lag = 1;			/* the first fract lags are decoded by lag 1 coefficients */

  fir1 = 1;			/* this is the row(-index) in the Decofir array */
  /* which currently contains decofir1 */

  decofir0 = NULL;		/* the first fract lags are not decoded from below */
  ntaps0 = 0;

  decofir1 = deco->Decofir[fir1];
  __user__getfilter (decofir1, &ntaps1, deco->code, lag, deco->nbits);

  /* 
   * Loop though all input profiles. Each input profile produces one or
   * two output profiles.
   */

  n_profiles = ((deco->maxlag) + 1) * nfract - 1;

  for (m = 1; m <= n_profiles; m++)
    {
      if (m % nfract == 0)
	{

	  lag = m / nfract;

	  /* Handle an integer lag profile

	   * First update the two decoding filters. The filter that earlier
	   * was for decoding from above, now becomes for decoding from below,
	   * and a new filter is computed for decoding from above.
	   */

	  decofir0 = deco->Decofir[fir1];
	  ntaps0 = ntaps1;
	  fir1 = 1 - fir1;
	  decofir1 = deco->Decofir[fir1];
	  if (m == (deco->nbits - 1) * nfract)
	    {
	      /* 
	       * Beyond the last integer lag, we will only decode from below 
	       */

	      decofir1 = NULL;
	    }
	  else
	    {
	      __user__getfilter (decofir1, &ntaps1, deco->code, lag + 1,
				 deco->nbits);
	    }

	  /*
	   *  At integer lag boundary, the decoding filter shortens by
	   *  nfract, so the number of leading zeros drops by nfract, and
	   *  so the rawproflen which includes those zeros also
	   *  suddenly decreases by that amount.
	   */

	  rawproflen -= nfract;

	  /* 
	   *  The integer lag is decoded with that lag's coefficients
	   */

	  __user__filter (out_profile,
			  decofir0,
			  raw_profile,
			  &decoproflen,
			  ntaps0, rawproflen, deco->accum, deco->nfract);

	  *flops += decoproflen * 2 * (ntaps0 + (deco->accum > 0));

	  raw_profile += rawproflen;
	  rawproflen--;
	  out_profile += decoproflen;
	  outcount += decoproflen;
	}
      else
	{
	  /* 
	   * Handle fractional lag profile.
	   *
	   * First decode it from below
	   */
	  if (decofir0)
	    {
	      __user__filter (out_profile,
			      decofir0,
			      raw_profile,
			      &decoproflen,
			      ntaps0, rawproflen, deco->accum, deco->nfract);

	      *flops += decoproflen * 2 * (ntaps0 + (deco->accum > 0));

	      out_profile += decoproflen;
	      outcount += decoproflen;

	    }

	  /* 
	   * Then decode it from above. We need to skip over nfract
	   * "unnecessary" zeros in front of the profile. Those zeros
	   * are needed only for decoding from below, where the decoding
	   * filter is that much longer.
	   */

	  if (decofir1)
	    {
	      __user__filter (out_profile,
			      decofir1,
			      raw_profile + nfract,
			      &decoproflen,
			      ntaps1, rawproflen - nfract, deco->accum,
			      deco->nfract);

	      *flops += decoproflen * 2 * (ntaps1 + (deco->accum > 0));

	      out_profile += decoproflen;
	      outcount += decoproflen;

	    }

	  raw_profile += rawproflen;
	  rawproflen--;

	}

      incount += rawproflen + 1;
    }				/*loop back for the next input profile */

  assert (outcount == deco->ndeco);

  return 0;
}				/*ft_deco_gup */


/****************************************************************************
///////////// f_deco ////////////////////////////////////////////////////////
*****************************************************************************
** 
** A DECOlist based decoder for FRACTIONAL LAGS,  WITHOUT tail decoding.
** Also it is assumed that lagincrement from profile to profile is 1.
**
** 8-Jun-1998 Jm
** 22-Aug-2000 Jm: For KST.
****************************************************************************/
int
__user__f_deco (DECO * deco, DSPCMPLX * rawSA,	/* this points to the beginning of the CURRENT raw data */
		DSPCMPLX * dataSA, ulong * flops, char *msg)
{
  int k;
  int ntaps;
  int rawproflen;
  int rawproflen0;
  long decoproflen;
  int lag;
  int nfract;
  DSPCMPLX *raw_profile;
  DSPCMPLX *raw_profile0;
  DSPCMPLX *out_profile;

  nfract = deco->nfract;

  raw_profile0 = rawSA;		/* + deco->rawsa; */
  rawproflen0 = deco->zprflen1;
  out_profile = dataSA + deco->rmsa;

  for (lag = deco->fstlag; lag <= deco->maxlag; lag++)
    {
      /*
       * Note that "lag" here means a "full lag", lag=timelag/BAUDlength.
       * Each full lag is used to decode a group of 2*nfract-1 profiles,
       * surrounding the full lag profile. For instance, with nfract = 3,
       * lag 6 coefficients are used to decode the profiles 
       * 5+1/3  5+2/3  6  6+1/3  6+2/3
       */
      raw_profile = raw_profile0;
      rawproflen = rawproflen0;

      __user__getfilter (deco->Decofir[0], &ntaps, deco->code, lag,
			 deco->nbits);

      for (k = 1; k <= 2 * nfract - 1; k++)
	{

	  __user__filter (out_profile,
			  deco->Decofir[0],
			  raw_profile,
			  &decoproflen, ntaps, rawproflen, deco->accum,
			  deco->nfract);

	  *flops += decoproflen * 2 * (ntaps + (deco->accum > 0));

	  raw_profile += rawproflen;
	  rawproflen--;
	  out_profile += decoproflen;

	  if (k == nfract)
	    {
	      /*
	       * Raw_profile now points to the beginning of the first
	       * fractional profile after the full lag profile
	       * (to the "6+1/3" profile when lag=6, nfract=3).
	       * We need to decode that profile also with the NEXT full lag
	       * coefficients, so we save its start address and length now.
	       */
	      raw_profile0 = raw_profile;
	      rawproflen0 = rawproflen;
	    }

	}			/*loop back for next k */

    }				/*loop back to the next full lag */

  return 0;
}				/*f_deco */


/***************************************************************************
//////////////////////// filter ////////////////////////////////////////////
****************************************************************************
**
**  y0 points to the beginning of the complex output array
**  h0 points to the beginning of the real coefficient array
**  x0 points to the beginning of the complex input data array
**  nx number of complex points in the x array
**  nh number of real points in the h array
**  ny number of complex points in the y array
**  decim input decimation factor >= 1.  Use decim = 1 for integer lags.
**  complex data ordered as [real,imag,real,imag...].
**
**  Algorithm: sliding average with h, taking only every decim'th x
**  for a fiven output point. In the example, decim = 3 : 
**  
**          1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
**          x x x x x x x x x x x x x x x x
**          a _ _ b _ _ c                     ==> 1'st output point y[0]
**            a _ _ b _ _ c                   ==> 2'nd output point y[1]
**                         ...                ...
**                            a _ _ b _ _ c   ==> last output point y[ny-1]
**
**          y[0] = a*x1 + b*x3 + c*x7
**          y[1] = a*x2 + b*x4 + c*x8
**              etc
** 
**  a=h(1), b=h(2), c=h(3).
**
**  NOTE. This would be a normal input decimating filter IF we would now
**  compute only every decim'th output point. However, we want to
**  get finer spatial gate SEPARATION (though not spatial RESOLUTION), 
**  and therefore keep all available output points.
**  
**  9-Jun-1998 Jm
**  22-Aug-2000 Jm: For KST. Using floats instead of doubles for ceofficients.
**                  No integer <--> float conversions needed anymore.
**************************************************************************/
void
__user__filter (DSPCMPLX * yy0,
		float *h0,
		DSPCMPLX * x0, long *ny, long nh, long nx, int accum,
		int decim)
{
  DSPCMPLX *x;
  DSPCMPLX acc;
  float *h;
  int n, k;

  /*  
   * Number of output points 
   */

  *ny = nx - (nh - 1) * decim;

  /*  
   * The filter 
   */

  for (n = 0; n < *ny; n++)
    {

      x = x0;
      h = h0;

      /*  Dot product of the FIR and a segment of input data,
       *  using only every decim'th of the input points.
       */

      acc.re = *h * (x->re);
      acc.im = *h * (x->im);
      x += decim;
      h++;

      for (k = 1; k < nh; k++)
	{
	  acc.re += *h * (x->re);
	  acc.im += *h * (x->im);
	  h++;
	  x += decim;
	}

      if (accum > 0)
	{
	  yy0->re += acc.re;
	  yy0->im += acc.im;
	  yy0++;
	}
      else
	{
	  *yy0++ = acc;
	}

      x0++;
    }				/*loop back to compute next output point */
  return;
}				/*filter */

/***************************************************************************
/////////////// getfilter //////////////////////////////////////////////////
****************************************************************************
** Return the decoder filter for an alternating code lag L profile.
** The filter is essentially the the code's autocorrelation for that lag.
**
**  For example, if the code sign sequence in order of transmission is
**      code = [s0 s2 s3 s4 s5 s6 s7 s7]    ( nbits = 8 )
**  and lag = 2, then
**      h = [s0*s2 s1*s3 ... s5*s7]         ( nh = 6 )
**
** 8-Jun-1998 Jm.
** 22-Aug-2000 Jm : For KST, filter returned in floats, not doubles.
**************************************************************************/
void
__user__getfilter (float *h, int *nh, int *code, int lag, int nbits)
{
  int *c1, *c2;
  int i;

  c1 = code;
  c2 = c1 + lag;

  *nh = nbits - lag;

  for (i = 0; i < *nh; i++)
    {
      *h++ = (float) (*c1++) * (float) (*c2++);
    }
  return;
}				/* getfilter */

int
__user__decolen (int nsamples, int maxlag, int codelen, int nfract,
		 int taildeco)
{
  int n, niprof, nn, nraw, p1len, p2len, F = nfract, fstilag = 1;
  int lastilag = maxlag / nfract, N = nsamples;

  if (taildeco == DSP_TAILDECO)
    {
      niprof = lastilag - fstilag + 1;
      p1len = N - fstilag * F;
      p2len = N - lastilag * F;
      n = (2 * F - 1) * niprof * (p1len + p2len) / 2;
      if (lastilag < codelen - 1)
	{
	  nraw = N - lastilag * F;
	  nn = ((F - 1) * (2 * nraw - F)) / 2;
	  n = n + nn;
	}
    }
  else
    {
      F = nfract;
      n = (lastilag - fstilag + 1) * (2 * F - 1) * (N - (codelen - 1) * F);
      if (lastilag < codelen - 1)
	{
	  nn = (F - 1) * (2 * N + F * (3 - 2 * codelen)) / 2;
	  n = n + nn;
	}
    }

  return n;
}

void
__user__lag_profiling (DSPCMPLXSHORT * data, DSPCMPLX * res, ulong msize,
		       ulong max_lag, ulong lag_incr, int zer_away, int loops)
{
  int j, n, k, m, lp;
  float a, b, c, d, tmpr, tmpi;
  DSPCMPLXSHORT *datat;

  msize /= loops;
  for (lp = 0; lp < loops; lp++)
    {
      datat = data + lp * msize;
      for (j = 0; j < 2 * (msize - max_lag + 1); j += 2)
	{
	  a = (float) (datat + j / 2)->re;
	  b = (float) (datat + j / 2)->im;
	  for (n = -2 * msize + j - 2 * lag_incr * zer_away, k = 0;
	       k < max_lag * 2; k += 2 * lag_incr)
	    {
	      c = (float) (datat + (j + k) / 2)->re;
	      d = (float) (datat + (j + k) / 2)->im;
	      n += 2 * msize - (k - 2 * lag_incr) * zer_away;
	      tmpr = a * c + b * d;
	      tmpi = a * d - b * c;
	      (res + n / 2)->re += tmpr;
	      (res + n / 2)->im += tmpi;
	    }
	}
      for (m = 0, j = 2 * (msize - max_lag + 1); j < 2 * msize; j += 2)
	{
	  a = (float) (datat + j / 2)->re;
	  b = (float) (datat + j / 2)->im;
	  for (m++, n = -2 * msize + j - 2 * lag_incr * zer_away, k = 0;
	       k < (max_lag - m) * 2; k += 2 * lag_incr)
	    {
	      c = (float) (datat + (j + k) / 2)->re;
	      d = (float) (datat + (j + k) / 2)->im;
	      n += 2 * msize - (k - 2 * lag_incr) * zer_away;
	      tmpr = a * c + b * d;
	      tmpi = a * d - b * c;
	      (res + n / 2)->re += tmpr;
	      (res + n / 2)->im += tmpi;
	    }
	}
    }
}


void *
__user__decoder_thread (DECO * deco_send)
{
  static const int success = 0;
  ulong k, nlags, vec_len, nr_ipp, nr_code_sets, sub_int;
  ulong nsamples;
  int i, j, lag, nzeros;
  DSPCMPLXSHORT *in_data, *in;
  DSPCMPLX *lag_res1, *lag_res2, *to, *profsa;

  j = deco_send->code_nr;
  nlags = deco_send->maxlag * deco_send->nfract + deco_send->nfract - 1;
  vec_len = deco_send->vec_len;
  nr_ipp = deco_send->nr_ipp;
  sub_int = deco_send->sub_int;
  nr_code_sets = deco_send->nr_code_sets;
  in = deco_send->in;
  lag_res1 = deco_send->lag_res1;
  lag_res2 = deco_send->lag_res2;
  nsamples =
    deco_send->zprflen1 - (deco_send->nbits - 1) * deco_send->nfract + 1;

  bzero (lag_res1, sizeof (DSPCMPLX) * (nlags + 1) * vec_len);
  bzero (lag_res2, sizeof (DSPCMPLX) * (nlags + 1) * deco_send->zprflen1);

  for (k = 0; k < (nr_ipp / nr_code_sets) / sub_int; k++)
    {
      for (i = 0; i < sub_int; i++)
	{
	  in_data =
	    in + j * vec_len * sub_int +
	    k * vec_len * sub_int * nr_code_sets + i * vec_len;
	  __user__lag_profiling (in_data, lag_res1, nsamples, nlags + 1, 1,
				 0, 1);
	}
    }

  to = lag_res2;
  profsa = lag_res1 + nsamples;
  for (lag = 1; lag <= nlags; lag++)
    {
      /* Prepend zeros in front of the lag profile for tail decoding */
      nzeros =
	deco_send->nfract * (deco_send->nbits - 1 -
			     (lag / deco_send->nfract));
      bzero (to, nzeros * sizeof (DSPCMPLX));
      to += nzeros;
      /* Copy the profile proper */
      memcpy (to, profsa, (nsamples - lag) * sizeof (DSPCMPLX));
      /* Move to beginning of next profile (skipping trailing zeros) */
      to += nsamples - lag;
      profsa += nsamples;
    }

  if (j == nr_code_sets - 1)
    deco_send->accum = 2;
  else if (j == 0)
    deco_send->accum = 0;
  else
    deco_send->accum = 1;
  deco_send->code = deco_send->dpar + 7 + j * deco_send->nbits;

  bzero (deco_send->out_res, deco_send->ndeco * sizeof (DSPCMPLX));
  __user__decodoer (deco_send, lag_res2, deco_send->out_res, deco_send->info);
  return ((void *) &success);
}

void *
__user__decoder_int_thread (DECO * calc)
{
  static const int success = 0;
  int k;
  DSPCMPLX *to, *res;
  ulong nsamples, samples, max_nr_thread, thread_nr, j;

  thread_nr = calc->thread_nr;
  max_nr_thread = calc->nr_threads;
  nsamples = calc->ndeco;
  samples = (int) (nsamples / max_nr_thread);
  if (thread_nr == (max_nr_thread - 1))
    {
      samples =
	(int) (nsamples / max_nr_thread) + nsamples -
	((int) (nsamples / max_nr_thread)) * max_nr_thread;
    }
  to = calc->out + thread_nr * (int) (nsamples / max_nr_thread);
/*  printf
    ("thread nr %d nsamples %d max_nr_thread %d  samples %d \n", thread_nr, nsamples, max_nr_thread, samples);*/
  res = calc->out_res + thread_nr * (int) (nsamples / max_nr_thread);

  for (k = 0; k < max_nr_thread; k++)
    {
      for (j = 0; j < samples; j++)
	{
	  (to + j)->re += (res + j + k * nsamples)->re;
	  (to + j)->im += (res + j + k * nsamples)->im;
	}
    }
  return ((void *) &success);
}

void
matface (int *par, int *nin, double *in_r, double *in_i,
	 int *nout, double *out_r, double *out_i, double *upar)
{
  ulong nb = 7 + par[0] * par[6];
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
  alter_dec (nb, par, 0, fpar, in, out, upar);
  for (i = 0; i < *nout; i++)
    {
      out_r[i] = out[i].re;
      out_i[i] = out[i].im;
    }
  free (in);
  free (out);
}
