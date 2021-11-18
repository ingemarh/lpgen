/****************************************************************************
plwin.c -- decodump routine for decoding fft
*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <math.h>

void print_tid (double diff_time, int nth, char *decoder);
double get_tid ();

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

#ifdef AOCL
#include "/home/kstdev/amd/aocl/2.0/amd-fftw/include/fftw3.h"
#elif MKL
//#include <mkl_dfti.h>
#include "fftw3_mkl.h"
//fftw3_mkl.verbose = 0;
#else
#include <fftw3.h>
#endif
#define CPU_SPEED 1.e9
#ifdef __linux__
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#define NANOSEC 1000000000
unsigned long
gethrtime (void)
{
  struct timespec ts;
  if (clock_gettime (CLOCK_MONOTONIC_RAW, &ts) != 0)
    return (-1);
  return ((ts.tv_sec * NANOSEC) + ts.tv_nsec);
}
#endif
#define cl_min(a,b)	(((a) < (b)) ? (a) : (b))
#include <pthread.h>
struct pth
{
  int nth;			/* ew: Number of threads */
  ulong th;			/* ew: Is increased when starting serveral threads, leads to different values beeing set when calling the fft functions */
  ulong k1;			/* ew: Is used for segmenting itne intervals between k1 and k2 when calling the fft-functions (i.e. clutter, plwin_loop, ...)  */
  /* ew: Is in clutter used for segmenting the intervals of pth.in that are used */
  ulong k2;			/* ew: See comment on k1 */
  ulong nr_win;			/* Number of different subsets in the code (number of data windows). ew - data window equals one coded pulse */
  int *windows;			/* ew:   */
  ulong win_len;		/* Number of bits in the code */
  ulong nr_gates;		/* ew:   */
  ulong res_mult;
  ulong nr_rep;
  ulong nout;
  ulong noud;
  DSPCMPLXSHORT *in;
  DSPCMPLX *out;
  ulong nr_undec;
  ulong undec1;
  ulong undec2;
  ulong frac;
  int nex;
  int fft_len;
  ulong nr_lags;
  float *out_t;
  ulong nr_samp;
  ulong nr_pp;
  ulong do_pp;
  ulong nr_blk;
  ulong nr_clutter;
  int notch;
  int fft_clutter;
  ulong nr_ugates;
  ulong lower_tail;
  ulong upper_tail;
  ulong nr_virtsamp;
  ulong debug;
  int fft_len_2;
  int fft_dlayer;
  ulong dl_short;
  ulong dl_long;
  ulong nr_dgates;
  DSPCMPLX *fir_samp;
  fftwf_complex *inf;
  fftwf_plan p;
  fftwf_plan p_2;
  fftwf_plan pd;
  fftwf_plan pc;
  fftwf_plan pb;
  fftwf_plan pb_2;
  fftwf_plan pdb;
  fftwf_plan pcb;
};

void *plwin_clutter (struct pth *);
void *plwin_loop (struct pth *);
void *pldec_loop (struct pth *);
void *pldlayer (struct pth *);



/* nbits is the number of integer values supplied by the user */
/* in is a pointer to the complex input data vector */
/* out is a pointer to complex output data vector */
int
decoder_6 (ulong nbits, int *par, ulong fbits, float *fpar,
	   DSPCMPLXSHORT * in, DSPCMPLX * out, double *upar)
/* decoder_6 FFT version 3 fftw */
{
/* par: npar fft_len nr_gates nr_undec nr_win nr_rep nr_samp win_len frac res_mult nr_cal*/
/* the rest is windows (nbits=nfft*nr_win+npar) */
/* 
* ew: 	nbits is set from prface; nb=p[0]+p[4]*p[7], is parameters + (number pulses in pulse group (aka "data windwos" aka "subsets in code"))*(length of bitcode)
*   	*par is parameters for experiment and analysis
*		fbits is set from prface to a hard 0
* 		*fpar always seems to be set to NULL
* 		*in is set by prface as in=(DSPCMPLXSHORT *)malloc(nin*sizeof(DSPCMPLXSHORT)), 	nin=(p[6]+p[10])*p[5];
*			that is ((samples per subcycle (aka pulse))+(extra samples in the power profile))*(total number of subcycles)
*			The input signal is mixed to remove the carrier frequency, the spectrum is around 0 at input to this file
* 		*out is set by prface as out=(DSPCMPLX *)malloc(nout*sizeof(DSPCMPLX)), 
*			where nout=((p[2]+1)*p[8]+p[2]*p[1]/2+p[3]*p[8])*p[9]+(p[6]+p[10])*p[9]+p[19]*(p[17]+p[18]+1);
*		*upar is set by prface as upar=(double *)malloc(20*sizeof(double)); upar[0]=0;
*/

  float *out_tf;
  ulong i, j, k, k2, npar, n_deco, maxthreads, ninf, float_notch =
    0, nth_clutt;
  struct pth pth;		/* ew: Instance of struct pth, named pth */
#ifdef __linux__
  //unsigned long long int        start[]={0,0,0,0,0}; /* ew: Array for timtestams for evalutaing program runtime */
  double start[] = { 0, 0, 0, 0, 0 };	/* ew: Array for timtestams for evalutaing program runtime */

#else
  hrtime_t start[] = { 0, 0, 0, 0, 0 };
#endif
  //start[0]=gethrtime(); 
  start[0] = get_tid ();

  /* ew: Comments on pth below from file plwin_par-fileinfo.txt */
  /* Total number of lines in par-file: 24+nr_win*win_len */
  npar = par[0];		/* Number of parameters in the parameter file. ew - Does not include bit code */
  pth.fft_len = par[1];		/* Size of fft */
  pth.nr_gates = par[2];	/* Number of range gates */
  pth.nr_undec = par[3];	/* Number of undecoded gates */
  pth.nr_win = par[4];		/* Number of different subsets in the code (number of data windows). ew - data window equals one coded pulse */
  /* Number of pulses per pulese group */
  pth.nr_rep = par[5];		/* Total number of subcycles in one integration period.
				 * how many pulses in total, each pulse group is repeted 
				 * = nr pulses per group * reptetitions of the pulse group */
  pth.nr_samp = par[6];		/* Number of samples per subcycle */
  pth.win_len = par[7];		/* Number of bits in the code */
  pth.frac = par[8];		/* Number of samples per code bit */
  pth.res_mult = par[9];	/* Number of subsets of code (normally 1) */
  pth.nr_pp = pth.nr_samp + par[10];	/* par[10] = nr_pp - nr_samp            Number of extra samples in the power profile */
  /* nr_pp = Number of samples in the power profile = samples per subcycle (aka data window, pulse) + extra samples */
  pth.lower_tail = par[11] * pth.frac;	/* lower_tail / frac    Number of padded zeros at lower edge, code bit units */
  pth.upper_tail = par[12] * pth.frac;	/* upper_tail / frac    Number of padded zeros at upper edge, code bit units */
  pth.undec1 = par[13];		/* Start gate for undecoded */
  pth.undec2 = par[14];		/* Stop gate for undecoded */
  pth.debug = par[15];		/* If >0: Display debugging messages, If >1: Display more debugging messages */
  if (pth.debug > 1)
    printf ("%ld %ld\n", nbits, fbits);
  maxthreads = par[16];		/* Maximum number of threads in the calcuclations */
  pth.dl_short = par[17];	/* Number of short D-layer lags */
  pth.dl_long = par[18];	/* Number of long D-layer lags */
  pth.nr_dgates = par[19];	/* Number of D-layer gates (including lower tail) */
  pth.nr_blk = pth.nr_pp;	/* Offset in input vector slices */
  if (par[20] > 0)
    pth.nr_blk += par[20];	/* Offset in input vector slices */
  else
    in += par[20];		/* Offset in input data vector */
  if (pth.debug > 1)
    printf ("%ld %d\n", pth.nr_blk, par[20]);

  pth.do_pp = par[21];		/* If non-zero: Calculate power profile */
  pth.nr_clutter = par[22];	/* Number of points to de-clutter */
  if (par[23] > -1)
    pth.notch = par[23];	/* Filter width for clutter reduction (actual width: 2*notch+1). If negative: Use value provided through EROS */
  else
    pth.notch = (int) upar[0];	/* ew: upar[0] is set to 0 from prface before calling decode6, --> pth.notch = 0 */
  if (pth.notch < 0)
    float_notch = 1;



  /* 
   * ew: description of output format from file plwin_par-fileinfo
   * 
   * Output format:
   * Short description                    Size (number of elements)
   * ------------------------------------------------------------
   * Normal short lags                    [(nr_gates+1)*frac]
   * Normal long lags                     [nr_gates*fft_len/2]
   * Normal undecoded lags                        [nr_undec*frac]
   * Normal power profile                 [do_pp*nr_pp]
   * Pulse-to-pulse lags                  [nr_dgates*dl_short]
   * Pulseset-to-pulseset lags            [nr_dgates*dl_long]
   * Pulse-to-pulse coherent profile      [nr_dgates]
   * 
   * The first four are repeated if res_mult>1.
   * 
   * Ordering of lags:
   * [lag1 range1], [lag1 range2] , ..., [lag1 rangeMAX], [lag2 range1], [lag2 range2], ...
   */


  /* ew: calculate length of output from calculations */
  pth.fft_dlayer = 2 * (pth.dl_long + 1) * (pth.dl_short + 1);
  pth.fft_clutter = pth.nr_rep / pth.nr_win;	/* ew: This is the number of repetitions of the transmitted pulse-group */
  pth.nex = pth.fft_len - pth.win_len * pth.frac;
  pth.nr_virtsamp = pth.nr_samp + pth.lower_tail + pth.upper_tail;
  pth.nr_ugates = pth.nr_virtsamp / pth.frac;
  pth.nr_lags = (pth.fft_len + 1) / 2;
/* Dec ffts done as a proper ACF (nfft>=2*sample)*/
  /*if(nex>0) pth.fft_len_2=2*pth.fft_len/pth.win_len;
     else */
  pth.fft_len_2 = 2 * pth.frac;


  /* ew: nout is total length of output from calculations? */
  pth.nout = pth.fft_len * pth.nr_gates;	/* First high res spectra */
  pth.nout += pth.fft_len_2 * pth.nr_ugates;	/* Then low res undecoded ones */
  pth.nout += pth.fft_len_2 * (pth.nr_gates + 1);	/* And decoded low res */
  if (pth.do_pp > 0)
    pth.nout += pth.nr_pp;	/* Last a power profile */
  if (pth.debug > 1)
    printf ("%ld %d %ld %d\n", pth.nout, pth.nex, pth.frac, pth.fft_len_2);



  /* ew: Setup multiple threads for the calculations */
  pthread_attr_t sched_glob;
  void *retval;
  pthread_t *thread = NULL;
  struct pth *ptth = NULL;	/* ew: Pointer to memory location with memory for (number of threds) x (size of struct pth) */
  pth.nth = sysconf (_SC_NPROCESSORS_ONLN);	/* Number of threads set from system capacity */
  if (pth.nth > maxthreads)
    pth.nth = maxthreads;	/* Number of threads set to a lower value i limited lower by input parameter */
  if (pth.nth > pth.nr_win)
    pth.nth = pth.nr_win;	/* Number of threads set to a lower value i limited lower by input parameter */
  if (pth.nth > 1)
    {
      pthread_attr_init (&sched_glob);
      pthread_attr_setscope (&sched_glob, PTHREAD_SCOPE_SYSTEM);
      thread = (pthread_t *) malloc (pth.nth * sizeof (pthread_t));
      ptth = (struct pth *) malloc (pth.nth * sizeof (struct pth));
    }


  /* ew: Find length of noud. Unclear what it is for? */
  if (pth.nr_undec > pth.nr_ugates)
    pth.nr_undec = pth.nr_ugates;
  pth.noud = (pth.nr_gates + 1) * pth.frac;	/* First the short lags */
  pth.noud += pth.nr_gates * pth.nr_lags;	/* Then the long lags */
  pth.noud += pth.nr_undec * pth.frac;	/* Finally some undecoded ones */
  if (pth.do_pp > 0)
    pth.noud += pth.nr_pp;	/* Last a power profile */
  if (pth.debug > 1)
    printf ("%ld %d\n", pth.noud, pth.nth);

  /* ew: allocate memory for out_t */
  n_deco = pth.res_mult * pth.nout;
  pth.out_t = (float *) calloc (pth.nth * n_deco, sizeof (float));



  /* ew: Set pointers in pth to bitcode, input, and output. */
  pth.windows = par + npar;	/* ew: par i pointer to input parameters, npar i number of params exluding bit code --> windows becomes pointer to start of bit-code? */
  pth.in = in;			/* ew: In declared as pointer to DSPCMPLXSHORT. Could probably be pointer to input data signal to be precessed. Pointer offsetted according to input params */
  pth.out = out;		/* Pointer to output after processing */



  /* ew: calculate ninf. ninf is used for setting memory allocation of variable inf which is used as temporary memory in different functions */
  ninf = 2 * pth.fft_len + pth.fft_len_2;
  if (pth.nr_dgates > 0)
    {
      pth.fir_samp =
	calloc (pth.fft_dlayer * pth.nr_dgates, sizeof (DSPCMPLX));
      if (pth.fft_dlayer > ninf)
	ninf = pth.fft_dlayer;
      if (2 * pth.fft_dlayer > pth.nth * ninf)
	ninf = 2 * pth.fft_dlayer;	/* ew: Use the largest needed size of inf */
    }
  if (pth.nr_clutter > 0 && 2 * pth.fft_clutter > pth.nth * ninf)
    ninf = 2 * pth.fft_clutter;	/* ew: Use the largest needed size of inf  */


  /* ew: Set up specifications for all the ffts that should be executed. Are stored in the pth struct which are pointed to when running the different functions */
  pth.inf = fftwf_malloc (pth.nth * ninf * sizeof (fftwf_complex));	/* ew: Allocate memory for inf from ninf*(number of threads). inf declared as fftwf_complex *inf, is used by fftwf lib. */
  /* ew: Different segments of memory allocated to inf is used as both input and output of the fft calculations */
  /* ew: Could be that output from one fft can be used as input for another? Depends on segments? */
  /* ew: Ingemar -> inf is most likely not used for moving results between functions */

  /* ew: Usage of inf memory allocation (see the fftwf_execute() function calls)
   *
   *                                                      INPUT POINTER                                   N                                               OUTPUT POINTER
   * plwin_clutter - pc           inf + th*fft_clutter*2                  fft_clutter                             [input] + fft_clutter                   |       th is increased from 0 up to nth in loop for starting serveral threads 
   * plwin_clutter - pcb          [output pc]                                             fft_clutter                             [input pc]
   *
   * plwin_loop   - p                     inf+th*(fft_len*2+fft_len_2)    fft_len                                 [input p]+fft_len       
   * plwin_loop   - p_2           [start p] + fft_len*2                   fft_len_2                               [output p]
   *
   * ... todo the rest. Also checking if/how the memory is reset between the functions if applicable */



  /* Define plans for fft:s that can later be executed */
  /* ew: example from manual, fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE) */
  /* ew: N = size of transform, in and out are pointers, (can be equal for in-place transformation) */
  /* ew: FFTW_ESTIMATE is a flag that means no testing of optimisation is done */
  /* ew:Order of execution of fft:s; clutter -> plwin_loop -> pldec_loop --> pldlayer */

  /* ew:                                                                                          OBS, inputs and outputs defined here are overidden when the ffts are executed */
  /* ew:                                                  N                                       in                                              out                                                     direction               flags                                                                                                                     */
  pth.p = fftwf_plan_dft_1d (pth.fft_len, pth.inf, pth.inf + pth.fft_len, FFTW_FORWARD, FFTW_ESTIMATE);	/* ew: executed in plwin_loop */
  pth.p_2 = fftwf_plan_dft_1d (pth.fft_len_2, pth.inf + pth.fft_len * 2, pth.inf + pth.fft_len, FFTW_FORWARD, FFTW_ESTIMATE);	/* ew: executed in plwin_loop */

  pth.pb = fftwf_plan_dft_1d (pth.fft_len, pth.inf, pth.inf + pth.fft_len, FFTW_BACKWARD, FFTW_ESTIMATE);	/* ew: executed in pldec_loop */
  pth.pb_2 = fftwf_plan_dft_1d (pth.fft_len_2, pth.inf, pth.inf + pth.fft_len, FFTW_BACKWARD, FFTW_ESTIMATE);	/* ew: executed in pldec_loop */

  if (pth.nr_dgates > 0)
    {
      pth.pd = fftwf_plan_dft_1d (pth.fft_dlayer, pth.inf, pth.inf + pth.fft_dlayer, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);	/* executed in pldlayer */
      pth.pdb = fftwf_plan_dft_1d (pth.fft_dlayer, pth.inf, pth.inf + pth.fft_dlayer, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);	/* executed in pldlayer */
    }

  /* ew: Run function for clutter */
  if (pth.nr_clutter > 0)
    {
      if (float_notch)
	nth_clutt = 1;		/* ew: Number of threads for cluttering, use one thread if filter width is set negative, =  filter widh from EROS */
      else
	nth_clutt = 1; //pth.nth;	/* ew: otherwise use standard n threads */

      /*ew                                            N                                       in                                              out                                                     direction               flags           */
      pth.pc = fftwf_plan_dft_1d (pth.fft_clutter, pth.inf, pth.inf + pth.fft_clutter, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);	/* executed in plwin_clutter */
      pth.pcb = fftwf_plan_dft_1d (pth.fft_clutter, pth.inf, pth.inf + pth.fft_clutter, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);	/* executed in plwin_clutter */

      for (k = 0, i = 0; i < nth_clutt; i++, k = k2)
	{			/* ew: run plwin_clutter pointed at pth in one or many threads. some of the fftwf_plan_dft_1d defined in pth will likely be used in plwin_clutter */
	  k2 = pth.nr_clutter * (i + 1) / nth_clutt;	/* ew: splits up the number of clutters to be calcualted into segments defined by k1 and k2 */
	  /* ew: k is put to k2 at transition to next loop iteration. */
	  /* ew: pth.k1 and pth.k2 are set in a similar way before each calling of the function for fft evaluation */
	  pth.th = i;
	  pth.k1 = k;
	  pth.k2 = k2;
	  if (nth_clutt > 1)
	    {
	      memcpy (ptth + i, &pth, sizeof (struct pth));
	      pthread_create (thread + i, &sched_glob,
			      (void *(*)(void *)) plwin_clutter,
			      (void *) (ptth + i));
	    }
	  else
	    plwin_clutter (&pth);
	}
      if (nth_clutt > 1)	/* ew: Collect results if more than one thread was used for the clutter processing */
	for (i = 0; i < nth_clutt; i++)
	  {
	    pthread_join (thread[i], &retval);
	    if (i == 0)
	      pth.notch = ptth[i].notch;
	    else
	      pth.notch += ptth[i].notch;
	    if (float_notch && pth.debug > 0)
	      printf ("%d", ptth[i].notch);
	  }
      upar[1] = (float) pth.notch / (float) nth_clutt;
    }

  /* ew: Run function plwin_loop */
  //start[1]=gethrtime(); /* ew: record timestamt to use for processing time evaluation */
  start[1] = get_tid ();	/* ew: record timestamt to use for processing time evaluation */
  for (k = 0, i = 0; i < pth.nth; i++, k = k2)
    {				/* ew: run plwin_loop pointed at pth in one or many threads. some of the fftwf_plan_dft_1d defined in pth will likely be used */
      k2 = pth.nr_win * (i + 1) / pth.nth;
      pth.th = i;
      pth.k1 = k;
      pth.k2 = k2;
      if (pth.nth > 1)
	{
	  memcpy (ptth + i, &pth, sizeof (struct pth));
	  pthread_create (thread + i, &sched_glob,
			  (void *(*)(void *)) plwin_loop,
			  (void *) (ptth + i));
	}
      else
	plwin_loop (&pth);
    }

  for (i = 0; i < pth.nth; i++)
    {				/* ew: Collect results if more than one thread was used for the plwin_loop processing */
      if (pth.nth > 1)
	{
	  pthread_join (thread[i], &retval);
	}
      if (i > 0)
	{
	  out_tf = pth.out_t + i * n_deco;
	  for (j = 0; j < n_deco; j++)
	    pth.out_t[j] += out_tf[j];
	}
    }


  /* ew: Run function pldec_loop */
  //start[2]=gethrtime(); /* ew: record timestamt to use for processing time evaluation */
  start[2] = get_tid ();	/* ew: record timestamt to use for processing time evaluation */
  for (k = 0; k < pth.res_mult; k++)
    {				/* ew: run pldec_loop pointed at pth in one or many threads. some of the fftwf_plan_dft_1d defined in pth will likely be used */
      pth.k1 = k;
      pth.th = 0;
      if (pth.nth > 1 && pth.nth >= pth.res_mult)
	{
	  pth.th = k;
	  memcpy (ptth + k, &pth, sizeof (struct pth));
	  pthread_create (thread + k, &sched_glob,
			  (void *(*)(void *)) pldec_loop,
			  (void *) (ptth + k));
	}
      else
	pldec_loop (&pth);
    }

  for (k = 0; k < pth.res_mult; k++)	/* ew: Collect results if more than one thread was used for the plwin_loop processing */
    if (pth.nth > 1)
      pthread_join (thread[k], &retval);



  /* ew: Run function pldlayer */
  if (pth.nr_dgates > 0)
    {
      //start[3]=gethrtime(); /* ew: record timestamt to use for processing time evaluation */
      start[3] = get_tid ();	/* ew: record timestamt to use for processing time evaluation */
      pth.out += (pth.noud * pth.res_mult);	/* ew: Move pointer to output, based on noud. Still unclear what noud is? */
      memset (pth.out + pth.nr_dgates * pth.dl_short, 0,
	      pth.nr_dgates * (pth.dl_long + 1) * sizeof (DSPCMPLX));
      for (k = 0, i = 0; i < pth.nth; i++, k = k2)
	{			/* ew: run pldlayer pointed at pth in one or many threads. some of the fftwf_plan_dft_1d defined in pth will likely be used */
	  k2 = pth.nr_dgates * (i + 1) / pth.nth;
	  pth.th = i;
	  pth.k1 = k;
	  pth.k2 = k2;
	  if (pth.nth > 1)
	    {
	      memcpy (ptth + i, &pth, sizeof (struct pth));
	      pthread_create (thread + i, &sched_glob,
			      (void *(*)(void *)) pldlayer,
			      (void *) (ptth + i));
	    }
	  else
	    pldlayer (&pth);
	}
      for (i = 0; i < pth.nth; i++)	/* ew: Collect results if more than one thread was used for the plwin_loop processing */
	if (pth.nth > 1)
	  pthread_join (thread[i], &retval);
    }


  /* ew: Free memory and cleanup */
  if (pth.debug > 1)
    printf ("%ld\n",
	    pth.noud * pth.res_mult + pth.nr_dgates * (pth.dl_short +
						       pth.dl_long + 1));
  fftwf_destroy_plan (pth.p);
  fftwf_destroy_plan (pth.p_2);
  fftwf_destroy_plan (pth.pb);
  fftwf_destroy_plan (pth.pb_2);
  if (pth.nr_dgates > 0)
    {
      fftwf_destroy_plan (pth.pd);
      fftwf_destroy_plan (pth.pdb);
      free (pth.fir_samp);
    }
  if (pth.nr_clutter > 0)
    {
      fftwf_destroy_plan (pth.pc);
      fftwf_destroy_plan (pth.pcb);
    }
  fftwf_free (pth.inf);
  free (pth.out_t);
  if (pth.nth > 1)
    {
      free (thread);
      free (ptth);
    }
  if (pth.debug > 0)
    {
      //start[4]=gethrtime();
      start[4] = get_tid ();
      //printf("fd_alt: %d threads: %.2g s",pth.nth,(start[4]-start[0])/CPU_SPEED);
      print_tid (start[4] - start[0], pth.nth, "plwin fd_alt");
      if (float_notch)
	printf (" Notch: %.2g", upar[1]);
      if (pth.debug > 1)
	for (j = 0; j < 4; j++)
	  printf (" %.2f", (start[j + 1] - start[j]) / CPU_SPEED);
      printf ("\n");
    }
  return 0;
}				/* decoder_6 */




void *
plwin_loop (struct pth *pth)
/* ew: Main function for the lag profile inversion */
{
  ulong k, win, j, i, k1, k2, l, th, j2;
  int *wind;
  DSPCMPLX *of, *of_2, *to;
  DSPCMPLXSHORT *in;
  DSPCMPLX *add, *ia;
  fftwf_complex *inf, *outf, *inf_2;
  float *out_t, *out_tf, *out0, *out0_t;

  k1 = pth->k1;
  k2 = pth->k2;
  th = pth->th;
  inf = pth->inf + th * (pth->fft_len * 2 + pth->fft_len_2);
  outf = inf + pth->fft_len;
  inf_2 = inf + pth->fft_len * 2;
  memset (inf + pth->win_len * pth->frac, 0,
	  pth->nex * sizeof (fftwf_complex));
  memset (inf_2, 0, pth->fft_len_2 * sizeof (fftwf_complex));

  of = (DSPCMPLX *) inf;
  to = (DSPCMPLX *) outf;
  of_2 = (DSPCMPLX *) inf_2;
  out0 = (float *) malloc (pth->fft_len_2 * pth->nr_ugates * sizeof (float));

  add = (DSPCMPLX *) malloc (pth->nr_virtsamp * sizeof (DSPCMPLX));	/* ew: allocate memory "add". nr_virtsamp = pth.nr_samp+pth.lower_tail+pth.upper_tail; */
  /* ew: nr_samp is samples per data window */
  memset (add, 0, pth->lower_tail * sizeof (DSPCMPLX));	/* ew: zeros padded at beginning of memory "add" */
  memset (add + pth->nr_samp + pth->lower_tail, 0, pth->upper_tail * sizeof (DSPCMPLX));	/* ew: zeros padded at end of memory "add"  */

  for (win = k1; win < k2; win++)
    {				/* ew: (fixed intendent for this for loop) */
      /* ew: k1 and k2 segments into intervals based on number of threads */
      /* ew: if only one thread k1 = 0 and k2 = nr_win; nr_win is number of windows per pulse group */

      wind = pth->windows + win * pth->win_len;	/* ew: pth.windows is pointer to start of bit-code */
      /* ew: wind is hence pointer to start of code corresponding to current number on "win" */
      memset (out0, 0, pth->fft_len_2 * pth->nr_ugates * sizeof (float));	/* ew: Initiate all of out0 to zero */
      out_t = pth->out_t + (th * pth->res_mult + win % pth->res_mult) * pth->nout;	/* ew: unclear what this is, TODO. Could be pointer to output memory set per thread? */

      for (k = win; k < pth->nr_rep; k += pth->nr_win)
	{			/* ew: k goes thought the current data window number in all of the pulse groups */

	  in = pth->in + k * pth->nr_blk;	/* ew: nr_blk is samples per data window. "in" becomes pointer to input data for current data window "k" */
	  ia = add + pth->lower_tail;	/* ew: ia is pointer to memory "add" but pointing to after the initial padded zeros */

	  for (j = 0; j < pth->nr_samp; j++)
	    {			/* ew: nr_samp is samples per data window. */
	      ia[j].re = in[j].re;	/* ew: This loop will set the middle of memory add (between initial and ending padding of zero), */
	      /* to the samples for the current data window from "k" in the for-loop */
	      ia[j].im = in[j].im;
	    }

	  if (pth->do_pp > 0)
	    {
	      out_tf =
		out_t + pth->fft_len * pth->nr_gates +
		pth->fft_len_2 * (pth->nr_gates + 1 + pth->nr_ugates);
	      for (j = 0; j < pth->nr_samp; j++)
		out_tf[j] += (ia[j].re * ia[j].re + ia[j].im * ia[j].im);
	      for (j = pth->nr_samp; j < pth->nr_pp; j++)
		out_tf[j] +=
		  ((float) in[j].re * (float) in[j].re +
		   (float) in[j].im * (float) in[j].im);
	    }

	  out_tf = out_t;
	  for (i = 0; i < pth->nr_gates; i++)
	    {
	      ia = add + i * pth->frac;
	      if (pth->frac > 1)
		for (l = 0; l < pth->win_len;)
		  if (wind[l] == -1)
		    for (j = l++ * pth->frac; j < l * pth->frac; j++)
		      {
			of[j].re = -ia[j].re;
			of[j].im = -ia[j].im;
		      }
		  else
		    for (j = l++ * pth->frac; j < l * pth->frac; j++)
		      {
			of[j].re = ia[j].re;
			of[j].im = ia[j].im;
		      }
	      else
		for (j = 0; j < pth->win_len; j++)
		  {
		    of[j].re = wind[j] * ia[j].re;
		    of[j].im = wind[j] * ia[j].im;
		  }
	      if (pth->nr_dgates > 0 && i < pth->nr_dgates)
		{
		  ia = pth->fir_samp + i * pth->fft_dlayer + k;
		  for (j = 0; j < pth->fft_len - pth->nex; j++)
		    {
		      ia->re += of[j].re;
		      ia->im += of[j].im;
		    }
		}

	      fftwf_execute_dft (pth->p, inf, outf);
	      for (j = 0; j < pth->fft_len; j++)
		out_tf[j] += (to[j].re * to[j].re + to[j].im * to[j].im);

	      out_tf += pth->fft_len;
	    }
	  l = pth->lower_tail / pth->frac;
	  out0_t = out0 + l * pth->fft_len_2;
	  for (i = l; i < pth->nr_ugates - pth->upper_tail / pth->frac; i++)
	    {
	      ia = add + i * pth->frac;
	      for (j = 0; j < pth->frac; j++)
		{
		  of_2[j].re = ia[j].re;
		  of_2[j].im = ia[j].im;
		}
	      fftwf_execute_dft (pth->p_2, inf_2, outf);
	      for (j = 0; j < pth->fft_len_2; j++)
		out0_t[j] += (to[j].re * to[j].re + to[j].im * to[j].im);
	      out0_t += pth->fft_len_2;
	    }
	}

      for (j = 0; j < pth->fft_len_2 * pth->nr_ugates; j++)
	out_tf[j] += out0[j];
      out0_t = out0;
      out_tf += (pth->fft_len_2 * pth->nr_ugates);
      j2 = (pth->nr_gates + 1) * pth->fft_len_2;

      for (i = 0; i < pth->win_len - 1; i++)
	{
	  if (wind[i] * wind[i + 1] == -1)
	    for (j = 0; j < j2; j++)
	      out_tf[j] -= out0_t[j];
	  else
	    for (j = 0; j < j2; j++)
	      out_tf[j] += out0_t[j];
	  out0_t += pth->fft_len_2;
	}
    }

  free (out0);
  free (add);
  return NULL;
}



void *
pldec_loop (struct pth *pth)
{
  ulong k, j, i, l, ustep;
  fftwf_complex *inf, *outf;
  DSPCMPLX *of, *to, *too, *add;
  float *out_tf, *out_tt;
  int skip0 = 1;

  k = pth->k1;
  inf = pth->inf + pth->th * pth->fft_len * 2;
  outf = inf + pth->fft_len;

  of = (DSPCMPLX *) inf;
  to = (DSPCMPLX *) outf;
  memset (of, 0, pth->fft_len * sizeof (DSPCMPLX));
  if (pth->nr_undec > 0)
    ustep = (pth->undec2 - pth->undec1 + 1) / pth->nr_undec;
  out_tt = pth->out_t + pth->nr_gates * pth->fft_len + k * pth->nout;
  add =
    pth->out + (pth->nr_gates + 1) * pth->frac +
    pth->nr_gates * pth->nr_lags + k * pth->noud;
  for (i = 0; i < pth->nr_undec; i++)
    {
      for (l = 0; l < ustep; l++)
	{
	  out_tf = out_tt + (i * ustep + pth->undec1 + l) * pth->fft_len_2;
	  for (j = 0; j < pth->fft_len_2; j++)
	    of[j].re += out_tf[j] / pth->fft_len_2 / ustep;
	}
      fftwf_execute_dft (pth->pb_2, inf, outf);
      for (j = 0; j < pth->frac; j++)
	{
	  add[j * pth->nr_undec].re = to[j].re;
	  add[j * pth->nr_undec].im = to[j].im;
	}
      add++;
      memset (of, 0, pth->fft_len_2 * sizeof (DSPCMPLX));
    }
  out_tt = pth->out_t + k * pth->nout;
  add = pth->out + k * pth->noud + (pth->nr_gates + 1) * pth->frac;
  too = to + skip0;
  for (i = 0; i < pth->nr_gates; i++)
    {
      out_tf = out_tt + i * pth->fft_len;
      for (j = 0; j < pth->fft_len; j++)
	of[j].re = out_tf[j] / pth->fft_len;
      fftwf_execute_dft (pth->pb, inf, outf);
      for (j = 0; j < pth->nr_lags; j++)
	{
	  add[j * pth->nr_gates].re = too[j].re;
	  add[j * pth->nr_gates].im = too[j].im;
	}
      out_tf = out_tt + i * pth->fft_len_2 + pth->nr_gates * pth->fft_len;
      for (l = 1; l < pth->win_len; l++)
	for (j = 0; j < pth->fft_len_2; j++)
	  out_tf[j] += out_tf[j + l * pth->fft_len_2];
      for (j = 0; j < pth->fft_len_2; j++)
	of[j].re = out_tf[j] / pth->fft_len_2;
      fftwf_execute_dft (pth->pb_2, inf, outf);
      for (j = 0; j < pth->frac - skip0; j++)
	{
	  add[j * pth->nr_gates].re -= too[j].re;
	  add[j * pth->nr_gates].im -= too[j].im;
	}
      add++;
    }
  add = pth->out + k * pth->noud;
  out_tf =
    out_tt + pth->nr_gates * pth->fft_len + pth->nr_ugates * pth->fft_len_2;
  for (i = 0; i < pth->nr_gates + 1; i++)
    {
      for (j = 0; j < pth->fft_len_2; j++)
	of[j].re = out_tf[j] / pth->fft_len_2;
      fftwf_execute_dft (pth->pb_2, inf, outf);
      for (j = 0; j < pth->frac; j++)
	{
	  add[j * (pth->nr_gates + 1)].re = to[j].re;
	  add[j * (pth->nr_gates + 1)].im = to[j].im;
	}
      add++;
      out_tf += pth->fft_len_2;
    }
  add =
    pth->out + (pth->nr_gates + 1) * pth->frac +
    pth->nr_gates * pth->nr_lags + pth->nr_undec * pth->frac + k * pth->noud;
  if (pth->do_pp > 0)
    for (j = 0; j < pth->nr_pp; j++)
      {
	add[j].re = out_tf[j];
	add[j].im = 0;
      }
  return NULL;
}

void *
pldlayer (struct pth *pth)
{
  ulong k, j, i;
  fftwf_complex *outf;
  DSPCMPLX *to, *too, *add, *out_t, *out_s, *samp;

  outf = pth->inf + pth->th * (pth->fft_dlayer);
  to = (DSPCMPLX *) outf;

  out_t = pth->out + pth->nr_dgates * pth->dl_short;
  out_s = out_t + pth->nr_dgates * pth->dl_long;
  for (i = pth->k1; i < pth->k2; i++)
    {
      add = pth->fir_samp + i * pth->fft_dlayer;
      samp = out_s + i;
      for (j = 0; j < pth->nr_rep; j++)
	{
	  samp->re += add[j].re;
	  samp->im += add[j].im;
	}
      fftwf_execute_dft (pth->pd, (fftwf_complex *) add, outf);
      for (j = 0; j < pth->fft_dlayer; j++)
	{
	  add[j].re =
	    (to[j].re * to[j].re + to[j].im * to[j].im) / pth->fft_dlayer;
	  add[j].im = 0;
	}
      fftwf_execute_dft (pth->pdb, (fftwf_complex *) add, outf);
      add = pth->out + i;
      for (j = 0; j < pth->dl_short; j++)
	{
	  add[j * pth->nr_dgates].re = to[j + 1].re;
	  add[j * pth->nr_dgates].im = to[j + 1].im;
	}
      add = out_t + i;
      for (k = 0; k < pth->dl_long; k++)
	{
	  too = to + (k + 1) * (pth->dl_short + 1);
	  for (j = 0; j < (pth->dl_short + 1); j++)
	    {
	      add->re += too[j].re;
	      add->im += too[j].im;
	    }
	  add += pth->nr_dgates;
	}
    }
  return NULL;
}

void *
plwin_clutter (struct pth *pth)
/* ew: Function for filtering out ground reflections.
* This is done by looking for zero doppler shift, then applying a notch filter */
{
  ulong k, j, i, win;
  int notch;
  DSPCMPLXSHORT *in;
  DSPCMPLX *clx, *cly;
  float clf, clm, *clt;

  clx = (DSPCMPLX *) pth->inf + pth->th * pth->fft_clutter * 2;	/* ew: th is increased from 0 up to nth in loop for starting serveral threads */
  /* ew: leads to different settings on clx. clx is a pointer, so the pointer that changes */
  /* ew: clx takes the adress from poitner pth.inf, and increases it with +pth.th*pth.ff_clutter */
  cly = clx + pth->fft_clutter;
  i = pth->nr_rep / pth->nr_win;	/* ew: This is the number of pulse groups */
  notch = pth->notch;

  if (notch < 0)
    {				/* ew: This might never happen based on setting of notch at start of decode6?! if par->notch is -1 pth - notch is set to upar[0], which is set to hard 0 from prface */
      /* If decode6 is run through e.g. pyface this can be entered. prface is simplified interface, upars can be changed during operation */

      clt = (float *) calloc (pth->fft_clutter, sizeof (float));	/* ew: Allocate memory for clt; cumsum of absolute value ^2 of spectrum from fft  */

      /* ew: inf (and clx defined via inf) is still just an allocated memory but without contents, needs to be filled before fft can be performed */

      for (j = pth->k1; j < pth->k2; j++)	/* ew: j goes from k1 to k2, number of points to declutter (nr_clutter), segmented into intervals based on number of threads */
	/* ew: example, one tread --> j goes from 0 to nr_clutter */

	for (win = 0; win < pth->nr_win; win++)
	  {			/* ew: win goes from 0 to nr_win. One window is one coded pulse. The loop is over windows in one pulse group */
	    for (k = win; k < pth->nr_rep; k += pth->nr_win)
	      {			/* ew: k goes from current window number into the same window in each pulse group */

		in = pth->in + k * pth->nr_blk + j;	/* ew: nr_blk = Number of samples in the power profile = samples per subcycle (aka data window, pulse) + extra samples */
		/* ew: means that "in" will go through the first nr_clutter number of points in all of the data windows in all the pulse groups */

		clx[k / pth->nr_win].re = in->re;	/* ew: [(index of data window)/(windows per pulese group)] --> [number of pulse-group] */
		/* ew: The corresponding sample from each pulse group will be placed as a series into clx */
		/* ew: This will be the input for the fft. The fft will be repeated for each of the first nr_clutter samples for all windows */
		clx[k / pth->nr_win].im = in->im;	/* ew: Imaginary part set the same way as real part above */
	      }
	    memset (clx + i, 0, (pth->fft_clutter - i) * sizeof (DSPCMPLX));	/* ew: Set zeros in clx after the entered values (?) */

	    /* ew: Execute fft */
	    fftwf_execute_dft (pth->pc, (fftwf_complex *) clx, (fftwf_complex *) cly);	/* ew: input pointer clx, output pointer cly */
	    for (k = 0; k < pth->fft_clutter; k++)
	      clt[k] += (cly[k].re * cly[k].re + cly[k].im * cly[k].im);	/* clt = cumulative sum of (re**2 + im**2), of output spectrum from fft. Increases over the for loops */
	  }

      /* ew: All fft:s has been performed and the cumsum clt has been added for each fft evaluation */
      /* ew: Next step is finding the parameters for the notch (?) */
      notch = 0;
      clf = 2. * clt[0];
      clm = 0;
      j = pth->fft_clutter - 1;
      for (k = 1; k < pth->fft_clutter; k++)
	clm += clt[k];		/* ew: clm is cumulative sum of clt[1:] */
      while (notch < -pth->notch + 1 && clf / 2. > clm / j)
	{
	  notch++;
	  clf = clt[notch] + clt[pth->fft_clutter - notch];
	  clm -= clf;
	  j -= 2;
	  clf = 2. * cl_min (clt[notch], clt[pth->fft_clutter - notch]);
	}
      notch--;
      free (clt);
      pth->notch = notch;
    }


  /* ew: Previous if statment might never happen depending on upar[0]? */
  /* ew: ->jumps directly here. If previous if-statment is executed, notch will change value and this statement can maybe still be executed if notch >=0 */
  /* ew: Hypotheis, this if-statement applys the notch filter. Last if-statement finds value for width of notch filter if program is told to do so, but seems width 0 is set when running from prface */
  /* ew: Ingemar -> hypotheis above correct */
  if (notch >= 0)
    {

      clf = 1. / sqrt (1. - (float) (2 * notch + 1) / (float) i) / (float) pth->fft_clutter;	/* ew: Ingemar -> applying notch removes energy from spectrum. This factor conserves energy in signal. Applyed in time domain */

      for (j = pth->k1; j < pth->k2; j++)	/* ew: j goes from k1 to k2, number of points to declutter (nr_clutter), segmented into intervals based on number of threads */
	/* ew: example, one tread --> j goes from 0 to nr_clutter */

	for (win = 0; win < pth->nr_win; win++)
	  {			/* ew: win goes from 0 to nr_win. One window is one coded pulse. The loop is over windows in one pulse group */

	    for (k = win; k < pth->nr_rep; k += pth->nr_win)
	      {			/* ew: k goes from current window number into the same window in each pulse group */

		in = pth->in + k * pth->nr_blk + j;	/* ew: nr_blk = Number of samples in the power profile = samples per subcycle (aka data window, pulse) + extra samples */
		/* ew: means that "in" will go through the first nr_clutter number of points in all of the data windows in all the pulse groups */

		clx[k / pth->nr_win].re = in->re;	/* ew: [(index of data window)/(windows per pulese group)] --> [number of pulse-group] */
		/* ew: The corresponding sample from each pulse group will be placed as a series into clx */
		/* ew: By taking the same sample in each group the coding is the same and can be disregarded from /Ingemar */
		/* ew: This will be the input for the fft. The fft will be repeated for each of the first nr_clutter samples for all windows */
		clx[k / pth->nr_win].im = in->im;	/* ew: Imaginary part set the same way as real part above */
	      }
	    memset (clx + i, 0, (pth->fft_clutter - i) * sizeof (DSPCMPLX));

	    /* ew: Execute fft */
	    fftwf_execute_dft (pth->pc, (fftwf_complex *) clx,
			       (fftwf_complex *) cly);

	    /* ew: Apply the notch filter onto the calculated spectrum */
	    memset (cly, 0, (notch + 1) * sizeof (DSPCMPLX));	/* ew: The first values of the spectrum are set to 0. Number of zeros depends on notch fitler width */

	    /* ew: Set zeros in cly after the length of the fft-spectrum minus length of notch (?) */
	    memset (cly + pth->fft_clutter - notch, 0,
		    notch * sizeof (DSPCMPLX));

	    /* ew: Do backwards fft from spectrum to time series on the notch filtered spectrum */
	    fftwf_execute_dft (pth->pcb, (fftwf_complex *) cly,
			       (fftwf_complex *) clx);

	    /* ew: Adjust the input data based on time series after notch filter. Input data is adjusted in svereal steps during the execution of the for-loops */
	    for (k = win; k < pth->nr_rep; k += pth->nr_win)
	      {
		in = pth->in + k * pth->nr_blk + j;
		in->re = rint (clx[k / pth->nr_win].re * clf);
		in->im = rint (clx[k / pth->nr_win].im * clf);
	      }
	  }
    }

  return NULL;
}

void
matface (int *par, int *nin, double *in_r, double *in_i,
	 int *nout, double *out_r, double *out_i, double *upar)
{
  ulong nb = par[0] + par[4] * par[7];
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
  decoder_6 (nb, par, 0, fpar, in, out, upar);
  for (i = 0; i < *nout; i++)
    {
      out_r[i] = out[i].re;
      out_i[i] = out[i].im;
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
  ulong nb = par[0] + par[4] * par[7];
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
  free (in);
}

void
prface (int *p)
{
  ulong nb = p[0] + p[4] * p[7];
  ulong nout =
    ((p[2] + 1) * p[8] + p[2] * p[1] / 2 + p[3] * p[8]) * p[9] + (p[6] +
								  p[10]) *
    p[9] + p[19] * (p[17] + p[18] + 1);
  ulong nin = (p[6] + p[10]) * p[5];
  float *fpar = NULL;
  double *upar;
  int i, *par;
  DSPCMPLXSHORT *in;
  DSPCMPLX *out;
  in = (DSPCMPLXSHORT *) malloc (nin * sizeof (DSPCMPLXSHORT));
  out = (DSPCMPLX *) malloc (nout * sizeof (DSPCMPLX));
  upar = (double *) malloc (20 * sizeof (double));
  par = (int *) malloc (nb * sizeof (int));
  upar[0] = 0;
  for (i = 0; i < p[0]; i++)
    {
      par[i] = p[i];
    }
  for (i = p[0]; i < nb; i++)
    {
      par[i] = 2 * i % 2 - 1;
    }
  if (par[15] > 0)
    printf ("No pars:%ld\nNo in :%3.0f Msamples\nNo out:%3.0f Mentries\n", nb,
	    nin / 1.e6, nout / 1.e6);
  decoder_6 (nb, par, 0, fpar, in, out, upar);
  free (in);
  free (out);
  free (par);
}

int
main (int argc, char *argv[])
{
  int i, np = 24, par[np];
  FILE *fid = NULL;
  if (argc != 2 || (fid = fopen (argv[1], "r")) == NULL)
    {
      fprintf (stderr, "Usage: plwin parfile\n");
      exit (-1);
    }
  for (i = 0; i < np; i++)
    {
      if (fscanf (fid, "%d\n", &par[i]) != 1)
	exit (1);
    }
  fclose (fid);
  prface (par);
}
