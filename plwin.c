/****************************************************************************
plwin.c -- decodump routine for decoding fft
*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <math.h>

typedef struct {
	float	re;
	float	im;
}	DSPCMPLX;

typedef struct {
	int16_t	re;
	int16_t	im;
}	DSPCMPLXSHORT;

#include <fftw3.h>
#define CPU_SPEED 1.e9
#include <sys/time.h>
#ifndef NANOSEC /* gcc case*/
#include <time.h>
#define NANOSEC 1000000000 
typedef	unsigned long long int hrtime_t;
unsigned long gethrtime(void)
{
	struct timespec ts;
	if (clock_gettime(CLOCK_MONOTONIC_RAW, &ts) != 0) return (-1);
	return ((ts.tv_sec * NANOSEC) + ts.tv_nsec);
}
#endif
#define cl_min(a,b)	(((a) < (b)) ? (a) : (b))
#include <pthread.h>
struct pth {
	int nth;
	ulong th;
	ulong k1;
	ulong k2;
	ulong nr_win;
	int *windows;
	ulong win_len;
	ulong nr_gates;
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

void *plwin_clutter(struct pth *);
void *plwin_loop(struct pth *);
void *pldec_loop(struct pth *);
void *pldlayer(struct pth *);

/* nbits is the number of integer values supplied by the user */
/* in is a pointer to the complex input data vector */
/* out is a pointer to complex output data vector */
int decoder_6(ulong nbits,int *par,ulong fbits,float *fpar,DSPCMPLXSHORT *in,
	DSPCMPLX *out,double *upar)
/* decoder_6 FFT version 3 fftw */
{
/* par: npar fft_len nr_gates nr_undec nr_win nr_rep nr_samp win_len frac res_mult nr_cal*/
/* the rest is windows (nbits=nfft*nr_win+npar) */

	float	*out_tf;
	ulong	i,j,k,k2,npar,n_deco,maxthreads,ninf,float_notch=0,nth_clutt;
	struct pth	pth;
	hrtime_t	start[]={0,0,0,0,0};
	start[0]=gethrtime();

	npar=par[0];
	pth.fft_len=par[1];
	pth.nr_gates=par[2];
	pth.nr_undec=par[3];
	pth.nr_win=par[4];
	pth.nr_rep=par[5];
	pth.nr_samp=par[6];
	pth.win_len=par[7];
	pth.frac=par[8];
	pth.res_mult=par[9];
	pth.nr_pp=pth.nr_samp+par[10];
	pth.lower_tail=par[11]*pth.frac;
	pth.upper_tail=par[12]*pth.frac;
	pth.undec1=par[13];
	pth.undec2=par[14];
	pth.debug=par[15];
	if(pth.debug>1)
		printf("%ld %ld\n",nbits,fbits);
	maxthreads=par[16];
	pth.dl_short=par[17];
	pth.dl_long=par[18];
	pth.nr_dgates=par[19];
	pth.nr_blk=pth.nr_pp;
	if(par[20]>0)
		pth.nr_blk+=par[20];
	else
		in+=par[20];
	if(pth.debug>1)
		printf("%ld %ld\n",pth.nr_blk,par[20]);
	pth.do_pp=par[21];
	pth.nr_clutter=par[22];
	if(par[23]>-1) pth.notch=par[23];
	else pth.notch=(int)upar[0];
	if(pth.notch<0)	float_notch=1;
	pth.fft_dlayer=2*(pth.dl_long+1)*(pth.dl_short+1);
	pth.fft_clutter=pth.nr_rep/pth.nr_win;
	pth.nex=pth.fft_len-pth.win_len*pth.frac;
	pth.nr_virtsamp=pth.nr_samp+pth.lower_tail+pth.upper_tail;
	pth.nr_ugates=pth.nr_virtsamp/pth.frac;
	pth.nr_lags=(pth.fft_len+1)/2;
/* Dec ffts done as a proper ACF (nfft>=2*sample)*/
	/*if(nex>0) pth.fft_len_2=2*pth.fft_len/pth.win_len;
	else*/
		pth.fft_len_2=2*pth.frac;
	pth.nout=pth.fft_len*pth.nr_gates;	/* First high res spectra*/
	pth.nout+=pth.fft_len_2*pth.nr_ugates;	/* Then low res undecoded ones*/
	pth.nout+=pth.fft_len_2*(pth.nr_gates+1);/* And decoded low res*/
	if(pth.do_pp>0) pth.nout+=pth.nr_pp;	/* Last a power profile*/
	if(pth.debug>1)
		printf("%ld %d %ld %d\n",pth.nout,pth.nex,pth.frac,pth.fft_len_2);

	pthread_attr_t	sched_glob;
	void	*retval;
	pthread_t	*thread=NULL;
	struct pth	*ptth=NULL;
	pth.nth=sysconf(_SC_NPROCESSORS_ONLN);
	if(pth.nth>maxthreads) pth.nth=maxthreads;
	if(pth.nth>pth.nr_win) pth.nth=pth.nr_win;
	if(pth.nth>1) {
		pthread_attr_init(&sched_glob);
		pthread_attr_setscope(&sched_glob,PTHREAD_SCOPE_SYSTEM);
		thread=(pthread_t *)malloc(pth.nth*sizeof(pthread_t));
		ptth=(struct pth *)malloc(pth.nth*sizeof(struct pth));
	}

	if(pth.nr_undec>pth.nr_ugates) pth.nr_undec=pth.nr_ugates;
	pth.noud=(pth.nr_gates+1)*pth.frac;	/* First the short lags*/
	pth.noud+=pth.nr_gates*pth.nr_lags;	/* Then the long lags*/
	pth.noud+=pth.nr_undec*pth.frac;	/* Finally some undecoded ones*/
	if(pth.do_pp>0) pth.noud+=pth.nr_pp;	/* Last a power profile*/
	if(pth.debug>1)
		printf("%ld %d\n",pth.noud,pth.nth);
	n_deco=pth.res_mult*pth.nout;
	pth.out_t=(float *)calloc(pth.nth*n_deco,sizeof(float));

	pth.windows=par+npar;
	pth.in=in;
	pth.out=out;
	ninf=2*pth.fft_len+pth.fft_len_2;
	if(pth.nr_dgates>0) {
		pth.fir_samp=calloc(pth.fft_dlayer*pth.nr_dgates,sizeof(DSPCMPLX));
		if(pth.fft_dlayer>ninf) ninf=pth.fft_dlayer;
		if(2*pth.fft_dlayer>pth.nth*ninf) ninf=2*pth.fft_dlayer;
	}
	if(pth.nr_clutter>0 && 2*pth.fft_clutter>pth.nth*ninf)
		ninf=2*pth.fft_clutter;
	pth.inf=fftwf_malloc(pth.nth*ninf*sizeof(fftwf_complex));
	pth.p=fftwf_plan_dft_1d(pth.fft_len,pth.inf,pth.inf+pth.fft_len,FFTW_FORWARD,FFTW_ESTIMATE);
	pth.p_2=fftwf_plan_dft_1d(pth.fft_len_2,pth.inf+pth.fft_len*2,pth.inf+pth.fft_len,FFTW_FORWARD,FFTW_ESTIMATE);
	pth.pb=fftwf_plan_dft_1d(pth.fft_len,pth.inf,pth.inf+pth.fft_len,FFTW_BACKWARD,FFTW_ESTIMATE);
	pth.pb_2=fftwf_plan_dft_1d(pth.fft_len_2,pth.inf,pth.inf+pth.fft_len,FFTW_BACKWARD,FFTW_ESTIMATE);
	if(pth.nr_dgates>0) {
		pth.pd=fftwf_plan_dft_1d(pth.fft_dlayer,pth.inf,pth.inf+pth.fft_dlayer,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		pth.pdb=fftwf_plan_dft_1d(pth.fft_dlayer,pth.inf,pth.inf+pth.fft_dlayer,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
	}
	if(pth.nr_clutter>0) {
		if(float_notch) nth_clutt=1;
		else nth_clutt=pth.nth;
		pth.pc=fftwf_plan_dft_1d(pth.fft_clutter,pth.inf,pth.inf+pth.fft_clutter,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		pth.pcb=fftwf_plan_dft_1d(pth.fft_clutter,pth.inf,pth.inf+pth.fft_clutter,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		for(k=0,i=0;i<nth_clutt;i++,k=k2) {
			k2=pth.nr_clutter*(i+1)/nth_clutt;
			pth.th=i; pth.k1=k; pth.k2=k2;
			if(nth_clutt>1) {
				memcpy(ptth+i,&pth,sizeof(struct pth));
				pthread_create(thread+i,&sched_glob,(void *(*)(void *)) plwin_clutter,(void *)(ptth+i));
			} else plwin_clutter(&pth);
		}
		if(nth_clutt>1)
			for(i=0;i<nth_clutt;i++) {
				pthread_join(thread[i],&retval);
				if(i==0) pth.notch=ptth[i].notch;
				else pth.notch+=ptth[i].notch;
				if(float_notch && pth.debug>0) printf("%d",ptth[i].notch);
			}
		upar[1]=(float)pth.notch/(float)nth_clutt;
	}

	start[1]=gethrtime();
	for(k=0,i=0;i<pth.nth;i++,k=k2) {
		k2=pth.nr_win*(i+1)/pth.nth;
		pth.th=i; pth.k1=k; pth.k2=k2;
		if(pth.nth>1) {
			memcpy(ptth+i,&pth,sizeof(struct pth));
			pthread_create(thread+i,&sched_glob,(void *(*)(void *)) plwin_loop,(void *)(ptth+i));
		} else plwin_loop(&pth);
	}

	for(i=0;i<pth.nth;i++) {
		if(pth.nth>1) {
			pthread_join(thread[i],&retval);
		}
		if(i>0) {
			out_tf=pth.out_t+i*n_deco;
			for(j=0;j<n_deco;j++) pth.out_t[j]+=out_tf[j];
	}	}
	start[2]=gethrtime();
	for(k=0;k<pth.res_mult;k++) {
		pth.k1=k; pth.th=0;
		if(pth.nth>1 && pth.nth>=pth.res_mult) {
			pth.th=k;
			memcpy(ptth+k,&pth,sizeof(struct pth));
			pthread_create(thread+k,&sched_glob,(void *(*)(void *)) pldec_loop,(void *)(ptth+k));
		} else pldec_loop(&pth);
	}

	for(k=0;k<pth.res_mult;k++)
		if(pth.nth>1) pthread_join(thread[k],&retval);
	if(pth.nr_dgates>0) {
		start[3]=gethrtime();
		pth.out+=(pth.noud*pth.res_mult);
		memset(pth.out+pth.nr_dgates*pth.dl_short,0,pth.nr_dgates*(pth.dl_long+1)*sizeof(DSPCMPLX));
		for(k=0,i=0;i<pth.nth;i++,k=k2) {
			k2=pth.nr_dgates*(i+1)/pth.nth;
			pth.th=i; pth.k1=k; pth.k2=k2;
			if(pth.nth>1) {
				memcpy(ptth+i,&pth,sizeof(struct pth));
				pthread_create(thread+i,&sched_glob,(void *(*)(void *)) pldlayer,(void *)(ptth+i));
			} else pldlayer(&pth);
		}
		for(i=0;i<pth.nth;i++)
			if(pth.nth>1) pthread_join(thread[i],&retval);
	}
	if(pth.debug>1)
		printf("%ld %d %ld %d\n",pth.nout,pth.nex,pth.frac,pth.fft_len_2);
	fftwf_destroy_plan(pth.p);
	fftwf_destroy_plan(pth.p_2);
	fftwf_destroy_plan(pth.pb);
	fftwf_destroy_plan(pth.pb_2);
	if(pth.nr_dgates>0) {
		fftwf_destroy_plan(pth.pd);
		fftwf_destroy_plan(pth.pdb);
		free(pth.fir_samp);
	}
	if(pth.nr_clutter>0) {
		fftwf_destroy_plan(pth.pc);
		fftwf_destroy_plan(pth.pcb);
	}
	fftwf_free(pth.inf);
	free(pth.out_t);
	if(pth.nth>1) {
		free(thread); free(ptth);
	}
	if(pth.debug>0) {
		start[4]=gethrtime();
		printf("fd_alt: %d threads: %.2g s",pth.nth,(start[4]-start[0])/CPU_SPEED);
		if(float_notch)
			printf(" Notch: %.2g",upar[1]);
		if(pth.debug>1)
			for(j=0;j<4;j++)
				printf(" %.2f",(start[j+1]-start[j])/CPU_SPEED);
		printf("\n");
	}
	return 0;
}	/* decoder_6 */


void *plwin_loop(struct pth *pth)
{
	ulong	k,win,j,i,k1,k2,l,th,j2;
	int	*wind;
	DSPCMPLX	*of,*of_2,*to;
	DSPCMPLXSHORT	*in;
	DSPCMPLX	*add,*ia;
	fftwf_complex	*inf,*outf,*inf_2;
	float	*out_t,*out_tf,*out0,*out0_t;

	k1=pth->k1; k2=pth->k2; th=pth->th;
	inf=pth->inf+th*(pth->fft_len*2+pth->fft_len_2);
	outf=inf+pth->fft_len;
	inf_2=inf+pth->fft_len*2;
	memset(inf+pth->win_len*pth->frac,0,pth->nex*sizeof(fftwf_complex));
	memset(inf_2,0,pth->fft_len_2*sizeof(fftwf_complex));

	of=(DSPCMPLX *)inf;
	to=(DSPCMPLX *)outf;
	of_2=(DSPCMPLX *)inf_2;
	out0=(float *)malloc(pth->fft_len_2*pth->nr_ugates*sizeof(float));
	add=(DSPCMPLX *)malloc(pth->nr_virtsamp*sizeof(DSPCMPLX));
	memset(add,0,pth->lower_tail*sizeof(DSPCMPLX));
	memset(add+pth->nr_samp+pth->lower_tail,0,pth->upper_tail*sizeof(DSPCMPLX));

for(win=k1;win<k2;win++) {
	wind=pth->windows+win*pth->win_len;
	memset(out0,0,pth->fft_len_2*pth->nr_ugates*sizeof(float));
	out_t=pth->out_t+(th*pth->res_mult+win%pth->res_mult)*pth->nout;

	for(k=win;k<pth->nr_rep;k+=pth->nr_win) {
		in=pth->in+k*pth->nr_blk;
		ia=add+pth->lower_tail;
		for(j=0;j<pth->nr_samp;j++) {
			ia[j].re=in[j].re;
			ia[j].im=in[j].im;
		}
		if(pth->do_pp>0) {
			out_tf=out_t+pth->fft_len*pth->nr_gates+pth->fft_len_2*(pth->nr_gates+1+pth->nr_ugates);
			for(j=0;j<pth->nr_samp;j++)
				out_tf[j]+=(ia[j].re*ia[j].re+ia[j].im*ia[j].im);
			for(j=pth->nr_samp;j<pth->nr_pp;j++)
				out_tf[j]+=((float)in[j].re*(float)in[j].re+(float)in[j].im*(float)in[j].im);
		}
		out_tf=out_t;
		for(i=0;i<pth->nr_gates;i++) {
			ia=add+i*pth->frac;
			if(pth->frac>1)
				for(l=0;l<pth->win_len;)
					if(wind[l]==-1)
						for(j=l++*pth->frac;j<l*pth->frac;j++) {
							of[j].re=-ia[j].re;
							of[j].im=-ia[j].im;
						}
					else
						for(j=l++*pth->frac;j<l*pth->frac;j++) {
							of[j].re=ia[j].re;
							of[j].im=ia[j].im;
						}
			else
				for(j=0;j<pth->win_len;j++) {
					of[j].re=wind[j]*ia[j].re;
					of[j].im=wind[j]*ia[j].im;
				}
			if(pth->nr_dgates>0 && i<pth->nr_dgates) {
				ia=pth->fir_samp+i*pth->fft_dlayer+k;
				for(j=0;j<pth->fft_len-pth->nex;j++) {
					ia->re+=of[j].re;
					ia->im+=of[j].im;
			}	}
			fftwf_execute_dft(pth->p,inf,outf);
			for(j=0;j<pth->fft_len;j++)
				out_tf[j]+=(to[j].re*to[j].re+to[j].im*to[j].im);
			out_tf+=pth->fft_len;
		}
		l=pth->lower_tail/pth->frac;
		out0_t=out0+l*pth->fft_len_2;
		for(i=l;i<pth->nr_ugates-pth->upper_tail/pth->frac;i++) {
			ia=add+i*pth->frac;
			for(j=0;j<pth->frac;j++) {
				of_2[j].re=ia[j].re;
				of_2[j].im=ia[j].im;
			}
			fftwf_execute_dft(pth->p_2,inf_2,outf);
			for(j=0;j<pth->fft_len_2;j++)
				out0_t[j]+=(to[j].re*to[j].re+to[j].im*to[j].im);
			out0_t+=pth->fft_len_2;
	}	}
	for(j=0;j<pth->fft_len_2*pth->nr_ugates;j++)
		out_tf[j]+=out0[j];
	out0_t=out0;
	out_tf+=(pth->fft_len_2*pth->nr_ugates);
	j2=(pth->nr_gates+1)*pth->fft_len_2;
	for(i=0;i<pth->win_len-1;i++) {
		if(wind[i]*wind[i+1]==-1)
			for(j=0;j<j2;j++) out_tf[j]-=out0_t[j];
		else
			for(j=0;j<j2;j++) out_tf[j]+=out0_t[j];
		out0_t+=pth->fft_len_2;
}	}

	free(out0);
	free(add);
	return NULL;
}

void *pldec_loop(struct pth *pth)
{
	ulong	k,j,i,l,ustep;
	fftwf_complex	*inf,*outf;
	DSPCMPLX	*of,*to,*too,*add;
	float	*out_tf,*out_tt;
	int	skip0=1;

	k=pth->k1;
	inf=pth->inf+pth->th*pth->fft_len*2;
	outf=inf+pth->fft_len;

	of=(DSPCMPLX *)inf;
	to=(DSPCMPLX *)outf;
	memset(of,0,pth->fft_len*sizeof(DSPCMPLX));
	if(pth->nr_undec>0) ustep=(pth->undec2-pth->undec1+1)/pth->nr_undec;
	out_tt=pth->out_t+pth->nr_gates*pth->fft_len+k*pth->nout;
	add=pth->out+(pth->nr_gates+1)*pth->frac+pth->nr_gates*pth->nr_lags+k*pth->noud;
	for(i=0;i<pth->nr_undec;i++) {
		for(l=0;l<ustep;l++) {
			out_tf=out_tt+(i*ustep+pth->undec1+l)*pth->fft_len_2;
			for(j=0;j<pth->fft_len_2;j++)
				of[j].re+=out_tf[j]/pth->fft_len_2/ustep;
		}
		fftwf_execute_dft(pth->pb_2,inf,outf);
		for(j=0;j<pth->frac;j++) {
			add[j*pth->nr_undec].re=to[j].re;
			add[j*pth->nr_undec].im=to[j].im;
		}
		add++;
		memset(of,0,pth->fft_len_2*sizeof(DSPCMPLX));
	}
	out_tt=pth->out_t+k*pth->nout;
	add=pth->out+k*pth->noud+(pth->nr_gates+1)*pth->frac;
	too=to+skip0;
	for(i=0;i<pth->nr_gates;i++) {
		out_tf=out_tt+i*pth->fft_len;
		for(j=0;j<pth->fft_len;j++)
			of[j].re=out_tf[j]/pth->fft_len;
		fftwf_execute_dft(pth->pb,inf,outf);
		for(j=0;j<pth->nr_lags;j++) {
			add[j*pth->nr_gates].re=too[j].re;
			add[j*pth->nr_gates].im=too[j].im;
		}
		out_tf=out_tt+i*pth->fft_len_2+pth->nr_gates*pth->fft_len;
		for(l=1;l<pth->win_len;l++)
			for(j=0;j<pth->fft_len_2;j++)
				out_tf[j]+=out_tf[j+l*pth->fft_len_2];
		for(j=0;j<pth->fft_len_2;j++)
			of[j].re=out_tf[j]/pth->fft_len_2;
		fftwf_execute_dft(pth->pb_2,inf,outf);
		for(j=0;j<pth->frac-skip0;j++) {
			add[j*pth->nr_gates].re-=too[j].re;
			add[j*pth->nr_gates].im-=too[j].im;
		}
		add++;
	}
	add=pth->out+k*pth->noud;
	out_tf=out_tt+pth->nr_gates*pth->fft_len+pth->nr_ugates*pth->fft_len_2;
	for(i=0;i<pth->nr_gates+1;i++) {
		for(j=0;j<pth->fft_len_2;j++)
			of[j].re=out_tf[j]/pth->fft_len_2;
		fftwf_execute_dft(pth->pb_2,inf,outf);
		for(j=0;j<pth->frac;j++) {
			add[j*(pth->nr_gates+1)].re=to[j].re;
			add[j*(pth->nr_gates+1)].im=to[j].im;
		}
		add++;
		out_tf+=pth->fft_len_2;
	}
	add=pth->out+(pth->nr_gates+1)*pth->frac+pth->nr_gates*pth->nr_lags+pth->nr_undec*pth->frac+k*pth->noud;
	if(pth->do_pp>0)
		for(j=0;j<pth->nr_pp;j++) {
			add[j].re=out_tf[j];
			add[j].im=0;
		}
	return NULL;
}

void *pldlayer(struct pth *pth)
{
	ulong	k,j,i;
	fftwf_complex	*outf;
	DSPCMPLX	*to,*too,*add,*out_t,*out_s,*samp;

	outf=pth->inf+pth->th*(pth->fft_dlayer);
	to=(DSPCMPLX *)outf;

	out_t=pth->out+pth->nr_dgates*pth->dl_short;
	out_s=out_t+pth->nr_dgates*pth->dl_long;
	for(i=pth->k1;i<pth->k2;i++) {
		add=pth->fir_samp+i*pth->fft_dlayer;
		samp=out_s+i;
		for(j=0;j<pth->nr_rep;j++) {
			samp->re+=add[j].re;
			samp->im+=add[j].im;
		}
		fftwf_execute_dft(pth->pd,(fftwf_complex *)add,outf);
		for(j=0;j<pth->fft_dlayer;j++) {
			add[j].re=(to[j].re*to[j].re+to[j].im*to[j].im)/pth->fft_dlayer;
			add[j].im=0;
		}
		fftwf_execute_dft(pth->pdb,(fftwf_complex *)add,outf);
		add=pth->out+i;
		for(j=0;j<pth->dl_short;j++) {
			add[j*pth->nr_dgates].re=to[j+1].re;
			add[j*pth->nr_dgates].im=to[j+1].im;
		}
		add=out_t+i;
		for(k=0;k<pth->dl_long;k++) {
			too=to+(k+1)*(pth->dl_short+1);
			for(j=0;j<(pth->dl_short+1);j++) {
				add->re+=too[j].re;
				add->im+=too[j].im;
			}
			add+=pth->nr_dgates;
	}	}
	return NULL;
}

void *plwin_clutter(struct pth *pth)
{
	ulong	k,j,i,win;
	int	notch;
	DSPCMPLXSHORT	*in;
	DSPCMPLX	*clx,*cly;
	float	clf,clm,*clt;

	clx=(DSPCMPLX *)pth->inf+pth->th*pth->fft_clutter*2;
	cly=clx+pth->fft_clutter;
	i=pth->nr_rep/pth->nr_win;
	notch=pth->notch;
	if(notch<0) {
		clt=(float *)calloc(pth->fft_clutter,sizeof(float));
		for(j=pth->k1;j<pth->k2;j++)
			for(win=0;win<pth->nr_win;win++) {
				for(k=win;k<pth->nr_rep;k+=pth->nr_win) {
					in=pth->in+k*pth->nr_blk+j;
					clx[k/pth->nr_win].re=in->re;
					clx[k/pth->nr_win].im=in->im;
				}
				memset(clx+i,0,(pth->fft_clutter-i)*sizeof(DSPCMPLX));
				fftwf_execute_dft(pth->pc,(fftwf_complex *)clx,(fftwf_complex *)cly);
				for(k=0;k<pth->fft_clutter;k++)
					clt[k]+=(cly[k].re*cly[k].re+cly[k].im*cly[k].im);
			}
		notch=0; clf=2.*clt[0]; clm=0; j=pth->fft_clutter-1;
		for(k=1;k<pth->fft_clutter;k++) clm+=clt[k];
		while(notch<-pth->notch+1 && clf/2.>clm/j) {
			notch++;
			clf=clt[notch]+clt[pth->fft_clutter-notch];
			clm-=clf; j-=2;
			clf=2.*cl_min(clt[notch],clt[pth->fft_clutter-notch]);
		}
		notch--;
		free(clt);
		pth->notch=notch;
	}
	if(notch>=0) {
		clf=1./sqrt(1.-(float)(2*notch+1)/(float)i)/(float)pth->fft_clutter;
		for(j=pth->k1;j<pth->k2;j++)
			for(win=0;win<pth->nr_win;win++) {
				for(k=win;k<pth->nr_rep;k+=pth->nr_win) {
					in=pth->in+k*pth->nr_blk+j;
					clx[k/pth->nr_win].re=in->re;
					clx[k/pth->nr_win].im=in->im;
				}
				memset(clx+i,0,(pth->fft_clutter-i)*sizeof(DSPCMPLX));
				fftwf_execute_dft(pth->pc,(fftwf_complex *)clx,(fftwf_complex *)cly);
				memset(cly,0,(notch+1)*sizeof(DSPCMPLX));
				memset(cly+pth->fft_clutter-notch,0,notch*sizeof(DSPCMPLX));
				fftwf_execute_dft(pth->pcb,(fftwf_complex *)cly,(fftwf_complex *)clx);
				for(k=win;k<pth->nr_rep;k+=pth->nr_win) {
					in=pth->in+k*pth->nr_blk+j;
					in->re=rint(clx[k/pth->nr_win].re*clf);
					in->im=rint(clx[k/pth->nr_win].im*clf);
	}		}	}
	
	return NULL;
}

void matface(int *par,int *nin,double *in_r,double *in_i,
	int *nout, double *out_r, double *out_i, double *upar) {
	ulong nb=par[0]+par[4]*par[7];
	float *fpar=NULL;
	int i;
	DSPCMPLXSHORT *in;
	DSPCMPLX *out;
	in=(DSPCMPLXSHORT *)malloc(*nin*sizeof(DSPCMPLXSHORT));
	for(i=0;i<*nin;i++) {
		in[i].re=in_r[i];
		in[i].im=in_i[i];
	}
	out=(DSPCMPLX *)malloc(*nout*sizeof(DSPCMPLX));
	decoder_6(nb,par,0,fpar,in,out,upar);
	for(i=0;i<*nout;i++) {
		out_r[i]=out[i].re;
		out_i[i]=out[i].im;
	}
	free(in); free(out);
}
