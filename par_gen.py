#!/usr/bin/env python3
def Same_Code(Code,Start):
	loops=1
	first_code=Code[Start]
	while loops+Start<len(Code):
		if Code[loops+Start] != first_code: break
		loops=loops+1
	return loops
		
def Print_t2ps(t_to_ps,t1,t2,a,f):
	f1='%.1f'
	f2='%.1f'
	if t1%1==0: f1='%.0f'
	if t2%1==0: f2='%.0f'
	form=(f1+'\t'+f2+'\t%d\t%d\n')
	t_to_ps.write(form%(t1,t2,a,f))

def AC_Trx(Baud_Len,Code,t_start,t_to_ps,frq):
	t1=t_start
	j=0
	while j<len(Code):
		length=Same_Code(Code,j)
		i=Code[j]
		t2=t1+length*Baud_Len
		Print_t2ps(t_to_ps,t1,t2,i,frq)
		t1=t2
		j=j+length

def acgen(code_len,code_tx,nr_codes):
	ac_code=[]
	code_fil='m%d.ac'%code_len	#File with transmit codes
	random_code='+' * code_len
	if code_len==32:
		random_code='++--+-++++--+----+--+-+++---+-++'	#Code for randomisation
	elif code_len==256:
		random_code='++++---+-+-+++--+-+-+-+-----+++---+--+++--+-+++----------++-+--++-+---+---++-+++-+-+-+++--+-+-++++++-+++++++++-+++++--++++++++--+-+--+++++++---+--+++-+-+-+++++++--+++++++----++++------++--+-+--+-+++++-++-+--++---++-+-++-+--+-++---+++++-++--++--+++-+-+-+-++'
	ac_code=[]
	lines=open(code_fil,'r').readlines()
	ll=range(len(lines))
	import numpy
	numpy.random.seed(0)
	lll=numpy.random.permutation(ll)
	for l in ll:
		line=lines[l]
		i=0
		for s in line:
			if i>code_tx-1:	break
			c=(44-ord(random_code[i]))*(44-ord(s)) #Go from + to 1
			ac_code.append(c)
			i=i+1
	if len(ac_code)<code_tx*nr_codes:
		ac_code.extend(ac_code)
		print(len(ac_code))
		for i in range(nr_codes//2,nr_codes):
			for s in range(1,code_tx,2):
				ac_code[i*code_tx+s]*=-1 # Make it strong
	return ac_code

def plwingen(nr_pulses,plasma_pulses,plasma_frac,code_tx,nr_fullgates,dspexp,site,ion_frac,tails,mthread,nr_codes,nr_loop,dshort,dlong,ndgat,clutts,toptail,lowtail,loops,ac_code):
	for j in range(nr_pulses+plasma_pulses):
		nbits=code_tx
		nr_gates=nr_fullgates
		if(j<nr_pulses):
			parf=open(dspexp+'_'+site+chr(ord('a')+j)+'.par','w')
			frac=ion_frac
			fft_len=2*(code_tx-code_tx%2)*frac
			nr_undec=0
			undec1=0
			undec2=0
			do_pp=0
			codestart=-nbits*j
		else:
			parf=open(dspexp+'_'+site+'p.par','w')
			tails=0
			frac=plasma_frac
			fft_len=(code_tx+code_tx%2)*frac
			fft_len=300
			nr_undec=2
			undec1=0
			undec2=nr_fullgates-1
			do_pp=1
			codestart=-nbits
		upper_tail=tails
		lower_tail=tails
		dl_short=0
		dl_long=0
		debug=1
		maxthreads=mthread
		res_mult=1
		nr_cal=0
		nr_samp=frac*nbits+(nr_gates-1)*frac
		nr_pp=nr_samp+nr_cal
		nr_win=nr_codes
		nr_rep=nr_win*nr_loop
		bufjump=-nr_samp*nr_rep
		clutter=0
		notch=0
		nr_dgates=0
		codestep=1
		if j==0:
			isamp=nr_samp
			lower_tail=lowtail
			do_pp=1
			bufjump=0
			dl_short=dshort
			dl_long=dlong
			nr_dgates=ndgat
			clutter=clutts
			if clutter>0: notch=-1
		if j==nr_pulses-1:
			upper_tail=toptail
		if j>nr_pulses-1:
			bufjump=0
		nr_gates+=(lower_tail+upper_tail)
		undec1=undec1+lower_tail
		undec2=undec2+lower_tail
		nr_out=((nr_gates+1)*frac+nr_gates*fft_len/2+nr_undec*frac)*res_mult
		if do_pp>0:
			nr_out+=(nr_pp*res_mult)
		if nr_dgates>0:
			nr_out+=(nr_dgates*(dl_short+dl_long+1))
		npar=24
		pars=[npar,fft_len,nr_gates,nr_undec,nr_win,nr_rep,nr_samp,nbits,frac,
			res_mult,nr_cal,lower_tail,upper_tail,undec1,undec2,debug,
			maxthreads,dl_short,dl_long,nr_dgates,bufjump,do_pp,clutter,
			notch]
		if len(pars)!=npar:
			error('Npar mismatch')
		for i in range(npar):
			parf.write('%d\n'%pars[i])
		for i in range(code_tx*loops):
			ii=i+codestart
			ii%=(code_tx*loops)
			parf.write('%d\n'%ac_code[ii])
		parlen=npar+code_tx*nr_codes
		parf.close()
	return isamp

def cluttgen(exp_name,loops,nr_loop,isamp,clutts,ion_frac):
	parfile = open(exp_name + 'c.par', 'w')
	pars = [loops, nr_loop * loops, isamp, clutts, -1, 1, 1, ion_frac]
	for i in range(len(pars)):
		parfile.write('%d\n' % pars[i])
	parfile.close

def acdecgen(exp_name,ac_code,code_tx,nr_loop,loops,isamp,ion_frac,ion_lag):
	parfile = open(exp_name + 'ac.par', 'w')
	pars = [code_tx, nr_loop * loops, isamp, ion_frac, ion_lag, 1, loops]
	for i in range(len(pars)):
		parfile.write('%d\n' % pars[i])
	for i in range(code_tx * loops):
		parfile.write('%d\n' % ac_code[i])
	parfile.close

def t2ps(cal_samp,samp_speed,loops,baud_len,ac_code,code_len,code_tx,start_tx,ipp,trx_frq,site,expver,start_samp,isamp,calstop):
	cal_length=cal_samp*samp_speed
	if expver<2:
		t_to_ps=open('t_to_ps.txt.%s'%(site),'w')
	else:
		t_to_ps=open('t_to_ps.txt_%d.%s'%(expver,site),'w')
	k=0
	for j in range(1,loops+1):
		k=k+1;
		AC_Trx(baud_len,ac_code[code_len*(j-1):code_len*(j-1)+code_tx],start_tx+ipp*(k-1),t_to_ps,trx_frq)
		t1=start_samp+ipp*(k-1)
		t2=start_samp+isamp*samp_speed+ipp*(k-1)
		Print_t2ps(t_to_ps,t1,t2,2,trx_frq)
		t1=calstop-cal_length+ipp*(k-1)
		t3=calstop+ipp*(k-1)
		if(j%2):
			Print_t2ps(t_to_ps,t2,t3,1,0)
		Print_t2ps(t_to_ps,t1,t3,2,trx_frq)
	t_to_ps.close

def par_gen(dspexp,site):
	print('Doing experiment files for site='+site)
	###### Setup section
	expver=1.00		#Experiment version
	exp_name=dspexp+'-'+site#Name of experiment
	cal_samp=10		#Number of calibration and bakground samples
	loops=128		#Number of loops in a complete cycle
	code_len=64		#Number of Bauds in each code
	nr_codes=128		#Number of codes
	code_tx=64		#Number of Bauds to send
	nr_loop=25		#Number of loops to get wanted integration time 
	plasma_pulses=0		#Number of pulses to correlate
	plasma_frac=5		#Plasma frac
	nr_pulses=2		#Number of pulses to correlate
	ion_frac=2
	mthread=20
	dshort=nr_codes-1
	dlong=15
	if site=='l':
		ipp=1250	#IPP length in us
		baud_len=4	#Baud length in us
		start_tx=50	#When to start transmitting in us
		trx_frq=500.3	#Transmit frequency used
		start_samp=int(baud_len*(code_tx+1./ion_frac))+149	#When to start sampling
		start_samp-=2	#Bug in py2 version
		#plasma_ch=[1,5]	#Plasma line channels
		#plasma_pulses=1	#Number of pulses to correlate
		ndgat=150	#Number of Dregion gates
		calstop=ipp-4
	else:
		baud_len=2.4	#Baud length in us
		start_tx=73	#When to start transmitting in us
		code_tx=61	#Number of Bauds to send
		ipp=1500	#IPP length in us
		start_samp=int(baud_len*(code_tx+1./ion_frac))+124	#When to start sampling
		start_samp-=1	#Bug in py2 version
		ndgat=250	#Number of Dregion gates
		nr_pulses=1	#Number of pulses to correlate
		calstop=ipp-3
		if site=='v':
			trx_frq=5	#Transmit frequency used  		
		else:
			trx_frq=12	#Transmit frequency used  		
	if site=='r':
		guard=72-9	#Guard time in us
		start_samp=baud_len-guard	# When to start sampling in us
		cal_samp=300 	#Number of calibration and bakground samples
		nr_pulses=1	#Number of pulses to correlate
		ion_frac=1
		tails=0
		toptail=0
		trx_frq=5	#Transmit frequency used  		
		nr_fullgates=int(2*guard/baud_len+1)
	else:
		tails=code_tx-2
		toptail=tails
	samp_speed=baud_len/ion_frac
	start_samp+=start_tx	# When to start sampling in us

	######
	nr_fullgates=int((calstop-start_samp)/baud_len-(4+cal_samp)/ion_frac)-(code_tx-1)
	lowtail=tails
	if site=='l':
		clutts=min(269,(nr_fullgates+code_tx-1)*ion_frac)
	elif site=='r':
		nr_fullgates=int(2*guard/baud_len+1)
		ndgat=nr_fullgates
		clutts=0
	else:
		nr_fullgates=int((calstop-start_samp)/baud_len-(8+cal_samp)/ion_frac)-(code_tx-1)
		clutts=min(150,(nr_fullgates+code_tx-1)*ion_frac)
	print(nr_fullgates)

	ac_code=acgen(code_len,code_tx,nr_codes)
	isamp=plwingen(nr_pulses,plasma_pulses,plasma_frac,code_tx,nr_fullgates,dspexp,site,ion_frac,tails,mthread,nr_codes,nr_loop,dshort,dlong,ndgat,clutts,toptail,lowtail,loops,ac_code)
	t2ps(cal_samp,samp_speed,loops,baud_len,ac_code,code_len,code_tx,start_tx,ipp,trx_frq,site,expver,start_samp,isamp,calstop)
	
if __name__ == "__main__":
	import getopt,sys,importlib
	site='v'		#Name of site
	try:
		opts,exp=getopt.getopt(sys.argv[1:],"s:")
	except getopt.GetoptError:
		sys.exit(2)
	for o,a in opts:
		if o=="-s":
			site=a
	print(exp,site)
	exp=importlib.import_module(exp[0])
	print(exp)
	exp.myexp(site)
