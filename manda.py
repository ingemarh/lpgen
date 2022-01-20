#!/usr/bin/env python3
def myexp(site):
	try:
		import lpgen.par_gen as par_gen
	except:
		import par_gen
	print('Doing experiment files for site='+site)
	###### Setup section
	dspexp='manda'		#Experiment version
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

	ac_code=par_gen.acgen(code_len,code_tx,nr_codes)
	isamp=par_gen.plwingen(nr_pulses,plasma_pulses,plasma_frac,code_tx,nr_fullgates,dspexp,site,ion_frac,tails,mthread,nr_codes,nr_loop,dshort,dlong,ndgat,clutts,toptail,lowtail,loops,ac_code)
	par_gen.t2ps(cal_samp,samp_speed,loops,baud_len,ac_code,code_len,code_tx,start_tx,ipp,trx_frq,site,dspexp,start_samp,isamp,calstop)
