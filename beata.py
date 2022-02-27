#!/usr/bin/env python3
def myexp(site):
	try:
		import lpgen.par_gen as par_gen
	except:
		import par_gen
	print('Doing experiment files for site='+site)
	###### Setup section
	dspexp='beata'		#Experiment version
	exp_name=dspexp+'_'+site#Name of experiment
	cal_samp=16		#Number of calibration and bakground samples
	loops=64		#Number of loops in a complete cycle
	code_len=32		#Number of Bauds in each code
	nr_codes=64		#Number of codes
	code_tx=32		#Number of Bauds to send
	nr_loop=15		#Number of loops to get wanted integration time 
	plasma_pulses=0		#Number of pulses to correlate
	plasma_frac=5		#Plasma frac
	plasma_frq=99		#Plasma frq
	plasma_gates=100	#No plasma gates
	nr_pulses=0		#Number of pulses to correlate
	ion_frac=2
	ion_lag=41		#Maxlag for the ion line
	mthread=20
	#dshort=nr_codes-1
	#dlong=15
	if site=='l':
		ipp=6250	#IPP length in us
		baud_len=50	#Baud length in us
		start_tx=50	#When to start transmitting in us
		trx_frq=500.3	#Transmit frequency used
		code_tx=30	#Number of Bauds to send
		start_samp=int(baud_len*(code_tx+1.5/ion_frac))+240	#When to start sampling
		plasma_ch=1	#Plasma line channels
		#plasma_pulses=1	#Number of pulses to correlate
		#ndgat=150	#Number of Dregion gates
	else:
		nr_loop=14	#Number of loops to get wanted integration time 
		baud_len=20	#Baud length in us
		start_tx=82	#When to start transmitting in us
		ipp=5580	#IPP length in us
		start_samp=int(baud_len*(code_tx+1./ion_frac))+305	#When to start sampling
		plasma_ch=3	#Plasma line channels
		#ndgat=250	#Number of Dregion gates
		if site=='v':
			trx_frq=5	#Transmit frequency used  		
		else:
			trx_frq=8	#Transmit frequency used  		
			plasma_frac=50
			plasma_gates=78+code_tx-1
	plasma_speed=baud_len/plasma_frac
	plasma_samp=plasma_gates*plasma_frac
	calstop=ipp-5
	if site=='r':
		guard=100	#Guard time in us
		start_samp=baud_len-guard	# When to start sampling in us
		cal_samp=225 	#Number of calibration and bakground samples
		nr_pulses=1	#Number of pulses to correlate
		ion_frac=1
		ion_lag=31
		tails=0
		toptail=0
		trx_frq=5	#Transmit frequency used  		
		nr_fullgates=int(2*guard/baud_len+1)
		calstop=ipp-50
	else:
		tails=code_tx-2
		toptail=tails
	samp_speed=baud_len/ion_frac
	start_samp+=start_tx	# When to start sampling in us

	######
	nr_fullgates=int((calstop-start_samp)/baud_len)-(cal_samp+4)/ion_frac-(code_tx-1)
	if site=='l':
		clutts=460	#Clutter window us
	elif site=='r':
		clutts=0
		nr_fullgates=2*guard/baud_len+1
	else:
		clutts=0
	isamp=ion_frac*(code_tx+nr_fullgates-1)
	print(isamp,clutts)

	ac_code=par_gen.acgen(code_len,code_tx,nr_codes)
	par_gen.acdecgen(exp_name,ac_code,code_tx,nr_loop,loops,isamp,ion_frac,ion_lag)
	par_gen.t2ps(cal_samp,samp_speed,loops,baud_len,ac_code,code_len,code_tx,start_tx,ipp,trx_frq,site,dspexp,start_samp,isamp,calstop,plasma_samp,plasma_speed)
