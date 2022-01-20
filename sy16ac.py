#!/usr/bin/env python3
def myexp(site):
	try:
		import lpgen.par_gen as par_gen
	except:
		import par_gen
	###### Setup section
	dspexp='sy16ac'		#Experiment version
	cal_samp=0		#Number of calibration and bakground samples
	loops=32		#Number of loops in a complete cycle
	code_len=16		#Number of Bauds in each code
	nr_codes=32		#Number of codes
	code_tx=16		#Number of Bauds to send
	nr_loop=330		#Number of loops to get wanted integration time 
	plasma_pulses=0		#Number of pulses to correlate
	plasma_frac=5		#Plasma frac
	ndgat=0
	dlong=0
	dlong=0
	nr_pulses=1		#Number of pulses to correlate
	ion_frac=3
	ion_lag=47		#Maxlag for the ion line
	mthread=20
	dshort=0
	ipp=16000*82	#IPP length in us
	baud_len=30	#Baud length in us
	start_tx=0	#When to start transmitting in us
	trx_frq=450	#Transmit frequency used
	start_samp=int(90/0.15)	#When to start sampling
	calstop=ipp
	if site=='r':
		guard=100	#Guard time in us
		start_samp=baud_len-guard	# When to start sampling in us
		nr_pulses=1	#Number of pulses to correlate
		ion_frac=1
		ion_lag=31
		tails=0
		toptail=0
		nr_fullgates=int(2*guard/baud_len+1)
		calstop=ipp-50
	else:
		site='h'
		isamp=738
		tails=code_tx-2
		lowtail=tails
		toptail=0
	samp_speed=baud_len/ion_frac
	start_samp+=start_tx	# When to start sampling in us
	print('Doing experiment files for site='+site)
	exp_name=dspexp+'_'+site#Name of experiment

	######
	nr_fullgates=isamp/ion_frac-(code_tx-1)
	if site=='r':
		clutts=0
		nr_fullgates=2*guard/baud_len+1
	else:
		clutts=4	#Clutter window us
	print(isamp,clutts)
	ac_code=par_gen.acgen(code_len,code_tx,nr_codes,'b')
	isamp=par_gen.plwingen(nr_pulses,plasma_pulses,plasma_frac,code_tx,nr_fullgates,dspexp,site,ion_frac,tails,mthread,nr_codes,nr_loop,dshort,dlong,ndgat,clutts,toptail,lowtail,loops,ac_code)
	print(isamp,clutts)

	par_gen.acdecgen(exp_name,ac_code,code_tx,nr_loop,loops,isamp,ion_frac,ion_lag)
	par_gen.t2ps(cal_samp,samp_speed,loops,baud_len,ac_code,code_len,code_tx,start_tx,ipp,trx_frq,site,dspexp,start_samp,isamp,calstop)
