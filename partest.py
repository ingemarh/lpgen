#!/usr/bin/env python3

def pltest(parfile,data,doplot=0):
    import numpy,lp
    par=numpy.intc(numpy.loadtxt(parfile))
    import scipy.io
    mat=scipy.io.loadmat(data,squeeze_me=True)
    nr=int(len(mat['d_raw'])/2)
    upar,dd_data=lp.plwin(par,mat['d_parbl'][42:61],mat['d_raw'][:nr-1])
    
    if doplot:
        from matplotlib import pyplot as p
        ld=len(dd_data)
        nd=range(ld)
        new=dd_data.real
        old=mat['d_data'][20:ld+20].real
        p.subplot(2,1,1)
        p.plot(nd,new,nd,old)
        p.subplot(2,1,2)
        p.plot(nd,new-old)
        p.show()

if __name__ == "__main__":
    import getopt,sys,importlib
    try:
        opts,exp=getopt.getopt(sys.argv[1:],"p")
    except getopt.GetoptError:
        sys.exit(2)
    p=0
    for o,a in opts:
        if o=="-s": site=a
        if o=="-p": p=1
    try:
        pltest(exp[0],exp[1],p)
    except:
        print('Usage: '+sys.argv[0]+' [manda_va.par] [data_file]')
