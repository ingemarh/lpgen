#!/usr/bin/env python3

def pltest(parfile,doplot=0,data=None):
    import numpy,lp
    par=numpy.intc(numpy.loadtxt(parfile))
    if data:
        import scipy.io
        mat=scipy.io.loadmat(data,squeeze_me=True)
        nr=int(len(mat['d_raw'])/2)
        upar=mat['d_parbl'][42:61]
        draw=mat['d_raw'][:nr]
    else:
        nr=0
        draw=None
        upar=numpy.zeros(20)
    upar,dd_data=lp.plwin(par,upar,draw)
    
    if doplot:
        from matplotlib import pyplot as p
        ld=len(dd_data)
        nd=range(ld)
        new=dd_data.real
        if nr:
            old=mat['d_data'][20:ld+20].real
            p.subplot(2,1,2)
            p.plot(nd,new-old)
            p.subplot(2,1,1)
        else:
            old=dd_data.imag
        p.plot(nd,new,nd,old)
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
    if len(exp)==1: exp.append(None)
    try:
        pltest(exp[0],p,exp[1])
    except:
        print('Usage: '+sys.argv[0]+' [manda_va.par] [data_file]')
