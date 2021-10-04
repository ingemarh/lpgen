#!/usr/bin/env python3

def plwin(par,d_parbl,d_raw):
    import numpy as np
    from numpy import ctypeslib as cl
    import ctypes
    #par[15]=3

    plc=cl.load_library("lib/plwin.so",".")
    #void pface(int *par,int nin,DSPCMPLXDBL *in_r,DSPCMPLX *out,double *upar)
    plc.pyface.argtypes=[cl.ndpointer(np.intc),ctypes.c_int,cl.ndpointer(np.complex128),
                        cl.ndpointer(np.complex64),cl.ndpointer(np.float_)]
    nout=((par[2]+1)*par[8]+par[2]*par[1]/2+par[3]*par[8])*par[9]+(par[6]+par[10])*par[9]+par[19]*(par[17]+par[18]+1)
    out=np.zeros(int(nout),np.complex64);
    plc.pyface(par,len(d_raw),d_raw,out,d_parbl[42:61])
    return d_parbl,out

def pltest(parfile,data,doplot=0):
    import numpy
    par=numpy.intc(numpy.loadtxt(parfile))
    import scipy.io
    mat=scipy.io.loadmat(data,squeeze_me=True)
    nr=int(len(mat['d_raw'])/2)
    d_parbl,dd_data=plwin(par,mat['d_parbl'],mat['d_raw'][:nr-1])
    
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
