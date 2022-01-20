#!/usr/bin/env python3
def lpath():
    import os
    return os.path.realpath(__file__)

def plwin(par,upar,d_raw=None):
    import numpy as np
    from numpy import ctypeslib as cl
    import ctypes
    #par[15]=3

    plc=cl.load_library("lib/plwin.so",lpath())
    #void pface(int *par,int nin,DSPCMPLXDBL *in_r,DSPCMPLX *out,double *upar)
    plc.pyface.argtypes=[cl.ndpointer(np.intc),ctypes.c_int,cl.ndpointer(np.complex128),
                        cl.ndpointer(np.complex64),cl.ndpointer(np.float_)]
    nout=((par[2]+1)*par[8]+par[2]*par[1]/2+par[3]*par[8])*par[9]+(par[6]+par[10])*par[9]+par[19]*(par[17]+par[18]+1)
    if d_raw is None:
        nin=par[6]*par[5]
        ad=10
        d_raw=np.random.randint(-ad,ad+1,size=nin)+np.random.randint(-ad,ad+1,size=nin)*1j
    out=np.zeros(int(nout),np.complex64);
    plc.pyface(par,len(d_raw),d_raw,out,upar)
    return upar,out,d_raw

def acpoints(nsam,maxl,frac,nbit):
    n=0
    for l in range(1,maxl+1):
        n+=(nsam-l)
        if l>frac and l%frac and l<(nbit-1)*frac:
            n+=(nsam-l)
    return n

def altdec(par,d_raw=None):
    import numpy as np
    from numpy import ctypeslib as cl
    import ctypes

    plc=cl.load_library("lib/alt_decoder.so",lpath())
    #void pface(int *par,int nin,DSPCMPLXDBL *in_r,DSPCMPLX *out)
    plc.pyface.argtypes=[cl.ndpointer(np.intc),ctypes.c_int,cl.ndpointer(np.complex128),
                        cl.ndpointer(np.complex64)]
    nout=acpoints(par[2],par[4],par[3],par[0])
    if par[6]<0: nout-=(nout*par[1]//par[6])
    if d_raw is None:
        nin=par[1]*par[2]
        d_raw=np.random.randint(-9,10,size=nin)+np.random.randint(-9,10,size=nin)*1j
    out=np.zeros(int(nout),np.complex64);
    plc.pyface(par,len(d_raw),d_raw,out)
    return out

def uniclut(par,upar,d_raw=None):
    import numpy as np
    from numpy import ctypeslib as cl
    import ctypes

    plc=cl.load_library("lib/clutter.so",lpath())
    #void pface(int *par,int nin,DSPCMPLXDBL *in_r,DSPCMPLX *out)
    plc.pyface.argtypes=[cl.ndpointer(np.intc),ctypes.c_int,cl.ndpointer(np.complex128),
                        cl.ndpointer(np.complex64),cl.ndpointer(np.float_)]
    nout=par[4]*par[7]-par[7]*(par[7]-1)//2
    if par[7]==0: nout=1
    if d_raw is None:
        nin=par[1]*par[2]
        d_raw=np.random.randint(-9,10,size=nin)+np.random.randint(-9,10,size=nin)*1j
    out=np.zeros(int(nout),np.complex64);
    plc.pyface(par,len(d_raw),d_raw,out,upar)
    return upar,out,d_raw
