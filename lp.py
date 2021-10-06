#!/usr/bin/env python3

def plwin(par,upar,d_raw):
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
    plc.pyface(par,len(d_raw),d_raw,out,upar)
    return upar,out
