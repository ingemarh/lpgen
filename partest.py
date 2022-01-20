#!/usr/bin/env python3

def pltest(parfile,doplot=0,writedata=0,data=None):
    import numpy,lp,scipy.io
    par=numpy.intc(numpy.loadtxt(parfile))
    if data:
        mat=scipy.io.loadmat(data,squeeze_me=True)
        ddata='d_data' in mat.keys()
        if 'd_raw' in mat.keys():
            nr=int(len(mat['d_raw'])/2)
            upar=mat['d_parbl'][42:61]
            draw=mat['d_raw'][:nr]
        elif 'data_raw' in mat.keys():
            upar=numpy.zeros(20)
            draw=mat['data_raw']
            #draw=draw.transpose().flatten().astype('complex128')
            draw=draw.flatten().astype('complex128')
            nr=len(draw)
        else:
            print('Unknown format')
            return
    else:
        data=parfile
        ddata=None
        nr=0
        draw=None
        upar=numpy.zeros(20)
    if 'a.par' in parfile:
        upar,dd_data,draw=lp.plwin(par,upar,draw)
        print('plwin passed')
    elif 'ac.par' in parfile:
        dd_data=lp.altdec(par,draw)
        print('altdec passed')
    
    if doplot:
        from matplotlib import pyplot as p
        ld=len(dd_data)
        nd=range(ld)
        new=dd_data.real
        if ddata:
            old=mat['d_data'][20:ld+20].real
            p.subplot(2,1,2)
            p.plot(nd,new-old)
            p.subplot(2,1,1)
        else:
            old=dd_data.imag
        p.plot(nd,new,nd,old)
        p.show()
    if writedata:
        mdic={'d_data':dd_data}
        mat=scipy.io.savemat(data+'out',mdict=mdic,do_compression=True,oned_as='column')
        print(data+'out'+' written')


if __name__ == "__main__":
    import getopt,sys,importlib
    try:
        opts,exp=getopt.getopt(sys.argv[1:],"pw")
    except getopt.GetoptError:
        sys.exit(2)
    p=w=0
    for o,a in opts:
        if o=="-p": p=1
        if o=="-w": w=1
        if o=="-h":
            print('Usage: '+sys.argv[0]+' [manda_va.par] [data_file]')
            exit(1)
    if len(exp)==1: exp.append(None)
    pltest(exp[0],p,w,exp[1])
