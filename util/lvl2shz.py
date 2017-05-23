#!/usr/bin/env python2.7

"""
Usage: 
    lvl2shz.py [--with-strain] ( [-] | <file.lvl> )

Options:
    --with-strain    interpolates a strain distribution

Description:
    lvl2shz.py converts a level definition to finely sampled strain volume file
		
    lvl2shz subsamples a strain volume into smaller cuboid volumes
    of length, thickness, and width starting from lo, to, and wo, respectively
    increasing geometrically with down-dip distance with increment 
    alphal, alphat, and alphaw (alphal>1 for increase).

    input segment file is a list of:
      xo     origin position vector [x1 (north);x2 (east);x3 (down)]
      L      total length
      W      total width
      T      total thickness
      strike strike angle in degrees
      dip    dip angle in degrees
      lo     approximative initial length of output segments
      wo     approximative initial width of output segments
      to     approximative initial thickness of output segments
      alphal geometric factor for length increase
      alphaw geometric factor for width increase
      alphat geometric factor for thickness increase

    output list of output segments in the format
    i,x1,x2,x3,length,width,thickness,strike,dip

@author: sbarbot
"""

from docopt import docopt
from numpy import append,array,pi,cos,sin,ceil,savetxt
from sys import argv,exit,stdin,stdout

def lvl2shz(index,x1o,x2o,x3o,L,W,T,strike,dip,lo,wo,to,alphal,alphaw,alphat,strain=None):
    
    d2r=pi/180
    # create wi
    Wc=W
    k=0
    w=array([0])
    while Wc>0:
        Wt=wo*alphaw**k
        if Wt > Wc/2:
            Wt = Wc
    
        wn=min(Wt,Wc)
        w=append(w, wn)
        k=k+1
        Wc=Wc-wn
        
    Nw=k

    # create ti
    Tc=T
    k=0
    t=array([0])
    while Tc>0:
        Tt=to*alphat**k
        if Tt > Tc/2:
            Tt = Tc
    
        tn=min(Tt,Tc)
        t=append(t, tn)
        k=k+1
        Tc=Tc-tn
        
    Nt=k

    # strike, dip, and normal direction vectors
    Sv=[ cos(strike*d2r), sin(strike*d2r), 0]
    Dv=[-cos(dip*d2r)*sin(strike*d2r), cos(dip*d2r)*cos(strike*d2r), sin(dip*d2r)]
    Nv=[-sin(dip*d2r)*sin(strike*d2r), sin(dip*d2r)*cos(strike*d2r),-cos(dip*d2r)]

    # loop in dip direction
    for k in range(Nw):
        lt=lo*alphal**k
        Nl=int(ceil(L/lt))
        lt=float(L)/float(Nl)

        tt=to*alphat**k
        Nt=int(ceil(T/tt))
        tt=float(T)/float(Nt)
        
        # loop in normal direction
        for j in range(Nt):
            # loop in strike direction
            for i in range(Nl):
                x1=x1o+i*lt*Sv[0]+((float(j)+0.5)*tt-0.5*T)*Nv[0]+sum(w.take(range(k+1)))*Dv[0]
                x2=x2o+i*lt*Sv[1]+((float(j)+0.5)*tt-0.5*T)*Nv[1]+sum(w.take(range(k+1)))*Dv[1]
                x3=x3o+i*lt*Sv[2]+((float(j)+0.5)*tt-0.5*T)*Nv[2]+sum(w.take(range(k+1)))*Dv[2]
                index=index+1

                if strain[0] is None:
                    savetxt(stdout,[[index,x1,x2,x3,lt,w[k+1],t[j+1],strike,dip]],delimiter=" ",fmt="%4i %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %5.2f %5.2f")
                else:
                    savetxt(stdout,[[index,strain,x1,x2,x3,lt,w[k+1],t[k+1],strike,dip]],delimiter=" ",fmt="%4i %+10.3e %+10.3e %+10.3e %10.3e %10.3e %10.3e %8.4f %8.4f %8.4f %8.3f %8.3f")

    return index
    
def main(): 
	arguments = docopt(__doc__, version='1.0')
	offset=0
	isWithStrain=arguments['--with-strain']
	if isWithStrain: 
		offset+=6
	
	if 1+offset==len(argv):
		fid=stdin
	else:
		fname=argv[1+offset]
		#print fname, len(argv)
		#if not path.isfile(fname):
		#    raise ValueError("invalid file name: " + fname)
		fid=open(fname, 'r')
 
	if isWithStrain:
		print '#  n       e11     e12    e13    e22    e23     e33        x1        x2        x3    length     width   thickness     strike   dip'
	else:
		print '#  n        x1        x2        x3    length     width thickness    strike   dip'
        
	k=0
        for line in iter(fid.readlines()):
            if '#'==line[0]:
                continue
            numbers=map(float, line.split())
            e11=None
            e12=None
            e13=None
            e22=None
            e23=None
            e33=None
            if 15==len(numbers):
                i,x1,x2,x3,length,width,thickness,strike,dip,Lo,Wo,To,al,aw,at=numbers
            elif 21==len(numbers):
                if isWithStrain:
                    i,e11,e12,e13,e22,e23,e33,x1,x2,x3,length,width,thickness,strike,dip,Lo,Wo,To,al,aw,at=numbers
                else:
                    i,e11,e12,e13,e22,e23,e33,x1,x2,x3,length,width,thickness,strike,dip,Lo,Wo,To,al,aw,at=numbers
            else:
                ValueError("invalid number of column in input file "+fname)

            k=lvl2shz(k,x1,x2,x3,length,width,thickness,strike,dip,Lo,Wo,To,al,aw,at,strain=[e11,e12,e13,e22,e23,e33])

if __name__ == "__main__":
	main()
