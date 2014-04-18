#!/usr/bin/env python

"""
Usage:
    seg2flt.py [--with-slip] ( [-] | <file.seg>)

Option:
    --with-slip interpolates a slip distribution

Description:
    seg2flt.py converts a segment definition to finely sampled fault file
		
    seg2flt subsamples a fault patch into smaller segments
    of length and width starting from lo and wo, respectively
    increasing geometrically with down-dip distance with increment 
    alphal and alphaw (alphal>1 for increase).

    input segment file is a list of:
      xo     origin position vector [x1 (north);x2 (east);x3 (down)]
      L      total length
      W      total width
      strike strike angle in degrees
      dip    dip angle in degrees
      rake   rake angle of slip
      lo     approximative initial length of output segments
      wo     approximative initial width of output segments
      alpha1 geometric factor for length increase
      alphaw geometric factor for width increase

    output list of output segments in the format
    i,x1,x2,x3,length,width,strike,dip,rake

@author: sbarbot
"""

from docopt import docopt
from numpy import append,array,pi,cos,sin,ceil,savetxt
from sys import argv,exit,stdin,stdout

def seg2flt(index,x1o,x2o,x3o,L,W,strike,dip,rake,lo,wo,alphal,alphaw,slip=None):
    
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

    # strike and dip direction normal vectors
    Sv=[ cos(strike*d2r), sin(strike*d2r), 0]
    Dv=[-cos(dip*d2r)*sin(strike*d2r), cos(dip*d2r)*cos(strike*d2r), sin(dip*d2r)]
    
    # loop in dip direction
    for j in range(Nw):
        lt=lo*alphal**j
        Nl=int(ceil(L/lt))

        lt=L/Nl
        
        # loop in strike direction
        for i in range(Nl):
            x1=x1o+i*lt*Sv[0]+sum(w.take(range(j+1)))*Dv[0]
            x2=x2o+i*lt*Sv[1]+sum(w.take(range(j+1)))*Dv[1]
            x3=x3o+i*lt*Sv[2]+sum(w.take(range(j+1)))*Dv[2]
            index=index+1
            if slip is None:
                savetxt(stdout,[[index,x1,x2,x3,lt,w[j+1],strike,dip,rake]],delimiter=" ",fmt="%4i %8.4f %8.4f %8.4f %8.3f %8.3f %8.2f %5.2f %4.1f")
            else:
                savetxt(stdout,[[index,slip,x1,x2,x3,lt,w[j+1],strike,dip,rake]],delimiter=" ",fmt="%4i %+10.3e %8.4f %8.4f %8.4f %8.3f %8.3f %8.2f %5.2f %4.1f")

    return index
    
def main(): 
	arguments = docopt(__doc__, version='1.0')
	offset=0
	isWithSlip=arguments['--with-slip']
	if isWithSlip: 
		offset+=1
	
	if 1+offset==len(argv):
		fid=stdin
	else:
		fname=argv[1+offset]
		#print fname, len(argv)
		#if not path.isfile(fname):
		#    raise ValueError("invalid file name: " + fname)
		fid=open(fname, 'r')
 
	if isWithSlip:
		print '# nb       slip       x1       x2       x3   length    width   strike   dip  rake'
	else:
		print '# nb       x1       x2       x3   length    width   strike   dip  rake'
        
	k=0
	for line in iter(fid.readlines()):
		if '#'==line[0]:
			continue
		numbers=map(float, line.split())
		s=None
		if 13==len(numbers):
			i,x1,x2,x3,length,width,strike,dip,rake,Lo,Wo,al,aw=numbers
		elif 14==len(numbers):
			if isWithSlip:
				i,s,x1,x2,x3,length,width,strike,dip,rake,Lo,Wo,al,aw=numbers
			else:
				i,x1,x2,x3,length,width,strike,dip,rake,drake,Lo,Wo,al,aw=numbers
		else:
			ValueError("invalid number of column in input file "+fname)

		k=seg2flt(k,x1,x2,x3,length,width,strike,dip,rake,Lo,Wo,al,aw,slip=s)

if __name__ == "__main__":
	main()
