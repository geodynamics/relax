#!/usr/bin/env python

# computes the norm of the residuals between relax prediction and GPS time series

"""
Usage:
    obsres.py [options] <wdir>...
    
Options:
    -b --bounds=<min/max>         time interval to consider.
    -n --network=<opts.dat>       file containing list of stations names.
    -r --range=<mag1/mag2/dmag>   power exponents of the interval of time scales [default: 0/0/1].
    --relax                       use time series of postseismic deformation.
    --vscale                      automatically scales the amplitude of model deformation.
    -w --weight=<wn/we/wd>        relative weight of north, east and down components [default: 1/1/1].
    --ddir=<data/path>            data directory.
 
Description:
    obsres.py computes the weighted sum of the squares of residuals between Relax and GPS data time 
    series. 

    The --range option is for testing overall scaling of the model predictions in the time domain.
    For --range=-1/1/1, the program evaluates the residuals with the modified time t/tm, with tm 
    taking values of 10^(-1), 10^0 and 10^1. If a lower residual is found, the Relax model can be
    changed with all the values of gammadot0 modified to tm*gammadot0.
    
Examples:
    1) obsres.py --ddir=../gps/GPS_Nepal_Tibet --weight=0/0/1 G64_H{20,30,40}_g{0.1,1,10}
    2) obsres.py --bounds=0/7.1 --ddir=./data --range=0/0/1 --network=opts.dat --weight=0/0/1 wdir
 
"""

import numpy as np
import os
import sys
import math
from docopt import docopt

def main():

    # default parameters
	wn=1
	we=1
	wd=1
	ddir='./'
	exp1=0
	exp2=0
	exp3=1
	network='opts.dat'
	index=1
	t1=-np.inf
	t2=+np.inf
	isrelax=False
	isvscale=False

	arguments = docopt(__doc__, version='1.0')
	str1,str2 = arguments['--bounds'].split('/')	 
	if str1 and str2:
		t1=float(str1)
    	t2=float(str2)
	ddir = arguments['--ddir']
	str1,str2,str3 = arguments['--range'].split('/')	 
	if str1 and str2 and str3:
		exp1=float(str1)
		exp2=float(str2)
		exp3=float(str3)
	isrelax=arguments['--relax']
	isvscale=arguments['--vscale']
	str1,str2,str3 = arguments['--weight'].split('/')	
	if str1 and str2 and str3:
		wn=float(str1)
		we=float(str2)
		wd=float(str3)
	network = arguments['--network']
        args = arguments['<wdir>']
	
	exprange=1+int((exp2-exp1)/exp3+0.5)

	print '# obsres.py '+" ".join(sys.argv[1:])
	if 1==exprange:
		print '# '+'model'.ljust(max(map(len,args)))+'  residuals'
	else:
		print '# '+'model'.ljust(max(map(len,args)))+'  residuals time_scale coverage'
	sys.stdout.flush()

	# loop over the models
	for i in xrange(len(args)):
		wdir=args[i]
		fname=wdir+'/'+network
                try:
		    f=file(fname,'r')
                except IOError as e:
                    print >> sys.stderr, '# obsres.py: could not find '+fname+', skipping.'
                    continue

		name=np.loadtxt(f,comments='#',unpack=True,dtype='S4',usecols=[3])
		coverage=np.zeros(len(name))
		#print name

		for e in xrange(exprange):
			n=exp1+float(e)*exp3
			tscale=pow(10,n)

			# initialize norm of residuals
			norm=0
			norm0=0
                        count=0
			index=0
			for s in name:
                            # load data (test upper and lower case)
                            fname=ddir+'/'+s.upper()+'.txt'
			
                            try:
                                with open(fname) as f: pass
                            except IOError as e:
                                fname=ddir+'/'+s.lower()+'.txt'
                            try:
                                with open(fname) as f: pass
                            except IOError as e:
                                fname=ddir+'/'+s.upper()+'.dat'
                            try:
                                with open(fname) as f: pass
                            except IOError as e:
                                fname=ddir+'/'+s.lower()+'.dat'
                            try:
                                with open(fname) as f: pass
                            except IOError as e:
                                # skipping station s
                                print >> sys.stderr, '# obsres.py: could not find '+fname+', skipping.'
                                continue
			
                            f=file(fname,'r')
                            try:
                                tr,nr,er,dr,sn,se,sd=np.loadtxt(f,comments='#',unpack=True,usecols=[0,1,2,3,4,5,6])
                                tr=np.atleast_1d(tr)
                                nr=np.atleast_1d(nr)
                                er=np.atleast_1d(er)
                                dr=np.atleast_1d(dr)
                                sn=np.atleast_1d(sn)
                                se=np.atleast_1d(se)
                                sd=np.atleast_1d(sd)
                            except:
                                print >> sys.stderr, '# obsres.py: error loading file '+s+'. skipping.'
                                coverage[index]=1
                                index+=1
                                continue

                            pos=[]
                            pos=np.logical_and(np.logical_and(np.isfinite(nr+er+dr),tr>=t1),tr<=t2)
                            tr=tr[pos]
                            nr=nr[pos]
                            er=er[pos]
                            dr=dr[pos]
                            sn=sn[pos]
                            se=se[pos]
                            sd=sd[pos]

                            if 0==len(tr):
                                print >> sys.stderr, '# obsres.py: skipping station '+s+' because of insufficient data coverage.'
                                continue

                            # load model
                            fname=wdir+'/'+s+'.txt'
                            f=file(fname,'r')
                            tm,nm,em,dm=np.loadtxt(f,comments='#',unpack=True,usecols=[0,1,2,3])
                            tm=np.atleast_1d(tm)
                            nm=np.atleast_1d(nm)
                            em=np.atleast_1d(em)
                            dm=np.atleast_1d(dm)
                            tm=tm/tscale

                            # time interval
                            tmax=min([max(tm),max(tr)])
                            pos=[]
                            pos=tr<=tmax
                            coverage[index]=float(sum(pos))/float(len(pos))
                            index+=1
                            tr=tr[pos]
                            nr=nr[pos]
                            er=er[pos]
                            dr=dr[pos]
                            sn=sn[pos]
                            se=se[pos]
                            sd=sd[pos]

                            if 0==len(tr):
                                print >> sys.stderr, '# obsres.py: skipping station '+s+' because of insufficient model coverage. Try reducing scaling.'
                                continue

                            if isrelax:
                                nm-=nm[0]
                                em-=em[0]
                                dm-=dm[0]

                            #print [tr,tm]
                            nm=np.interp(tr,tm,nm)
                            em=np.interp(tr,tm,em)
                            dm=np.interp(tr,tm,dm)
		
                            if isvscale:
                                if 0<len(tr):
                                    scale=(np.std(nm)/np.std(nr)*wn+np.std(em)/np.std(er)*we+np.std(dm)/np.std(dr)*wd)/(wn+we+wd)
                                    nm/=scale
                                    em/=scale
                                    dm/=scale

                            dnorm0=sum(pow(nr/sn,2)*wn+pow(er/se,2)*we+pow(dr/sd,2)*wd)
                            dnorm=sum(pow((nr-nm)/sn,2)*wn+pow((er-em)/se,2)*we+pow((dr-dm)/sd,2)*wd)
                            count+=wn*float(len(nm))+we*float(len(em))+wd*float(len(em))
                            norm+=dnorm
                            norm0+=dnorm0
                            #print '{0:8.2e}'.format(dnorm)

                        if 1==exprange:
                            print '{0}    {1:8.2e}'.format(wdir.ljust(max(map(len,args))+1),norm/count)
                        else:
                            print '{0}    {1:8.2e}  {2:9.2e} {3:8.2e}'.format(wdir.ljust(max(map(len,args))+1),norm/count,tscale,sum(coverage)/float(len(coverage)))
                            #print 1. - norm/norm0
                        sys.stdout.flush()

if __name__ == "__main__":
    main()
