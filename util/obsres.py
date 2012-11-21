#!/sw/bin/python

# computes the norm of the residuals between relax prediction and GPS time series

import getopt
import numpy as np
import os
import sys
import math

def usage():
    print 'obsres.py computes residuals between Relax and GPS data time series'
    print ''
    print 'usage: obsres.py --bounds=0/7.1 --ddir=./data --range=0/0/1 \\'
    print '                 --network=opts.dat --weight=0/0/1 wdir'
    print ''
    print 'options:'
    print '  -b --bounds:  time interval to consider'
    print '  -n --network: file containing list of stations names'
    print '  -r --range:   power exponents of the interval of time scales'
    print '  --relax:      use time series of postseismic deformation'
    print '  --vscale:     automatically scales the amplitude of model deformation'
    print '  -w --weight:  relative weight of north, east and down components'
    print ''
    print 'example:'
    print './obsres.py --ddir=../gps/GPS_Nepal_Tibet --weight=0/0/1 G64_H{20,30,40}_g{0.1,1,10}'
    print ''

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
 
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:d:n:r:w:", ["help","bounds=","ddir=","network=","relax","range=","weight="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print >> sys.stderr, 'obsres.py:', str(err) # will print something like "option -a not recognized"
        print ''
        usage()
        sys.exit(2)

    if 0==len(args):
	usage()
	sys.exit(2)

    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-b", "--bounds"):
            str1, str2 = a.split('/')
            t1=float(str1)
            t2=float(str2)
        elif o in ("-d", "--ddir"):
            ddir = a
        elif o in ("-r", "--range"):
            str1, str2, str3 = a.split('/')
            exp1=float(str1)
            exp2=float(str2)
	    exp3=float(str3)
        elif o == "--relax":
	    isrelax=True
        elif o == "--vscale":
	    isvscale=True
        elif o in ("-w", "--weight"):
            str1, str2, str3 = a.split('/')
            wn=float(str1)
            we=float(str2)
            wd=float(str3)
        elif o in ("-n", "--network"):
            network = a
        else:
            print >> sys.stderr, 'obsres.py: unhandled option:', o, a
            assert False, "unhandled option"

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
        f=file(fname,'r')
        name=np.loadtxt(f,comments='#',unpack=True,dtype='S4',usecols=[3])
	coverage=np.zeros(len(name))
        #print name

        for e in xrange(exprange):
            n=exp1+float(e)*exp3
            tscale=pow(10,n)

	    # initialize norm of residuals
	    norm=0
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
                        # skipping station s
                        print >> sys.stderr, 'obsres.py: could not fine '+fname+', skipping.'
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
                    print >> sys.stderr, 'obsres.py: error loading file '+s+'. skipping.'
		    coverage[index]=1
		    index+=1
		    continue
   		pos=np.logical_and(np.logical_and(np.isfinite(nr+er+dr),tr>=t1),tr<=t2)
		tr=tr[pos]
		nr=nr[pos]
		er=er[pos]
		dr=dr[pos]
		sn=sn[pos]
		se=se[pos]
		sd=sd[pos]
		pos=[]

		if 0==len(tr):
                    print >> sys.stderr, 'obsres.py: skipping station '+s+' because of insufficient data coverage.'
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
		pos=[]

		if 0==len(tr):
                    print >> sys.stderr, 'obsres.py: skipping station '+s+' because of insufficient model coverage. Try reducing scaling.'
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

		dnorm=sum(pow(nr-nm,2)/sn*wn+pow(er-em,2)/se*we+pow(dr-dm,2)/sd*wd)
		norm+=dnorm
		#print '{0:8.2e}'.format(dnorm)

    	    if 1==exprange:
	        print '{0}    {1:8.2e}'.format(wdir.ljust(max(map(len,args))+1),norm)
            else:
	        print '{0}    {1:8.2e}  {2:9.2e} {3:8.2e}'.format(wdir.ljust(max(map(len,args))+1),norm,tscale,sum(coverage)/float(len(coverage)))
	    sys.stdout.flush()

if __name__ == "__main__":
    main()

