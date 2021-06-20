      subroutine ctfft (data,n,ndim,isign,iform,work,nwork)             
c     cooley-tukey fast fourier transform in usasi basic fortran.       
c     multi-dimensional transform, dimensions of arbitrary size,        
c     complex or real data.  n points can be transformed in time        
c     proportional to n*log(n), whereas other methods take n**2 time.   
c     furthermore, less error is built up.  written by norman brenner   
c     of mit lincoln laboratory, june 1968.                             
c                                                                       
c     dimension data(n(1),n(2),...),transform(n(1),n(2),...),n(ndim)    
c     transform(k1,k2,...) = sum(data(j1,j2,...)*exp(isign*2*pi*sqrt(-1)
c     *((j1-1)*(k1-1)/n(1)+(j2-1)*(k2-1)/n(2)+...))), summed for all    
c     j1 and k1 from 1 to n(1), j2 and k2 from 1 to n(2), etc. for all  
c     ndim subscripts.  ndim must be positive and each n(idim) may be   
c     any integer.  isign is +1 or -1.  let ntot = n(1)*n(2)...         
c     ...*n(ndim).  then a -1 transform followed by a +1 one            
c     (or vice versa) returns ntot times the original data.             
c     iform = 1, 0 or -1, as data is complex, real or the               
c     first half of a complex array.  transform values are              
c     returned to array data.  they are complex, real or                
c     the first half of a complex array, as iform = 1, -1 or 0.         
c     the transform of a real array (iform = 0) dimensioned n(1) by n(2)
c     by ... will be returned in the same array, now considered to      
c     be complex of dimensions n(1)/2+1 by n(2) by ....  note that if   
c     iform = 0 or -1, n(1) must be even, and enough room must be       
c     reserved.  the missing values may be obtained by complex conju-   
c     gation.  the reverse transformation, of a half complex array      
c     dimensioned n(1)/2+1 by n(2) by ..., is accomplished setting iform
c     to -1.  in the n array, n(1) must be the true n(1), not n(1)/2+1. 
c     the transform will be real and returned to the input array.       
c     work is a one-dimensional complex array used for working storage. 
c     its length, nwork, need never be larger than the largest n(idim)  
c     and frequently may be much smaller.  fourt computes the minimum   
c     length working storage required and checks that nwork is at least 
c     as long.  this minimum length is ccomputed as shown below.        
c                                                                       
c     for example--                                                     
c     dimension data(1960),work(10)                                     
c     complex data,work                                                 
c     call fourt(data,1960,1,-1,+1,work,10)                             
c                                                                       
c     the multi-dimensional transform is broken down into one-dimen-    
c     sional transforms of length n(idim).  these are further broken    
c     down into transforms of length ifact(if), where these are the     
c     prime factors of n(idim).  for example, n(1) = 1960, ifact(if) =  
c     2, 2, 2, 5, 7 and 7.  the running time is proportional to ntot *  
c     sum(ifact(if)), though factors of two and three will run espe-    
c     cially fast.  naive transform programs will run in time ntot**2.  
c     arrays whose size ntot is prime will run much slower than those   
c     with composite ntot.  for example, ntot = n(1) = 1951 (a prime),  
c     running time will be 1951*1951, while for ntot = 1960, it will    
c     be 1960*(2+2+2+5+7+7), a speedup of eighty times.  naive calcul-  
c     ation will run both in the slower time.  if an array is of        
c     inconvenient length, simply add zeroes to pad it out.  the results
c     will be interpolated according to the new length (see below).     
c                                                                       
c     a fourier transform of length ifact(if) requires a work array     
c     of that length.  therefore, nwork must be as big as the largest   
c     prime factor.  further, work is needed for digit reversal--       
c     each n(idim) (but n(1)/2 if iform = 0 or -1) is factored symmetri-
c     cally, and nwork must be as big as the center factor.  (to factor 
c     symmetrically, separate pairs of identical factors to the flanks, 
c     combining all leftovers in the center.)  for example, n(1) = 1960 
c     =2*2*2*5*7*7=2*7*10*7*2, so nwork must at least max(7,10) = 10.   
c                                                                       
c     an upper bound for the rms relative error is given by gentleman   
c     and sande (3)-- 3 * 2**(-b) * sum(f**1.5), where 2**(-b) is the   
c     smallest bit in the floating point fraction and the sum is over   
c     the prime factors of ntot.                                        
c                                                                       
c     if the input data are a time series, with index j representing    
c     a time (j-1)*deltat, then the corresponding index k in the        
c     transform represents the frequency (k-1)*2*pi/(n*deltat), which   
c     by periodicity, is the same as frequency -(n-k+1)*2*pi/(n*deltat).
c     this is true for n = each n(idim) independently.                  
c                                                                       
c     references--                                                      
c     1.  cooley, j.w. and tukey, j.w., an algorithm for the machine    
c     calculation of complex fourier series.  math. comp., 19, 90,      
c     (april 1967), 297-301.                                            
c     2.  rader, c., et al., what is the fast fourier transform, ieee   
c     transactions on audio and electroacoustics, au-15, 2 (june 1967). 
c     (special issue on the fast fourier transform and its applications)
c     3.  gentleman, w.m. and sande, g., fast fourier transforms--      
c     for fun and profit.  1966 fall joint comp. conf., spartan books,  
c     washington, 1966.                                                 
c     4.  goertzel, g., an algorithm for the evaluation of finite       
c     trigonometric series.  am. math. mo., 65, (1958), 34-35.          
c     5.  singleton, r.c., a method for computing the fast fourier      
c     transform with auxiliary memory and limited high-speed storage.   
c     in (2).                                                           
      dimension data(*), n(1), work(*), ifsym(32), ifcnt(10), ifact(32) 
      if (iform) 10,10,40                                               
 10   if (n(1)-2*(n(1)/2)) 20,40,20                                     
 20   continue
c20   write (6,30) iform,(n(idim),idim=1,ndim)                          
c30   format ('error in fourt.  iform = ',i2,'(real or half-complex)'
c    $' but n(1) is not even./14h dimensions = ',20i5)                  
      return                                                            
 40   ntot=1                                                            
      do 50 idim=1,ndim                                                 
 50   ntot=ntot*n(idim)                                                 
      nrem=ntot                                                         
      if (iform) 60,70,70                                               
 60   nrem=1                                                            
      ntot=(ntot/n(1))*(n(1)/2+1)                                       
c     loop over all dimensions.                                         
 70   do 230 jdim=1,ndim                                                
      if (iform) 80,90,90                                               
 80   idim=ndim+1-jdim                                                  
      go to 100                                                         
 90   idim=jdim                                                         
      nrem=nrem/n(idim)                                                 
 100  ncurr=n(idim)                                                     
      if (idim-1) 110,110,140                                           
 110  if (iform) 120,130,140                                            
 120  call fixrl (data,n(1),nrem,isign,iform)                           
      ntot=(ntot/(n(1)/2+1))*n(1)                                       
 130  ncurr=ncurr/2                                                     
 140  if (ncurr-1) 190,190,150                                          
c     factor n(idim), the length of this dimension.                     
 150  call factr (ncurr,ifact,nfact)                                    
      ifmax=ifact(nfact)                                                
c     arrange the factors symmetrically for simpler digit reversal.     
      call smfac (ifact,nfact,isym,ifsym,nfsym,icent,ifcnt,nfcnt)       
      ifmax=max0(ifmax,icent)                                           
      if (ifmax-nwork) 180,180,160                                      
  160 continue
c 160 write (6,170) nwork,idim,ncurr,icent,(ifact(if),if=1,nfact)       
c 170 format (26h0error in fourt.  nwork = ,i4,20h is too small for n(, 
c    $i1,4h) = ,i5,17h, whose center = ,i4,31h, and whose prime factors 
c    $are--/(1x,20i5))                                                  
      return                                                            
 180  nprev=ntot/(n(idim)*nrem)                                         
c     digit reverse on symmetric factors, for example 2*7*6*7*2.        
      call symrv (data,nprev,ncurr,nrem,ifsym,nfsym)                    
c     digit reverse the asymmetric center, for example, on 6 = 2*3.     
      call asmrv (data,nprev*isym,icent,isym*nrem,ifcnt,nfcnt,work)     
c     fourier transform on each factor, for example, on 2,7,2,3,7 and 2.
      call cool (data,nprev,ncurr,nrem,isign,ifact,work)                
 190  if (iform) 200,210,230                                            
 200  nrem=nrem*n(idim)                                                 
      go to 230                                                         
 210  if (idim-1) 220,220,230                                           
 220  call fixrl (data,n(1),nrem,isign,iform)                           
      ntot=ntot/n(1)*(n(1)/2+1)                                         
 230  continue                                                          
      return                                                            
      end                                                               
      subroutine asmrv (data,nprev,n,nrem,ifact,nfact,work)             
c     shuffle the data array by reversing the digits of one index.      
c     the operation is the same as in symrv, except that the factors    
c     need not be symmetrically arranged, i.e., generally ifact(if) not=
c     ifact(nfact+1-if).  consequently, a work array of length n is     
c     needed.                                                           
      dimension data(*), work(*), ifact(1)                              
      if (nfact-1) 60,60,10                                             
 10   ip0=2                                                             
      ip1=ip0*nprev                                                     
      ip4=ip1*n                                                         
      ip5=ip4*nrem                                                      
      do 50 i1=1,ip1,ip0                                                
      do 50 i5=i1,ip5,ip4                                               
      iwork=1                                                           
      i4rev=i5                                                          
      i4max=i5+ip4-ip1                                                  
      do 40 i4=i5,i4max,ip1                                             
      work(iwork)=data(i4rev)                                           
      work(iwork+1)=data(i4rev+1)                                       
      ip3=ip4                                                           
      do 30 if=1,nfact                                                  
      ip2=ip3/ifact(if)                                                 
      i4rev=i4rev+ip2                                                   
      if (i4rev-ip3-i5) 40,20,20                                        
 20   i4rev=i4rev-ip3                                                   
 30   ip3=ip2                                                           
 40   iwork=iwork+ip0                                                   
      iwork=1                                                           
      do 50 i4=i5,i4max,ip1                                             
      data(i4)=work(iwork)                                              
      data(i4+1)=work(iwork+1)                                          
 50   iwork=iwork+ip0                                                   
 60   return                                                            
      end                                                               
      subroutine cool (data,nprev,n,nrem,isign,ifact,work)              
c     fourier transform of length n.  in place cooley-tukey method,     
c     digit-reversed to normal order, sande-tukey factoring (2).        
c     dimension data(nprev,n,nrem)                                      
c     complex data                                                      
c     data(i1,j2,i3) = sum(data(i1,i2,i3)*exp(isign*2*pi*i*((i2-1)*     
c     (j2-1)/n))), summed over i2 = 1 to n for all i1 from 1 to nprev,  
c     j2 from 1 to n and i3 from 1 to nrem.  the factors of n are given 
c     in any order in array ifact.  factors of two are done in pairs    
c     as much as possible (fourier transform of length four), factors of
c     three are done separately, and all factors five or higher         
c     are done by goertzel's algorithm (4).                             
      dimension data(*), work(*), ifact(1)                              
      twopi=6.283185307*float(isign)                                    
      ip0=2                                                             
      ip1=ip0*nprev                                                     
      ip4=ip1*n                                                         
      ip5=ip4*nrem                                                      
      if=0                                                              
      ip2=ip1                                                           
 10   if (ip2-ip4) 20,240,240                                           
 20   if=if+1                                                           
      ifcur=ifact(if)                                                   
      if (ifcur-2) 60,30,60                                             
 30   if (4*ip2-ip4) 40,40,60                                           
 40   if (ifact(if+1)-2) 60,50,60                                       
 50   if=if+1                                                           
      ifcur=4                                                           
 60   ip3=ip2*ifcur                                                     
      theta=twopi/float(ifcur)                                          
      sinth=sin(theta/2.)                                               
      rootr=-2.*sinth*sinth                                             
c     cos(theta)-1, for accuracy.                                       
      rooti=sin(theta)                                                  
      theta=twopi/float(ip3/ip1)                                        
      sinth=sin(theta/2.)                                               
      wstpr=-2.*sinth*sinth                                             
      wstpi=sin(theta)                                                  
      wr=1.                                                             
      wi=0.                                                             
      do 230 i2=1,ip2,ip1                                               
      if (ifcur-4) 70,70,210                                            
 70   if ((i2-1)*(ifcur-2)) 240,90,80                                   
 80   w2r=wr*wr-wi*wi                                                   
      w2i=2.*wr*wi                                                      
      w3r=w2r*wr-w2i*wi                                                 
      w3i=w2r*wi+w2i*wr                                                 
 90   i1max=i2+ip1-ip0                                                  
      do 200 i1=i2,i1max,ip0                                            
      do 200 i5=i1,ip5,ip3                                              
      j0=i5                                                             
      j1=j0+ip2                                                         
      j2=j1+ip2                                                         
      j3=j2+ip2                                                         
      if (i2-1) 140,140,100                                             
 100  if (ifcur-3) 130,120,110                                          
c     apply the phase shift factors                                     
 110  tempr=data(j3)                                                    
      data(j3)=w3r*tempr-w3i*data(j3+1)                                 
      data(j3+1)=w3r*data(j3+1)+w3i*tempr                               
      tempr=data(j2)                                                    
      data(j2)=wr*tempr-wi*data(j2+1)                                   
      data(j2+1)=wr*data(j2+1)+wi*tempr                                 
      tempr=data(j1)                                                    
      data(j1)=w2r*tempr-w2i*data(j1+1)                                 
      data(j1+1)=w2r*data(j1+1)+w2i*tempr                               
      go to 140                                                         
 120  tempr=data(j2)                                                    
      data(j2)=w2r*tempr-w2i*data(j2+1)                                 
      data(j2+1)=w2r*data(j2+1)+w2i*tempr                               
 130  tempr=data(j1)                                                    
      data(j1)=wr*tempr-wi*data(j1+1)                                   
      data(j1+1)=wr*data(j1+1)+wi*tempr                                 
 140  if (ifcur-3) 150,160,170                                          
c     do a fourier transform of length two                              
 150  tempr=data(j1)                                                    
      tempi=data(j1+1)                                                  
      data(j1)=data(j0)-tempr                                           
      data(j1+1)=data(j0+1)-tempi                                       
      data(j0)=data(j0)+tempr                                           
      data(j0+1)=data(j0+1)+tempi                                       
      go to 200                                                         
c     do a fourier transform of length three                            
 160  sumr=data(j1)+data(j2)                                            
      sumi=data(j1+1)+data(j2+1)                                        
      tempr=data(j0)-.5*sumr                                            
      tempi=data(j0+1)-.5*sumi                                          
      data(j0)=data(j0)+sumr                                            
      data(j0+1)=data(j0+1)+sumi                                        
      difr=rooti*(data(j2+1)-data(j1+1))                                
      difi=rooti*(data(j1)-data(j2))                                    
      data(j1)=tempr+difr                                               
      data(j1+1)=tempi+difi                                             
      data(j2)=tempr-difr                                               
      data(j2+1)=tempi-difi                                             
      go to 200                                                         
c     do a fourier transform of length four (from bit reversed order)   
 170  t0r=data(j0)+data(j1)                                             
      t0i=data(j0+1)+data(j1+1)                                         
      t1r=data(j0)-data(j1)                                             
      t1i=data(j0+1)-data(j1+1)                                         
      t2r=data(j2)+data(j3)                                             
      t2i=data(j2+1)+data(j3+1)                                         
      t3r=data(j2)-data(j3)                                             
      t3i=data(j2+1)-data(j3+1)                                         
      data(j0)=t0r+t2r                                                  
      data(j0+1)=t0i+t2i                                                
      data(j2)=t0r-t2r                                                  
      data(j2+1)=t0i-t2i                                                
      if (isign) 180,180,190                                            
 180  t3r=-t3r                                                          
      t3i=-t3i                                                          
 190  data(j1)=t1r-t3i                                                  
      data(j1+1)=t1i+t3r                                                
      data(j3)=t1r+t3i                                                  
      data(j3+1)=t1i-t3r                                                
 200  continue                                                          
      go to 220                                                         
c     do a fourier transform of length five or more                     
 210  call goert (data(i2),nprev,ip2/ip1,ifcur,ip5/ip3,work,wr,wi,rootr,
     $rooti)                                                            
 220  tempr=wr                                                          
      wr=wstpr*tempr-wstpi*wi+tempr                                     
 230  wi=wstpr*wi+wstpi*tempr+wi                                        
      ip2=ip3                                                           
      go to 10                                                          
 240  return                                                            
      end                                                               
      subroutine factr (n,ifact,nfact)                                  
c     factor n into its prime factors, nfact in number.  for example,   
c     for n = 1960, nfact = 6 and ifact(if) = 2, 2, 2, 5, 7 and 7.      
      dimension ifact(1)                                                
      if=0                                                              
      npart=n                                                           
      do 50 id=1,n,2                                                    
      idiv=id                                                           
      if (id-1) 10,10,20                                                
 10   idiv=2                                                            
 20   iquot=npart/idiv                                                  
      if (npart-idiv*iquot) 40,30,40                                    
 30   if=if+1                                                           
      ifact(if)=idiv                                                    
      npart=iquot                                                       
      go to 20                                                          
 40   if (iquot-idiv) 60,60,50                                          
 50   continue                                                          
 60   if (npart-1) 80,80,70                                             
 70   if=if+1                                                           
      ifact(if)=npart                                                   
 80   nfact=if                                                          
      return                                                            
      end                                                               
      subroutine fixrl (data,n,nrem,isign,iform)                        
c     for iform = 0, convert the transform of a doubled-up real array,  
c     considered complex, into its true transform.  supply only the     
c     first half of the complex transform, as the second half has       
c     conjugate symmetry.  for iform = -1, convert the first half       
c     of the true transform into the transform of a doubled-up real     
c     array.  n must be even.                                           
c     using complex notation and subscripts starting at zero, the       
c     transformation is--                                               
c     dimension data(n,nrem)                                            
c     zstp = exp(isign*2*pi*i/n)                                        
c     do 10 i2=0,nrem-1                                                 
c     data(0,i2) = conj(data(0,i2))*(1+i)                               
c     do 10 i1=1,n/4                                                    
c     z = (1+(2*iform+1)*i*zstp**i1)/2                                  
c     i1cnj = n/2-i1                                                    
c     dif = data(i1,i2)-conj(data(i1cnj,i2))                            
c     temp = z*dif                                                      
c     data(i1,i2) = (data(i1,i2)-temp)*(1-iform)                        
c 10  data(i1cnj,i2) = (data(i1cnj,i2)+conj(temp))*(1-iform)            
c     if i1=i1cnj, the calculation for that value collapses into        
c     a simple conjugation of data(i1,i2).                              
      dimension data(*)                                                 
      twopi=6.283185307*float(isign)                                    
      ip0=2                                                             
      ip1=ip0*(n/2)                                                     
      ip2=ip1*nrem                                                      
      if (iform) 10,70,70                                               
c     pack the real input values (two per column)                       
 10   j1=ip1+1                                                          
      data(2)=data(j1)                                                  
      if (nrem-1) 70,70,20                                              
 20   j1=j1+ip0                                                         
      i2min=ip1+1                                                       
      do 60 i2=i2min,ip2,ip1                                            
      data(i2)=data(j1)                                                 
      j1=j1+ip0                                                         
      if (n-2) 50,50,30                                                 
 30   i1min=i2+ip0                                                      
      i1max=i2+ip1-ip0                                                  
      do 40 i1=i1min,i1max,ip0                                          
      data(i1)=data(j1)                                                 
      data(i1+1)=data(j1+1)                                             
 40   j1=j1+ip0                                                         
 50   data(i2+1)=data(j1)                                               
 60   j1=j1+ip0                                                         
 70   do 80 i2=1,ip2,ip1                                                
      tempr=data(i2)                                                    
      data(i2)=data(i2)+data(i2+1)                                      
 80   data(i2+1)=tempr-data(i2+1)                                       
      if (n-2) 200,200,90                                               
 90   theta=twopi/float(n)                                              
      sinth=sin(theta/2.)                                               
      zstpr=-2.*sinth*sinth                                             
      zstpi=sin(theta)                                                  
      zr=(1.-zstpi)/2.                                                  
      zi=(1.+zstpr)/2.                                                  
      if (iform) 100,110,110                                            
 100  zr=1.-zr                                                          
      zi=-zi                                                            
 110  i1min=ip0+1                                                       
      i1max=ip0*(n/4)+1                                                 
      do 190 i1=i1min,i1max,ip0                                         
      do 180 i2=i1,ip2,ip1                                              
      i2cnj=ip0*(n/2+1)-2*i1+i2                                         
      if (i2-i2cnj) 150,120,120                                         
 120  if (isign*(2*iform+1)) 130,140,140                                
 130  data(i2+1)=-data(i2+1)                                            
 140  if (iform) 170,180,180                                            
 150  difr=data(i2)-data(i2cnj)                                         
      difi=data(i2+1)+data(i2cnj+1)                                     
      tempr=difr*zr-difi*zi                                             
      tempi=difr*zi+difi*zr                                             
      data(i2)=data(i2)-tempr                                           
      data(i2+1)=data(i2+1)-tempi                                       
      data(i2cnj)=data(i2cnj)+tempr                                     
      data(i2cnj+1)=data(i2cnj+1)-tempi                                 
      if (iform) 160,180,180                                            
 160  data(i2cnj)=data(i2cnj)+data(i2cnj)                               
      data(i2cnj+1)=data(i2cnj+1)+data(i2cnj+1)                         
 170  data(i2)=data(i2)+data(i2)                                        
      data(i2+1)=data(i2+1)+data(i2+1)                                  
 180  continue                                                          
      tempr=zr-.5                                                       
      zr=zstpr*tempr-zstpi*zi+zr                                        
 190  zi=zstpr*zi+zstpi*tempr+zi                                        
c     recursion saves time, at a slight loss in accuracy.  if available,
c     use double precision to compute zr and zi.                        
 200  if (iform) 270,210,210                                            
c     unpack the real transform values (two per column)                 
 210  i2=ip2+1                                                          
      i1=i2                                                             
      j1=ip0*(n/2+1)*nrem+1                                             
      go to 250                                                         
 220  data(j1)=data(i1)                                                 
      data(j1+1)=data(i1+1)                                             
      i1=i1-ip0                                                         
      j1=j1-ip0                                                         
 230  if (i2-i1) 220,240,240                                            
 240  data(j1)=data(i1)                                                 
      data(j1+1)=0.                                                     
 250  i2=i2-ip1                                                         
      j1=j1-ip0                                                         
      data(j1)=data(i2+1)                                               
      data(j1+1)=0.                                                     
      i1=i1-ip0                                                         
      j1=j1-ip0                                                         
      if (i2-1) 260,260,230                                             
 260  data(2)=0.                                                        
 270  return                                                            
      end                                                               
      subroutine goert(data,nprev,iprod,ifact,irem,work,wminr,wmini,    
     $ rootr,rooti)                                                     
c     phase-shifted fourier transform of length ifact by the goertzel   
c     algorithm (4).  ifact must be odd and at least 5.  further speed  
c     is gained by computing two transform values at the same time.     
c     dimension data(nprev,iprod,ifact,irem)                            
c     data(i1,1,j3,i5) = sum(data(i1,1,i3,i5) * w**(i3-1)), summed      
c     over i3 = 1 to ifact for all i1 from 1 to nprev, j3 from 1 to     
c     ifact and i5 from 1 to irem.                                      
c     w = wmin * exp(isign*2*pi*i*(j3-1)/ifact).                        
      dimension data(*), work(*)                                        
      ip0=2                                                             
      ip1=ip0*nprev                                                     
      ip2=ip1*iprod                                                     
      ip3=ip2*ifact                                                     
      ip5=ip3*irem                                                      
      if (wmini) 10,40,10                                               
c     apply the phase shift factors                                     
 10   wr=wminr                                                          
      wi=wmini                                                          
      i3min=1+ip2                                                       
      do 30 i3=i3min,ip3,ip2                                            
      i1max=i3+ip1-ip0                                                  
      do 20 i1=i3,i1max,ip0                                             
      do 20 i5=i1,ip5,ip3                                               
      tempr=data(i5)                                                    
      data(i5)=wr*tempr-wi*data(i5+1)                                   
 20   data(i5+1)=wr*data(i5+1)+wi*tempr                                 
      tempr=wr                                                          
      wr=wminr*tempr-wmini*wi                                           
 30   wi=wminr*wi+wmini*tempr                                           
 40   do 90 i1=1,ip1,ip0                                                
      do 90 i5=i1,ip5,ip3                                               
c     straight summation for the first term                             
      sumr=0.                                                           
      sumi=0.                                                           
      i3max=i5+ip3-ip2                                                  
      do 50 i3=i5,i3max,ip2                                             
      sumr=sumr+data(i3)                                                
 50   sumi=sumi+data(i3+1)                                              
      work(1)=sumr                                                      
      work(2)=sumi                                                      
      wr=rootr+1.                                                       
      wi=rooti                                                          
      iwmin=1+ip0                                                       
      iwmax=ip0*((ifact+1)/2)-1                                         
      do 80 iwork=iwmin,iwmax,ip0                                       
      twowr=wr+wr                                                       
      i3=i3max                                                          
      oldsr=0.                                                          
      oldsi=0.                                                          
      sumr=data(i3)                                                     
      sumi=data(i3+1)                                                   
      i3=i3-ip2                                                         
 60   tempr=sumr                                                        
      tempi=sumi                                                        
      sumr=twowr*sumr-oldsr+data(i3)                                    
      sumi=twowr*sumi-oldsi+data(i3+1)                                  
      oldsr=tempr                                                       
      oldsi=tempi                                                       
      i3=i3-ip2                                                         
      if (i3-i5) 70,70,60                                               
c     in a fourier transform the w corresponding to the point at k      
c     is the conjugate of that at ifact-k (that is, exp(twopi*i*        
c     k/ifact) = conj(exp(twopi*i*(ifact-k)/ifact))).  since the        
c     main loop of goertzels algorithm is indifferent to the imaginary  
c     part of w, it need be supplied only at the end.                   
 70   tempr=-wi*sumi                                                    
      tempi=wi*sumr                                                     
      sumr=wr*sumr-oldsr+data(i3)                                       
      sumi=wr*sumi-oldsi+data(i3+1)                                     
      work(iwork)=sumr+tempr                                            
      work(iwork+1)=sumi+tempi                                          
      iwcnj=ip0*(ifact+1)-iwork                                         
      work(iwcnj)=sumr-tempr                                            
      work(iwcnj+1)=sumi-tempi                                          
c     singleton's recursion, for accuracy and speed (5).                
      tempr=wr                                                          
      wr=wr*rootr-wi*rooti+wr                                           
 80   wi=tempr*rooti+wi*rootr+wi                                        
      iwork=1                                                           
      do 90 i3=i5,i3max,ip2                                             
      data(i3)=work(iwork)                                              
      data(i3+1)=work(iwork+1)                                          
 90   iwork=iwork+ip0                                                   
      return                                                            
      end                                                               
      subroutine smfac (ifact,nfact,isym,ifsym,nfsym,icent,ifcnt,nfcnt) 
c     rearrange the prime factors of n into a square and a non-         
c     square.  n = isym*icent*isym, where icent is square-free.         
c     isym = ifsym(1)*...*ifsym(nfsym), each a prime factor.            
c     icent = ifcnt(1)*...*ifcnt(nfcnt), each a prime factor.           
c     for example, n = 1960 = 14*10*14.  then isym = 14, icent = 10,    
c     nfsym = 2, nfcnt = 2, nfact = 6, ifsym(ifs) = 2, 7, ifcnt(ifc) =  
c     2, 5 and ifact(if) = 2, 7, 2, 5, 7, 2.                            
      dimension ifsym(1), ifcnt(1), ifact(1)                            
      isym=1                                                            
      icent=1                                                           
      ifs=0                                                             
      ifc=0                                                             
      if=1                                                              
 10   if (if-nfact) 20,40,50                                            
 20   if (ifact(if)-ifact(if+1)) 40,30,40                               
 30   ifs=ifs+1                                                         
      ifsym(ifs)=ifact(if)                                              
      isym=ifact(if)*isym                                               
      if=if+2                                                           
      go to 10                                                          
 40   ifc=ifc+1                                                         
      ifcnt(ifc)=ifact(if)                                              
      icent=ifact(if)*icent                                             
      if=if+1                                                           
      go to 10                                                          
 50   nfsym=ifs                                                         
      nfcnt=ifc                                                         
      nfsm2=2*nfsym                                                     
      nfact=2*nfsym+nfcnt                                               
      if (nfcnt) 80,80,60                                               
 60   nfsm2=nfsm2+1                                                     
      ifsym(nfsym+1)=icent                                              
      do 70 ifc=1,nfcnt                                                 
      if=nfsym+ifc                                                      
 70   ifact(if)=ifcnt(ifc)                                              
 80   if (nfsym) 110,110,90                                             
 90   do 100 ifs=1,nfsym                                                
      ifscj=nfsm2+1-ifs                                                 
      ifsym(ifscj)=ifsym(ifs)                                           
      ifact(ifs)=ifsym(ifs)                                             
      ifcnj=nfact+1-ifs                                                 
 100  ifact(ifcnj)=ifsym(ifs)                                           
 110  nfsym=nfsm2                                                       
      return                                                            
      end                                                               
      subroutine symrv (data,nprev,n,nrem,ifact,nfact)                  
c     shuffle the data array by reversing the digits of one index.      
c     dimension data(nprev,n,nrem)                                      
c     replace data(i1,i2,i3) by data(i1,i2rev,i3) for all i1 from 1 to  
c     nprev, i2 from 1 to n and i3 from 1 to nrem.  i2rev-1 is the      
c     integer whose digit representation in the multi-radix notation    
c     of factors ifact(if) is the reverse of the representation of i2-1.
c     for example, if all ifact(if) = 2, i2-1 = 11001, i2rev-1 = 10011. 
c     the factors must be symmetrically arranged, i.e., ifact(if) =     
c     ifact(nfact+1-if).                                                
      dimension data(*), ifact(1)                                       
      if (nfact-1) 80,80,10                                             
 10   ip0=2                                                             
      ip1=ip0*nprev                                                     
      ip4=ip1*n                                                         
      ip5=ip4*nrem                                                      
      i4rev=1                                                           
      do 70 i4=1,ip4,ip1                                                
      if (i4-i4rev) 20,40,40                                            
 20   i1max=i4+ip1-ip0                                                  
      do 30 i1=i4,i1max,ip0                                             
      do 30 i5=i1,ip5,ip4                                               
      i5rev=i4rev+i5-i4                                                 
      tempr=data(i5)
      tempi=data(i5+1)                                                  
      data(i5)=data(i5rev)                                              
      data(i5+1)=data(i5rev+1)                                          
      data(i5rev)=tempr                                                 
 30   data(i5rev+1)=tempi                                               
 40   ip3=ip4                                                           
      do 60 if=1,nfact                                                  
      ip2=ip3/ifact(if)                                                 
      i4rev=i4rev+ip2                                                   
      if (i4rev-ip3) 70,70,50                                           
 50   i4rev=i4rev-ip3                                                   
 60   ip3=ip2                                                           
 70   continue                                                          
 80   return                                                            
      end                                                               
