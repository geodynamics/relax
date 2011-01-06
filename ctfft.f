      subroutine ctfft (data,n,ndim,isign,iform,work,nwork)             fft   1
c     cooley-tukey fast fourier transform in usasi basic fortran.       fft   2
c     multi-dimensional transform, dimensions of arbitrary size,        fft   3
c     complex or real data.  n points can be transformed in time        fft   4
c     proportional to n*log(n), whereas other methods take n**2 time.   fft   5
c     furthermore, less error is built up.  written by norman brenner   fft   6
c     of mit lincoln laboratory, june 1968.                             fft   7
c                                                                       fft   8
c     dimension data(n(1),n(2),...),transform(n(1),n(2),...),n(ndim)    fft   9
c     transform(k1,k2,...) = sum(data(j1,j2,...)*exp(isign*2*pi*sqrt(-1)fft  10
c     *((j1-1)*(k1-1)/n(1)+(j2-1)*(k2-1)/n(2)+...))), summed for all    fft  11
c     j1 and k1 from 1 to n(1), j2 and k2 from 1 to n(2), etc. for all  fft  12
c     ndim subscripts.  ndim must be positive and each n(idim) may be   fft  13
c     any integer.  isign is +1 or -1.  let ntot = n(1)*n(2)...         fft  14
c     ...*n(ndim).  then a -1 transform followed by a +1 one            fft  15
c     (or vice versa) returns ntot times the original data.             fft  16
c     iform = 1, 0 or -1, as data is complex, real or the               fft  17
c     first half of a complex array.  transform values are              fft  18
c     returned to array data.  they are complex, real or                fft  19
c     the first half of a complex array, as iform = 1, -1 or 0.         fft  20
c     the transform of a real array (iform = 0) dimensioned n(1) by n(2)fft  21
c     by ... will be returned in the same array, now considered to      fft  22
c     be complex of dimensions n(1)/2+1 by n(2) by ....  note that if   fft  23
c     iform = 0 or -1, n(1) must be even, and enough room must be       fft  24
c     reserved.  the missing values may be obtained by complex conju-   fft  25
c     gation.  the reverse transformation, of a half complex array      fft  26
c     dimensioned n(1)/2+1 by n(2) by ..., is accomplished setting iformfft  27
c     to -1.  in the n array, n(1) must be the true n(1), not n(1)/2+1. fft  28
c     the transform will be real and returned to the input array.       fft  29
c     work is a one-dimensional complex array used for working storage. fft  30
c     its length, nwork, need never be larger than the largest n(idim)  fft  31
c     and frequently may be much smaller.  fourt computes the minimum   fft  32
c     length working storage required and checks that nwork is at least fft  33
c     as long.  this minimum length is ccomputed as shown below.        fft  34
c                                                                       fft  35
c     for example--                                                     fft  36
c     dimension data(1960),work(10)                                     fft  37
c     complex data,work                                                 fft  38
c     call fourt(data,1960,1,-1,+1,work,10)                             fft  39
c                                                                       fft  40
c     the multi-dimensional transform is broken down into one-dimen-    fft  41
c     sional transforms of length n(idim).  these are further broken    fft  42
c     down into transforms of length ifact(if), where these are the     fft  43
c     prime factors of n(idim).  for example, n(1) = 1960, ifact(if) =  fft  44
c     2, 2, 2, 5, 7 and 7.  the running time is proportional to ntot *  fft  45
c     sum(ifact(if)), though factors of two and three will run espe-    fft  46
c     cially fast.  naive transform programs will run in time ntot**2.  fft  47
c     arrays whose size ntot is prime will run much slower than those   fft  48
c     with composite ntot.  for example, ntot = n(1) = 1951 (a prime),  fft  49
c     running time will be 1951*1951, while for ntot = 1960, it will    fft  50
c     be 1960*(2+2+2+5+7+7), a speedup of eighty times.  naive calcul-  fft  51
c     ation will run both in the slower time.  if an array is of        fft  52
c     inconvenient length, simply add zeroes to pad it out.  the resultsfft  53
c     will be interpolated according to the new length (see below).     fft  54
c                                                                       fft  55
c     a fourier transform of length ifact(if) requires a work array     fft  56
c     of that length.  therefore, nwork must be as big as the largest   fft  57
c     prime factor.  further, work is needed for digit reversal--       fft  58
c     each n(idim) (but n(1)/2 if iform = 0 or -1) is factored symmetri-fft  59
c     cally, and nwork must be as big as the center factor.  (to factor fft  60
c     symmetrically, separate pairs of identical factors to the flanks, fft  61
c     combining all leftovers in the center.)  for example, n(1) = 1960 fft  62
c     =2*2*2*5*7*7=2*7*10*7*2, so nwork must at least max(7,10) = 10.   fft  63
c                                                                       fft  64
c     an upper bound for the rms relative error is given by gentleman   fft  65
c     and sande (3)-- 3 * 2**(-b) * sum(f**1.5), where 2**(-b) is the   fft  66
c     smallest bit in the floating point fraction and the sum is over   fft  67
c     the prime factors of ntot.                                        fft  68
c                                                                       fft  69
c     if the input data are a time series, with index j representing    fft  70
c     a time (j-1)*deltat, then the corresponding index k in the        fft  71
c     transform represents the frequency (k-1)*2*pi/(n*deltat), which   fft  72
c     by periodicity, is the same as frequency -(n-k+1)*2*pi/(n*deltat).fft  73
c     this is true for n = each n(idim) independently.                  fft  74
c                                                                       fft  75
c     references--                                                      fft  76
c     1.  cooley, j.w. and tukey, j.w., an algorithm for the machine    fft  77
c     calculation of complex fourier series.  math. comp., 19, 90,      fft  78
c     (april 1967), 297-301.                                            fft  79
c     2.  rader, c., et al., what is the fast fourier transform, ieee   fft  80
c     transactions on audio and electroacoustics, au-15, 2 (june 1967). fft  81
c     (special issue on the fast fourier transform and its applications)fft  82
c     3.  gentleman, w.m. and sande, g., fast fourier transforms--      fft  83
c     for fun and profit.  1966 fall joint comp. conf., spartan books,  fft  84
c     washington, 1966.                                                 fft  85
c     4.  goertzel, g., an algorithm for the evaluation of finite       fft  86
c     trigonometric series.  am. math. mo., 65, (1958), 34-35.          fft  87
c     5.  singleton, r.c., a method for computing the fast fourier      fft  88
c     transform with auxiliary memory and limited high-speed storage.   fft  89
c     in (2).                                                           fft  90
      dimension data(*), n(1), work(*), ifsym(32), ifcnt(10), ifact(32) fft  91
      if (iform) 10,10,40                                               fft  92
 10   if (n(1)-2*(n(1)/2)) 20,40,20                                     fft  93
 20   continue
c20   write (6,30) iform,(n(idim),idim=1,ndim)                          fft  94
c30   format ('error in fourt.  iform = ',i2,'(real or half-complex)'
c    $' but n(1) is not even./14h dimensions = ',20i5)                  fft  96
      return                                                            fft  97
 40   ntot=1                                                            fft  98
      do 50 idim=1,ndim                                                 fft  99
 50   ntot=ntot*n(idim)                                                 fft 100
      nrem=ntot                                                         fft 101
      if (iform) 60,70,70                                               fft 102
 60   nrem=1                                                            fft 103
      ntot=(ntot/n(1))*(n(1)/2+1)                                       fft 104
c     loop over all dimensions.                                         fft 105
 70   do 230 jdim=1,ndim                                                fft 106
      if (iform) 80,90,90                                               fft 107
 80   idim=ndim+1-jdim                                                  fft 108
      go to 100                                                         fft 109
 90   idim=jdim                                                         fft 110
      nrem=nrem/n(idim)                                                 fft 111
 100  ncurr=n(idim)                                                     fft 112
      if (idim-1) 110,110,140                                           fft 113
 110  if (iform) 120,130,140                                            fft 114
 120  call fixrl (data,n(1),nrem,isign,iform)                           fft 115
      ntot=(ntot/(n(1)/2+1))*n(1)                                       fft 116
 130  ncurr=ncurr/2                                                     fft 117
 140  if (ncurr-1) 190,190,150                                          fft 118
c     factor n(idim), the length of this dimension.                     fft 119
 150  call factr (ncurr,ifact,nfact)                                    fft 120
      ifmax=ifact(nfact)                                                fft 121
c     arrange the factors symmetrically for simpler digit reversal.     fft 122
      call smfac (ifact,nfact,isym,ifsym,nfsym,icent,ifcnt,nfcnt)       fft 123
      ifmax=max0(ifmax,icent)                                           fft 124
      if (ifmax-nwork) 180,180,160                                      fft 125
  160 continue
c 160 write (6,170) nwork,idim,ncurr,icent,(ifact(if),if=1,nfact)       fft 126
c 170 format (26h0error in fourt.  nwork = ,i4,20h is too small for n(, fft 127
c    $i1,4h) = ,i5,17h, whose center = ,i4,31h, and whose prime factors fft 128
c    $are--/(1x,20i5))                                                  fft 129
      return                                                            fft 130
 180  nprev=ntot/(n(idim)*nrem)                                         fft 131
c     digit reverse on symmetric factors, for example 2*7*6*7*2.        fft 132
      call symrv (data,nprev,ncurr,nrem,ifsym,nfsym)                    fft 133
c     digit reverse the asymmetric center, for example, on 6 = 2*3.     fft 134
      call asmrv (data,nprev*isym,icent,isym*nrem,ifcnt,nfcnt,work)     fft 135
c     fourier transform on each factor, for example, on 2,7,2,3,7 and 2.fft 136
      call cool (data,nprev,ncurr,nrem,isign,ifact,work)                fft 137
 190  if (iform) 200,210,230                                            fft 138
 200  nrem=nrem*n(idim)                                                 fft 139
      go to 230                                                         fft 140
 210  if (idim-1) 220,220,230                                           fft 141
 220  call fixrl (data,n(1),nrem,isign,iform)                           fft 142
      ntot=ntot/n(1)*(n(1)/2+1)                                         fft 143
 230  continue                                                          fft 144
      return                                                            fft 145
      end                                                               fft 146-
      subroutine asmrv (data,nprev,n,nrem,ifact,nfact,work)             asm   1
c     shuffle the data array by reversing the digits of one index.      asm   2
c     the operation is the same as in symrv, except that the factors    asm   3
c     need not be symmetrically arranged, i.e., generally ifact(if) not=asm   4
c     ifact(nfact+1-if).  consequently, a work array of length n is     asm   5
c     needed.                                                           asm   6
      dimension data(*), work(*), ifact(1)                              asm   7
      if (nfact-1) 60,60,10                                             asm   8
 10   ip0=2                                                             asm   9
      ip1=ip0*nprev                                                     asm  10
      ip4=ip1*n                                                         asm  11
      ip5=ip4*nrem                                                      asm  12
      do 50 i1=1,ip1,ip0                                                asm  13
      do 50 i5=i1,ip5,ip4                                               asm  14
      iwork=1                                                           asm  15
      i4rev=i5                                                          asm  16
      i4max=i5+ip4-ip1                                                  asm  17
      do 40 i4=i5,i4max,ip1                                             asm  18
      work(iwork)=data(i4rev)                                           asm  19
      work(iwork+1)=data(i4rev+1)                                       asm  20
      ip3=ip4                                                           asm  21
      do 30 if=1,nfact                                                  asm  22
      ip2=ip3/ifact(if)                                                 asm  23
      i4rev=i4rev+ip2                                                   asm  24
      if (i4rev-ip3-i5) 40,20,20                                        asm  25
 20   i4rev=i4rev-ip3                                                   asm  26
 30   ip3=ip2                                                           asm  27
 40   iwork=iwork+ip0                                                   asm  28
      iwork=1                                                           asm  29
      do 50 i4=i5,i4max,ip1                                             asm  30
      data(i4)=work(iwork)                                              asm  31
      data(i4+1)=work(iwork+1)                                          asm  32
 50   iwork=iwork+ip0                                                   asm  33
 60   return                                                            asm  34
      end                                                               asm  35-
      subroutine cool (data,nprev,n,nrem,isign,ifact,work)              coo   1
c     fourier transform of length n.  in place cooley-tukey method,     coo   2
c     digit-reversed to normal order, sande-tukey factoring (2).        coo   3
c     dimension data(nprev,n,nrem)                                      coo   4
c     complex data                                                      coo   5
c     data(i1,j2,i3) = sum(data(i1,i2,i3)*exp(isign*2*pi*i*((i2-1)*     coo   6
c     (j2-1)/n))), summed over i2 = 1 to n for all i1 from 1 to nprev,  coo   7
c     j2 from 1 to n and i3 from 1 to nrem.  the factors of n are given coo   8
c     in any order in array ifact.  factors of two are done in pairs    coo   9
c     as much as possible (fourier transform of length four), factors ofcoo  10
c     three are done separately, and all factors five or higher         coo  11
c     are done by goertzel's algorithm (4).                             coo  12
      dimension data(*), work(*), ifact(1)                              coo  13
      twopi=6.283185307*float(isign)                                    coo  14
      ip0=2                                                             coo  15
      ip1=ip0*nprev                                                     coo  16
      ip4=ip1*n                                                         coo  17
      ip5=ip4*nrem                                                      coo  18
      if=0                                                              coo  19
      ip2=ip1                                                           coo  20
 10   if (ip2-ip4) 20,240,240                                           coo  21
 20   if=if+1                                                           coo  22
      ifcur=ifact(if)                                                   coo  23
      if (ifcur-2) 60,30,60                                             coo  24
 30   if (4*ip2-ip4) 40,40,60                                           coo  25
 40   if (ifact(if+1)-2) 60,50,60                                       coo  26
 50   if=if+1                                                           coo  27
      ifcur=4                                                           coo  28
 60   ip3=ip2*ifcur                                                     coo  29
      theta=twopi/float(ifcur)                                          coo  30
      sinth=sin(theta/2.)                                               coo  31
      rootr=-2.*sinth*sinth                                             coo  32
c     cos(theta)-1, for accuracy.                                       coo  33
      rooti=sin(theta)                                                  coo  34
      theta=twopi/float(ip3/ip1)                                        coo  35
      sinth=sin(theta/2.)                                               coo  36
      wstpr=-2.*sinth*sinth                                             coo  37
      wstpi=sin(theta)                                                  coo  38
      wr=1.                                                             coo  39
      wi=0.                                                             coo  40
      do 230 i2=1,ip2,ip1                                               coo  41
      if (ifcur-4) 70,70,210                                            coo  42
 70   if ((i2-1)*(ifcur-2)) 240,90,80                                   coo  43
 80   w2r=wr*wr-wi*wi                                                   coo  44
      w2i=2.*wr*wi                                                      coo  45
      w3r=w2r*wr-w2i*wi                                                 coo  46
      w3i=w2r*wi+w2i*wr                                                 coo  47
 90   i1max=i2+ip1-ip0                                                  coo  48
      do 200 i1=i2,i1max,ip0                                            coo  49
      do 200 i5=i1,ip5,ip3                                              coo  50
      j0=i5                                                             coo  51
      j1=j0+ip2                                                         coo  52
      j2=j1+ip2                                                         coo  53
      j3=j2+ip2                                                         coo  54
      if (i2-1) 140,140,100                                             coo  55
 100  if (ifcur-3) 130,120,110                                          coo  56
c     apply the phase shift factors                                     coo  57
 110  tempr=data(j3)                                                    coo  58
      data(j3)=w3r*tempr-w3i*data(j3+1)                                 coo  59
      data(j3+1)=w3r*data(j3+1)+w3i*tempr                               coo  60
      tempr=data(j2)                                                    coo  61
      data(j2)=wr*tempr-wi*data(j2+1)                                   coo  62
      data(j2+1)=wr*data(j2+1)+wi*tempr                                 coo  63
      tempr=data(j1)                                                    coo  64
      data(j1)=w2r*tempr-w2i*data(j1+1)                                 coo  65
      data(j1+1)=w2r*data(j1+1)+w2i*tempr                               coo  66
      go to 140                                                         coo  67
 120  tempr=data(j2)                                                    coo  68
      data(j2)=w2r*tempr-w2i*data(j2+1)                                 coo  69
      data(j2+1)=w2r*data(j2+1)+w2i*tempr                               coo  70
 130  tempr=data(j1)                                                    coo  71
      data(j1)=wr*tempr-wi*data(j1+1)                                   coo  72
      data(j1+1)=wr*data(j1+1)+wi*tempr                                 coo  73
 140  if (ifcur-3) 150,160,170                                          coo  74
c     do a fourier transform of length two                              coo  75
 150  tempr=data(j1)                                                    coo  76
      tempi=data(j1+1)                                                  coo  77
      data(j1)=data(j0)-tempr                                           coo  78
      data(j1+1)=data(j0+1)-tempi                                       coo  79
      data(j0)=data(j0)+tempr                                           coo  80
      data(j0+1)=data(j0+1)+tempi                                       coo  81
      go to 200                                                         coo  82
c     do a fourier transform of length three                            coo  83
 160  sumr=data(j1)+data(j2)                                            coo  84
      sumi=data(j1+1)+data(j2+1)                                        coo  85
      tempr=data(j0)-.5*sumr                                            coo  86
      tempi=data(j0+1)-.5*sumi                                          coo  87
      data(j0)=data(j0)+sumr                                            coo  88
      data(j0+1)=data(j0+1)+sumi                                        coo  89
      difr=rooti*(data(j2+1)-data(j1+1))                                coo  90
      difi=rooti*(data(j1)-data(j2))                                    coo  91
      data(j1)=tempr+difr                                               coo  92
      data(j1+1)=tempi+difi                                             coo  93
      data(j2)=tempr-difr                                               coo  94
      data(j2+1)=tempi-difi                                             coo  95
      go to 200                                                         coo  96
c     do a fourier transform of length four (from bit reversed order)   coo  97
 170  t0r=data(j0)+data(j1)                                             coo  98
      t0i=data(j0+1)+data(j1+1)                                         coo  99
      t1r=data(j0)-data(j1)                                             coo 100
      t1i=data(j0+1)-data(j1+1)                                         coo 101
      t2r=data(j2)+data(j3)                                             coo 102
      t2i=data(j2+1)+data(j3+1)                                         coo 103
      t3r=data(j2)-data(j3)                                             coo 104
      t3i=data(j2+1)-data(j3+1)                                         coo 105
      data(j0)=t0r+t2r                                                  coo 106
      data(j0+1)=t0i+t2i                                                coo 107
      data(j2)=t0r-t2r                                                  coo 108
      data(j2+1)=t0i-t2i                                                coo 109
      if (isign) 180,180,190                                            coo 110
 180  t3r=-t3r                                                          coo 111
      t3i=-t3i                                                          coo 112
 190  data(j1)=t1r-t3i                                                  coo 113
      data(j1+1)=t1i+t3r                                                coo 114
      data(j3)=t1r+t3i                                                  coo 115
      data(j3+1)=t1i-t3r                                                coo 116
 200  continue                                                          coo 117
      go to 220                                                         coo 118
c     do a fourier transform of length five or more                     coo 119
 210  call goert (data(i2),nprev,ip2/ip1,ifcur,ip5/ip3,work,wr,wi,rootr,coo 120
     $rooti)                                                            coo 121
 220  tempr=wr                                                          coo 122
      wr=wstpr*tempr-wstpi*wi+tempr                                     coo 123
 230  wi=wstpr*wi+wstpi*tempr+wi                                        coo 124
      ip2=ip3                                                           coo 125
      go to 10                                                          coo 126
 240  return                                                            coo 127
      end                                                               coo 128-
      subroutine factr (n,ifact,nfact)                                  fac   1
c     factor n into its prime factors, nfact in number.  for example,   fac   2
c     for n = 1960, nfact = 6 and ifact(if) = 2, 2, 2, 5, 7 and 7.      fac   3
      dimension ifact(1)                                                fac   4
      if=0                                                              fac   5
      npart=n                                                           fac   6
      do 50 id=1,n,2                                                    fac   7
      idiv=id                                                           fac   8
      if (id-1) 10,10,20                                                fac   9
 10   idiv=2                                                            fac  10
 20   iquot=npart/idiv                                                  fac  11
      if (npart-idiv*iquot) 40,30,40                                    fac  12
 30   if=if+1                                                           fac  13
      ifact(if)=idiv                                                    fac  14
      npart=iquot                                                       fac  15
      go to 20                                                          fac  16
 40   if (iquot-idiv) 60,60,50                                          fac  17
 50   continue                                                          fac  18
 60   if (npart-1) 80,80,70                                             fac  19
 70   if=if+1                                                           fac  20
      ifact(if)=npart                                                   fac  21
 80   nfact=if                                                          fac  22
      return                                                            fac  23
      end                                                               fac  24-
      subroutine fixrl (data,n,nrem,isign,iform)                        fix   1
c     for iform = 0, convert the transform of a doubled-up real array,  fix   2
c     considered complex, into its true transform.  supply only the     fix   3
c     first half of the complex transform, as the second half has       fix   4
c     conjugate symmetry.  for iform = -1, convert the first half       fix   5
c     of the true transform into the transform of a doubled-up real     fix   6
c     array.  n must be even.                                           fix   7
c     using complex notation and subscripts starting at zero, the       fix   8
c     transformation is--                                               fix   9
c     dimension data(n,nrem)                                            fix  10
c     zstp = exp(isign*2*pi*i/n)                                        fix  11
c     do 10 i2=0,nrem-1                                                 fix  12
c     data(0,i2) = conj(data(0,i2))*(1+i)                               fix  13
c     do 10 i1=1,n/4                                                    fix  14
c     z = (1+(2*iform+1)*i*zstp**i1)/2                                  fix  15
c     i1cnj = n/2-i1                                                    fix  16
c     dif = data(i1,i2)-conj(data(i1cnj,i2))                            fix  17
c     temp = z*dif                                                      fix  18
c     data(i1,i2) = (data(i1,i2)-temp)*(1-iform)                        fix  19
c 10  data(i1cnj,i2) = (data(i1cnj,i2)+conj(temp))*(1-iform)            fix  20
c     if i1=i1cnj, the calculation for that value collapses into        fix  21
c     a simple conjugation of data(i1,i2).                              fix  22
      dimension data(*)                                                 fix  23
      twopi=6.283185307*float(isign)                                    fix  24
      ip0=2                                                             fix  25
      ip1=ip0*(n/2)                                                     fix  26
      ip2=ip1*nrem                                                      fix  27
      if (iform) 10,70,70                                               fix  28
c     pack the real input values (two per column)                       fix  29
 10   j1=ip1+1                                                          fix  30
      data(2)=data(j1)                                                  fix  31
      if (nrem-1) 70,70,20                                              fix  32
 20   j1=j1+ip0                                                         fix  33
      i2min=ip1+1                                                       fix  34
      do 60 i2=i2min,ip2,ip1                                            fix  35
      data(i2)=data(j1)                                                 fix  36
      j1=j1+ip0                                                         fix  37
      if (n-2) 50,50,30                                                 fix  38
 30   i1min=i2+ip0                                                      fix  39
      i1max=i2+ip1-ip0                                                  fix  40
      do 40 i1=i1min,i1max,ip0                                          fix  41
      data(i1)=data(j1)                                                 fix  42
      data(i1+1)=data(j1+1)                                             fix  43
 40   j1=j1+ip0                                                         fix  44
 50   data(i2+1)=data(j1)                                               fix  45
 60   j1=j1+ip0                                                         fix  46
 70   do 80 i2=1,ip2,ip1                                                fix  47
      tempr=data(i2)                                                    fix  48
      data(i2)=data(i2)+data(i2+1)                                      fix  49
 80   data(i2+1)=tempr-data(i2+1)                                       fix  50
      if (n-2) 200,200,90                                               fix  51
 90   theta=twopi/float(n)                                              fix  52
      sinth=sin(theta/2.)                                               fix  53
      zstpr=-2.*sinth*sinth                                             fix  54
      zstpi=sin(theta)                                                  fix  55
      zr=(1.-zstpi)/2.                                                  fix  56
      zi=(1.+zstpr)/2.                                                  fix  57
      if (iform) 100,110,110                                            fix  58
 100  zr=1.-zr                                                          fix  59
      zi=-zi                                                            fix  60
 110  i1min=ip0+1                                                       fix  61
      i1max=ip0*(n/4)+1                                                 fix  62
      do 190 i1=i1min,i1max,ip0                                         fix  63
      do 180 i2=i1,ip2,ip1                                              fix  64
      i2cnj=ip0*(n/2+1)-2*i1+i2                                         fix  65
      if (i2-i2cnj) 150,120,120                                         fix  66
 120  if (isign*(2*iform+1)) 130,140,140                                fix  67
 130  data(i2+1)=-data(i2+1)                                            fix  68
 140  if (iform) 170,180,180                                            fix  69
 150  difr=data(i2)-data(i2cnj)                                         fix  70
      difi=data(i2+1)+data(i2cnj+1)                                     fix  71
      tempr=difr*zr-difi*zi                                             fix  72
      tempi=difr*zi+difi*zr                                             fix  73
      data(i2)=data(i2)-tempr                                           fix  74
      data(i2+1)=data(i2+1)-tempi                                       fix  75
      data(i2cnj)=data(i2cnj)+tempr                                     fix  76
      data(i2cnj+1)=data(i2cnj+1)-tempi                                 fix  77
      if (iform) 160,180,180                                            fix  78
 160  data(i2cnj)=data(i2cnj)+data(i2cnj)                               fix  79
      data(i2cnj+1)=data(i2cnj+1)+data(i2cnj+1)                         fix  80
 170  data(i2)=data(i2)+data(i2)                                        fix  81
      data(i2+1)=data(i2+1)+data(i2+1)                                  fix  82
 180  continue                                                          fix  83
      tempr=zr-.5                                                       fix  84
      zr=zstpr*tempr-zstpi*zi+zr                                        fix  85
 190  zi=zstpr*zi+zstpi*tempr+zi                                        fix  86
c     recursion saves time, at a slight loss in accuracy.  if available,fix  87
c     use double precision to compute zr and zi.                        fix  88
 200  if (iform) 270,210,210                                            fix  89
c     unpack the real transform values (two per column)                 fix  90
 210  i2=ip2+1                                                          fix  91
      i1=i2                                                             fix  92
      j1=ip0*(n/2+1)*nrem+1                                             fix  93
      go to 250                                                         fix  94
 220  data(j1)=data(i1)                                                 fix  95
      data(j1+1)=data(i1+1)                                             fix  96
      i1=i1-ip0                                                         fix  97
      j1=j1-ip0                                                         fix  98
 230  if (i2-i1) 220,240,240                                            fix  99
 240  data(j1)=data(i1)                                                 fix 100
      data(j1+1)=0.                                                     fix 101
 250  i2=i2-ip1                                                         fix 102
      j1=j1-ip0                                                         fix 103
      data(j1)=data(i2+1)                                               fix 104
      data(j1+1)=0.                                                     fix 105
      i1=i1-ip0                                                         fix 106
      j1=j1-ip0                                                         fix 107
      if (i2-1) 260,260,230                                             fix 108
 260  data(2)=0.                                                        fix 109
 270  return                                                            fix 110
      end                                                               fix 111-
      subroutine goert(data,nprev,iprod,ifact,irem,work,wminr,wmini,    goe   1
     $ rootr,rooti)                                                     goe   2
c     phase-shifted fourier transform of length ifact by the goertzel   goe   3
c     algorithm (4).  ifact must be odd and at least 5.  further speed  goe   4
c     is gained by computing two transform values at the same time.     goe   5
c     dimension data(nprev,iprod,ifact,irem)                            goe   6
c     data(i1,1,j3,i5) = sum(data(i1,1,i3,i5) * w**(i3-1)), summed      goe   7
c     over i3 = 1 to ifact for all i1 from 1 to nprev, j3 from 1 to     goe   8
c     ifact and i5 from 1 to irem.                                      goe   9
c     w = wmin * exp(isign*2*pi*i*(j3-1)/ifact).                        goe  10
      dimension data(*), work(*)                                        goe  11
      ip0=2                                                             goe  12
      ip1=ip0*nprev                                                     goe  13
      ip2=ip1*iprod                                                     goe  14
      ip3=ip2*ifact                                                     goe  15
      ip5=ip3*irem                                                      goe  16
      if (wmini) 10,40,10                                               goe  17
c     apply the phase shift factors                                     goe  18
 10   wr=wminr                                                          goe  19
      wi=wmini                                                          goe  20
      i3min=1+ip2                                                       goe  21
      do 30 i3=i3min,ip3,ip2                                            goe  22
      i1max=i3+ip1-ip0                                                  goe  23
      do 20 i1=i3,i1max,ip0                                             goe  24
      do 20 i5=i1,ip5,ip3                                               goe  25
      tempr=data(i5)                                                    goe  26
      data(i5)=wr*tempr-wi*data(i5+1)                                   goe  27
 20   data(i5+1)=wr*data(i5+1)+wi*tempr                                 goe  28
      tempr=wr                                                          goe  29
      wr=wminr*tempr-wmini*wi                                           goe  30
 30   wi=wminr*wi+wmini*tempr                                           goe  31
 40   do 90 i1=1,ip1,ip0                                                goe  32
      do 90 i5=i1,ip5,ip3                                               goe  33
c     straight summation for the first term                             goe  34
      sumr=0.                                                           goe  35
      sumi=0.                                                           goe  36
      i3max=i5+ip3-ip2                                                  goe  37
      do 50 i3=i5,i3max,ip2                                             goe  38
      sumr=sumr+data(i3)                                                goe  39
 50   sumi=sumi+data(i3+1)                                              goe  40
      work(1)=sumr                                                      goe  41
      work(2)=sumi                                                      goe  42
      wr=rootr+1.                                                       goe  43
      wi=rooti                                                          goe  44
      iwmin=1+ip0                                                       goe  45
      iwmax=ip0*((ifact+1)/2)-1                                         goe  46
      do 80 iwork=iwmin,iwmax,ip0                                       goe  47
      twowr=wr+wr                                                       goe  48
      i3=i3max                                                          goe  49
      oldsr=0.                                                          goe  50
      oldsi=0.                                                          goe  51
      sumr=data(i3)                                                     goe  52
      sumi=data(i3+1)                                                   goe  53
      i3=i3-ip2                                                         goe  54
 60   tempr=sumr                                                        goe  55
      tempi=sumi                                                        goe  56
      sumr=twowr*sumr-oldsr+data(i3)                                    goe  57
      sumi=twowr*sumi-oldsi+data(i3+1)                                  goe  58
      oldsr=tempr                                                       goe  59
      oldsi=tempi                                                       goe  60
      i3=i3-ip2                                                         goe  61
      if (i3-i5) 70,70,60                                               goe  62
c     in a fourier transform the w corresponding to the point at k      goe  63
c     is the conjugate of that at ifact-k (that is, exp(twopi*i*        goe  64
c     k/ifact) = conj(exp(twopi*i*(ifact-k)/ifact))).  since the        goe  65
c     main loop of goertzels algorithm is indifferent to the imaginary  goe  66
c     part of w, it need be supplied only at the end.                   goe  67
 70   tempr=-wi*sumi                                                    goe  68
      tempi=wi*sumr                                                     goe  69
      sumr=wr*sumr-oldsr+data(i3)                                       goe  70
      sumi=wr*sumi-oldsi+data(i3+1)                                     goe  71
      work(iwork)=sumr+tempr                                            goe  72
      work(iwork+1)=sumi+tempi                                          goe  73
      iwcnj=ip0*(ifact+1)-iwork                                         goe  74
      work(iwcnj)=sumr-tempr                                            goe  75
      work(iwcnj+1)=sumi-tempi                                          goe  76
c     singleton's recursion, for accuracy and speed (5).                goe  77
      tempr=wr                                                          goe  78
      wr=wr*rootr-wi*rooti+wr                                           goe  79
 80   wi=tempr*rooti+wi*rootr+wi                                        goe  80
      iwork=1                                                           goe  81
      do 90 i3=i5,i3max,ip2                                             goe  82
      data(i3)=work(iwork)                                              goe  83
      data(i3+1)=work(iwork+1)                                          goe  84
 90   iwork=iwork+ip0                                                   goe  85
      return                                                            goe  86
      end                                                               goe  87-
      subroutine smfac (ifact,nfact,isym,ifsym,nfsym,icent,ifcnt,nfcnt) smf   1
c     rearrange the prime factors of n into a square and a non-         smf   2
c     square.  n = isym*icent*isym, where icent is square-free.         smf   3
c     isym = ifsym(1)*...*ifsym(nfsym), each a prime factor.            smf   4
c     icent = ifcnt(1)*...*ifcnt(nfcnt), each a prime factor.           smf   5
c     for example, n = 1960 = 14*10*14.  then isym = 14, icent = 10,    smf   6
c     nfsym = 2, nfcnt = 2, nfact = 6, ifsym(ifs) = 2, 7, ifcnt(ifc) =  smf   7
c     2, 5 and ifact(if) = 2, 7, 2, 5, 7, 2.                            smf   8
      dimension ifsym(1), ifcnt(1), ifact(1)                            smf   9
      isym=1                                                            smf  10
      icent=1                                                           smf  11
      ifs=0                                                             smf  12
      ifc=0                                                             smf  13
      if=1                                                              smf  14
 10   if (if-nfact) 20,40,50                                            smf  15
 20   if (ifact(if)-ifact(if+1)) 40,30,40                               smf  16
 30   ifs=ifs+1                                                         smf  17
      ifsym(ifs)=ifact(if)                                              smf  18
      isym=ifact(if)*isym                                               smf  19
      if=if+2                                                           smf  20
      go to 10                                                          smf  21
 40   ifc=ifc+1                                                         smf  22
      ifcnt(ifc)=ifact(if)                                              smf  23
      icent=ifact(if)*icent                                             smf  24
      if=if+1                                                           smf  25
      go to 10                                                          smf  26
 50   nfsym=ifs                                                         smf  27
      nfcnt=ifc                                                         smf  28
      nfsm2=2*nfsym                                                     smf  29
      nfact=2*nfsym+nfcnt                                               smf  30
      if (nfcnt) 80,80,60                                               smf  31
 60   nfsm2=nfsm2+1                                                     smf  32
      ifsym(nfsym+1)=icent                                              smf  33
      do 70 ifc=1,nfcnt                                                 smf  34
      if=nfsym+ifc                                                      smf  35
 70   ifact(if)=ifcnt(ifc)                                              smf  36
 80   if (nfsym) 110,110,90                                             smf  37
 90   do 100 ifs=1,nfsym                                                smf  38
      ifscj=nfsm2+1-ifs                                                 smf  39
      ifsym(ifscj)=ifsym(ifs)                                           smf  40
      ifact(ifs)=ifsym(ifs)                                             smf  41
      ifcnj=nfact+1-ifs                                                 smf  42
 100  ifact(ifcnj)=ifsym(ifs)                                           smf  43
 110  nfsym=nfsm2                                                       smf  44
      return                                                            smf  45
      end                                                               smf  46-
      subroutine symrv (data,nprev,n,nrem,ifact,nfact)                  sym   1
c     shuffle the data array by reversing the digits of one index.      sym   2
c     dimension data(nprev,n,nrem)                                      sym   3
c     replace data(i1,i2,i3) by data(i1,i2rev,i3) for all i1 from 1 to  sym   4
c     nprev, i2 from 1 to n and i3 from 1 to nrem.  i2rev-1 is the      sym   5
c     integer whose digit representation in the multi-radix notation    sym   6
c     of factors ifact(if) is the reverse of the representation of i2-1.sym   7
c     for example, if all ifact(if) = 2, i2-1 = 11001, i2rev-1 = 10011. sym   8
c     the factors must be symmetrically arranged, i.e., ifact(if) =     sym   9
c     ifact(nfact+1-if).                                                sym  10
      dimension data(*), ifact(1)                                       sym  11
      if (nfact-1) 80,80,10                                             sym  12
 10   ip0=2                                                             sym  13
      ip1=ip0*nprev                                                     sym  14
      ip4=ip1*n                                                         sym  15
      ip5=ip4*nrem                                                      sym  16
      i4rev=1                                                           sym  17
      do 70 i4=1,ip4,ip1                                                sym  18
      if (i4-i4rev) 20,40,40                                            sym  19
 20   i1max=i4+ip1-ip0                                                  sym  20
      do 30 i1=i4,i1max,ip0                                             sym  21
      do 30 i5=i1,ip5,ip4                                               sym  22
      i5rev=i4rev+i5-i4                                                 sym  23
      tempr=data(i5)
      tempi=data(i5+1)                                                  sym  25
      data(i5)=data(i5rev)                                              sym  26
      data(i5+1)=data(i5rev+1)                                          sym  27
      data(i5rev)=tempr                                                 sym  28
 30   data(i5rev+1)=tempi                                               sym  29
 40   ip3=ip4                                                           sym  30
      do 60 if=1,nfact                                                  sym  31
      ip2=ip3/ifact(if)                                                 sym  32
      i4rev=i4rev+ip2                                                   sym  33
      if (i4rev-ip3) 70,70,50                                           sym  34
 50   i4rev=i4rev-ip3                                                   sym  35
 60   ip3=ip2                                                           sym  36
 70   continue                                                          sym  37
 80   return                                                            sym  38
      end                                                               sym  39-
