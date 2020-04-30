	program tish
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ************** tish.f ****************
c Computation of SH synthetic seismograms 
c in transversely isotropic media for anisotropic PREM
c using modified DSM operators & modified source representation.
c Synthetics for shallow events can be computed.
c
c                                                 2002.10 K.Kawai
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ----------------------------<<constants>>----------------------------
	implicit none
	real*8 pi,lmaxdivf,shallowdepth
	integer maxnlay
	integer maxnzone,maxnr,maxlmax,ilog
	parameter ( pi=3.1415926535897932d0 )
	parameter ( maxnlay = 88300 )
	parameter ( maxnzone = 15 )
	parameter ( maxnr = 600 )
	parameter ( maxlmax = 80000 )
	parameter ( ilog = 0 )
	parameter ( lmaxdivf = 2.d4)
	parameter ( shallowdepth = 100.d0 )
c ----------------------------<<variables>>----------------------------
c variable for the trial function
	integer nnlayer,nlayer(maxnzone)
	integer l,m
	real*8 ra(maxnlay+maxnzone+1),gra(3),plm(3,0:3,maxnr)
	complex*16 bvec(3,-2:2,maxnr)
c variable for the structure
	integer nzone
	integer ndc,vnp
 	real*8 rmin,rmax
	real*8 vrmin(maxnzone),vrmax(maxnzone)
	real*8 rrho(4,maxnzone),vsv(4,maxnzone),vsh(4,maxnzone)
	real*8 qmu(maxnzone)
	real*8 vra(maxnlay+2*maxnzone+1)
	real*8 rho(maxnlay+2*maxnzone+1)
	real*8 ecL(maxnlay+2*maxnzone+1)
	real*8 ecN(maxnlay+2*maxnzone+1)
	real*8 gvra(3),grho(3),gecL(3),gecN(3)
	complex*16 coef(maxnzone)
c variable for the periodic range
	integer np,imin,imax
	real*8 tlen,omega,omegai
	complex*16 u(3,maxnr)
c variable for the source
	integer spn,ns
	real*8 r0,mt(3,3),spo,mu0,eqlat,eqlon
c variable for the station
	integer nr,ir
	real*8 theta(maxnr),phi(maxnr)
	real*8 lat(maxnr),lon(maxnr)
c variable for the matrix elements
	complex*16 a0( 2,maxnlay+1 ), a2( 2,maxnlay+1 )
	complex*16  a( 2,maxnlay+1 )
	real*8 t( 4*maxnlay )
	real*8 h1( 4*maxnlay ),h2( 4*maxnlay )
        real*8 h3( 4*maxnlay ),h4( 4*maxnlay )
	real*8 gt(8),gh1(8),gh2(8),gh3(8),gh4(8)
	complex*16 aa(4),ga(8),ga2(2,3),gdr(3)
	complex*16 g( maxnlay+1 )
c variable for the file
	character*80 output(maxnr)
c variable for grid spacing
	real*8 tmpr(maxnlay+1)
	real*8 gridpar(maxnzone),dzpar(maxnzone),vmin(maxnzone)
	real*8 re,ratc,ratl,maxamp
	integer kc,lsuf,ismall,llog
c variable for the stack point
	integer isp(maxnzone),jsp(maxnzone),ins
c other variables
	integer i,j,ii,jj,nn,lda,ier
	real*8 eps,work( 4*maxnlay ),lsq
	complex*16 dr(maxnlay+1),z(maxnlay+1)
	complex*16 cwork( 4*maxnlay )
	integer ltmp(2),iimax
c
	data lda/ 2 /
	data eps/ -1.d0 /
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *************** Inputting and computing the parameters ***************
c --- inputting parameter ---
	call pinput2( maxnlay,maxnzone,maxnr,re,ratc,ratl,
     &       tlen,np,omegai,imin,imax,
     &       nzone,vrmin,vrmax,rrho,vsv,vsh,qmu,
     &       r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output )
c --- computing the required parameters ---
c computing and checking the parameters
	rmin = vrmin(1)
	rmax = vrmax(nzone)
	ndc = nzone - 1
	do ir=1,nr
           theta(ir)= theta(ir) / 1.8d2 * pi
           phi(ir)= phi(ir)   / 1.8d2 * pi
        enddo
	if ( (r0.lt.rmin) .or. (r0.gt.rmax) )
     &       pause 'Location of the source is improper.'
c ************************** Files Handling **************************
	do ir=1,nr
	   open( unit=11,file=output(ir),status='unknown')
	   write(11,*) 'VERSION172'
	   write(11,*) tlen
	   write(11,*) np
	   write(11,*) omegai,lat(ir),lon(ir)
c         write(11,*) theta(ir)*1.8d2/pi,phi(ir)*1.8d2/pi
	   write(11,*) eqlat,eqlon,r0
	   close(11)
        enddo
	if(ilog.eq.1) then
	   open(unit=11,file='llog.log',status='unknown')
	   close(11)
	endif
	iimax = imax
	if( (rmax-r0).lt.shallowdepth) then ! option for shallow events
c computing of the number and the location of grid points
	   iimax = int(tlen * 2.d0)
	   call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax,
     &       iimax,1,tlen,
     &       vmin,gridpar,dzpar )
	   call calra ( maxnlay,maxnzone,
     &       nnlayer,
     &       gridpar,dzpar,nzone,vrmin,vrmax,
     &       rmin,rmax,nlayer,ra,re )
c --- checking the parameter
	   if ( nnlayer.gt.maxnlay )
     &       pause 'The number of grid points is too large.'
c computing the stack points
	   call calsp( ndc,nlayer,isp,jsp )
c computing the source location
	   call calspo( ndc,vrmax,nnlayer,r0,rmin,rmax,ra,
     &       isp,spo,spn )
c computing grids for source computations
	   call calgra( isp,ra,r0,spn,spo,gra )
c ******************* Computing the matrix elements *******************
c computing the structure grid points
	   call calstg( nzone,rrho,vsv,vsh,
     &       nnlayer,nlayer,ra,rmax,
     &       vnp,vra,rho,ecL,ecN)
	   call calgstg( spn,rrho,vsv,vsh,
     &       gra,gvra,rmax,grho,gecL,gecN,r0,mu0 )
	   do i=1,ndc+1
	      call calmatc( nlayer(i),vnp,vra,rho,2,0,0,
     &          ra( isp(i) ),t( jsp(i) ),work( jsp(i) ) )
	      call calmatc( nlayer(i),vnp,vra,ecL ,2,1,1,
     &          ra( isp(i) ),h1( jsp(i) ),work( jsp(i) ) )
	      call calmatc( nlayer(i),vnp,vra,ecL ,1,1,0,
     &          ra( isp(i) ),h2( jsp(i) ),work( jsp(i) ) )
	      call calmatc( nlayer(i),vnp,vra,ecL ,0,0,0,
     &          ra( isp(i) ),h3( jsp(i) ),work( jsp(i) ) )
	      call calmatc( nlayer(i),vnp,vra,ecN ,0,0,0,
     &          ra( isp(i) ),h4( jsp(i) ),work( jsp(i) ) )
	      call caltl( nlayer(i),vnp,vra,rho,
     &          ra( isp(i) ),work( jsp(i) ) )
	      call calt( nlayer(i),  t( jsp(i) ),  work( jsp(i) ),
     &          t( jsp(i) ) )
	      call calhl( nlayer(i),vnp,vra,ecL,
     &          ra( isp(i) ),work( jsp(i) ) )
	      call calt( nlayer(i), h3( jsp(i) ), work( jsp(i) ),
     &          h3( jsp(i) ) )
	      call calhl( nlayer(i),vnp,vra,ecN,
     &          ra( isp(i) ),work( jsp(i) ) )
	      call calt( nlayer(i), h4( jsp(i) ), work( jsp(i) ),
     &          h4( jsp(i) ) )
	   enddo
	   call calmatc( 2,3,gvra,grho,2,0,0,gra,gt, work )
	   call calmatc( 2,3,gvra,gecL ,2,1,1,gra,gh1,work )
	   call calmatc( 2,3,gvra,gecL ,1,1,0,gra,gh2,work )
	   call calmatc( 2,3,gvra,gecL ,0,0,0,gra,gh3,work )
	   call calmatc( 2,3,gvra,gecN ,0,0,0,gra,gh4,work )
	   call caltl( 2,3,gvra,grho,gra,work )
	   call calt( 2,  gt, work, gt )
	   call calhl( 2,3,gvra,gecL, gra,work )
	   call calt( 2, gh3, work, gh3 )
	   call calhl( 2,3,gvra,gecN, gra,work )
	   call calt( 2, gh4, work, gh4 )
c     
	   nn = nnlayer + 1
	   ns = isp(spn) + dint(spo)
	   ins = 4 * ns - 3
c
	   llog = 0
	   do ii=1,2	! omega-loop
	      if(ii.eq.1) then
		 if(imin.eq.0) then
		    i=1
		 else
		    i=imin
		 endif
	      endif
	      if(ii.eq.2) i=imax
              omega = 2.d0 * pi * dble(i) / tlen
              call callsuf(omega,nzone,vrmax,vsv,lsuf)
              call calcoef( nzone,omega,qmu,coef )
c
              call cmatinit( lda,nn,a0 )
              call cmatinit( lda,nn,a2 )
              do j=1,ndc+1
                 call cala0( nlayer(j),omega,omegai,
     &                t(jsp(j)), h1(jsp(j)),
     &                h2(jsp(j)), h3(jsp(j)),
     &                h4(jsp(j)),
     &                coef(j), cwork(jsp(j)) )
                 call overlap( nlayer(j),cwork(jsp(j)),
     &                a0( 1,isp(j) ) )
                 call cala2( nlayer(j),h4(jsp(j)),
     &                coef(j), cwork(jsp(j)) )
                 call overlap( nlayer(j),cwork(jsp(j)),
     &                a2( 1,isp(j) ) )
              enddo
c
              kc = 1
              ismall = 0
              maxamp = -1.d0
              ltmp(ii) = maxlmax
              do l=0,maxlmax	! l-loop
                 if( ismall.gt.20 ) then
                    if(ltmp(ii).gt.l) ltmp(ii) = l
                    exit
                 endif
c     
                 do jj=1,maxnlay+1 ! initialize
                    tmpr(jj) = 0.d0
                 enddo
                 lsq = dsqrt( dble(l)*dble(l+1) )
c computing the coefficient matrix elements
c --- Initializing the matrix elements
                 call cmatinit( lda,nn,a )
                 call cmatinit( lda,3,ga2 )
                 call cala( nn,l,lda,a0,a2,a )
                 call calga( 1,omega,omegai,l,
     &                t(ins),h1(ins),h2(ins),h3(ins),h4(ins),
     &                coef(spn),aa )
                 call calga( 2,omega,omegai,l,gt,gh1,gh2,gh3,gh4,
     &                coef(spn),ga )
                 call overlap( 2,ga,ga2 )
c
                 do m=-2,2	! m-loop
                    if ( ( m.ne.0 ).and.( iabs(m).le.iabs(l) ) ) then
                       call cvecinit( nn,g )
                       call calg2( l,m,spo,r0,mt,mu0,coef(spn),
     &                      ga,aa,ga2,gdr,g( isp(spn) ) )
                       if( mod(l,100).eq.0) then
                          if ( (m.eq.-2).or.(m.eq.-l) ) then
                             call dclisb0( a,nn,1,lda,g,eps,dr,z,ier)
                          else
                             call dcsbsub0( a,nn,1,lda,g,eps,dr,z,ier)
                          endif
                          do jj=1,nn !sum up c of the same l
                             tmpr(jj) = tmpr(jj) + cdabs(g(jj))
                          enddo
                       else
                          if ( (m.eq.-2).or.(m.eq.-l) ) then
                             call dclisb( a(1,kc),nn-kc+1,1,lda,ns-kc+1
     &                            ,g(kc),eps,dr,z,ier)
                          else
                             call dcsbsub( a(1,kc),nn-kc+1,1,lda,ns-kc+1
     &                            ,g(kc),eps,dr,z,ier)
                          endif
                       endif
c     
                       if( mod(l,100).eq.0) then
                          call calcutd(nzone,nlayer,tmpr,ratc,nn,ra,kc)
                       endif
c     
                       call calamp(g(nn),l,lsuf,maxamp,ismall,ratl)
                    endif
                 enddo          ! m-loop
              enddo             ! l-loop
	   enddo		! omega-loop
	   iimax = dble(max(ltmp(1),ltmp(2))) * tlen / lmaxdivf
	endif			! option for shallow events
c     
c computing of the number and the location of grid points
	call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax,
     &       iimax,1,tlen,
     &       vmin,gridpar,dzpar )
	call calra ( maxnlay,maxnzone,
     &       nnlayer,
     &       gridpar,dzpar,nzone,vrmin,vrmax,
     &       rmin,rmax,nlayer,ra,re )
c --- checking the parameter
	if ( nnlayer.gt.maxnlay )
     &       pause 'The number of grid points is too large.'
c computing the stack points
	call calsp( ndc,nlayer,isp,jsp )
c computing the source location
	call calspo( ndc,vrmax,nnlayer,r0,rmin,rmax,ra,
     &       isp,spo,spn )
c computing grids for source computations
	call calgra( isp,ra,r0,spn,spo,gra )
c ******************* Computing the matrix elements *******************
c computing the structure grid points
	call calstg( nzone,rrho,vsv,vsh,
     &       nnlayer,nlayer,ra,rmax,
     &       vnp,vra,rho,ecL,ecN)
	call calgstg( spn,rrho,vsv,vsh,
     &       gra,gvra,rmax,grho,gecL,gecN,r0,mu0 )
	do i=1,ndc+1
           call calmatc( nlayer(i),vnp,vra,rho,2,0,0,
     &          ra( isp(i) ),t( jsp(i) ),work( jsp(i) ) )
           call calmatc( nlayer(i),vnp,vra,ecL ,2,1,1,
     &          ra( isp(i) ),h1( jsp(i) ),work( jsp(i) ) )
           call calmatc( nlayer(i),vnp,vra,ecL ,1,1,0,
     &          ra( isp(i) ),h2( jsp(i) ),work( jsp(i) ) )
           call calmatc( nlayer(i),vnp,vra,ecL ,0,0,0,
     &          ra( isp(i) ),h3( jsp(i) ),work( jsp(i) ) )
           call calmatc( nlayer(i),vnp,vra,ecN ,0,0,0,
     &          ra( isp(i) ),h4( jsp(i) ),work( jsp(i) ) )
           call caltl( nlayer(i),vnp,vra,rho,
     &          ra( isp(i) ),work( jsp(i) ) )
           call calt( nlayer(i),  t( jsp(i) ),  work( jsp(i) ),
     &          t( jsp(i) ) )
           call calhl( nlayer(i),vnp,vra,ecL,
     &          ra( isp(i) ),work( jsp(i) ) )
           call calt( nlayer(i), h3( jsp(i) ), work( jsp(i) ),
     &          h3( jsp(i) ) )
           call calhl( nlayer(i),vnp,vra,ecN,
     &          ra( isp(i) ),work( jsp(i) ) )
           call calt( nlayer(i), h4( jsp(i) ), work( jsp(i) ),
     &          h4( jsp(i) ) )
        enddo
	call calmatc( 2,3,gvra,grho,2,0,0,gra,gt, work )
	call calmatc( 2,3,gvra,gecL ,2,1,1,gra,gh1,work )
	call calmatc( 2,3,gvra,gecL ,1,1,0,gra,gh2,work )
	call calmatc( 2,3,gvra,gecL ,0,0,0,gra,gh3,work )
	call calmatc( 2,3,gvra,gecN ,0,0,0,gra,gh4,work )
	call caltl( 2,3,gvra,grho,gra,work )
	call calt( 2,  gt, work, gt )
	call calhl( 2,3,gvra,gecL, gra,work )
	call calt( 2, gh3, work, gh3 )
	call calhl( 2,3,gvra,gecN, gra,work )
	call calt( 2, gh4, work, gh4 )
c     
c ******************** Computing the displacement *********************
	nn = nnlayer + 1
	ns = isp(spn) + dint(spo)
	ins = 4 * ns - 3
c
	llog = 0
	do i=imin,imax		! omega-loop
           call cmatinit( 3,nr,u )
           if ( i.ne.0 ) then
              omega = 2.d0 * pi * dble(i) / tlen
              call callsuf(omega,nzone,vrmax,vsv,lsuf)
              do ir=1,nr
                 call matinit( 3,4,plm(1,0,ir) )
              enddo
c	    if ( lmin(ipband).gt.0 ) then
c	      do 160 l=0,lmin(ipband)-1
c	        do 150 ir=1,nr
c	          call calbvec( l,theta(ir),phi(ir),
c     &	                        plm(1,0,ir),bvec(1,-2,ir) )
c  150	        continue
c  160	      continue
c	    endif
              call calcoef( nzone,omega,qmu,coef )
c
              call cmatinit( lda,nn,a0 )
              call cmatinit( lda,nn,a2 )
              do j=1,ndc+1
                 call cala0( nlayer(j),omega,omegai,
     &                t(jsp(j)), h1(jsp(j)),
     &                h2(jsp(j)), h3(jsp(j)),
     &                h4(jsp(j)),
     &                coef(j), cwork(jsp(j)) )
                 call overlap( nlayer(j),cwork(jsp(j)),
     &                a0( 1,isp(j) ) )
                 call cala2( nlayer(j),h4(jsp(j)),
     &                coef(j), cwork(jsp(j)) )
                 call overlap( nlayer(j),cwork(jsp(j)),
     &                a2( 1,isp(j) ) )
              enddo
c
              kc = 1
              ismall = 0
              maxamp = -1.d0
              llog = maxlmax
              do l=0,maxlmax	! l-loop
                 if( ismall.gt.20 ) then
                    if(llog.gt.l) llog = l
                    cycle
                 endif
c     
                 do jj=1,maxnlay+1 ! initialize
                    tmpr(jj) = 0.d0
                 enddo
                 lsq = dsqrt( dble(l)*dble(l+1) )
c ***** Computing the trial function *****
                 do ir=1,nr
                    call calbvec( l,theta(ir),phi(ir),
     &                   plm(1,0,ir),bvec(1,-2,ir) )
                 enddo
c computing the coefficient matrix elements
c --- Initializing the matrix elements
                 call cmatinit( lda,nn,a )
                 call cmatinit( lda,3,ga2 )
                 call cala( nn,l,lda,a0,a2,a )
                 call calga( 1,omega,omegai,l,
     &                t(ins),h1(ins),h2(ins),h3(ins),h4(ins),
     &                coef(spn),aa )
c	      call calga2( 2,omega,omegai,l,gt,gh1,gh2,gh3,
c     &                    coef(spn),ga )
                 call calga( 2,omega,omegai,l,gt,gh1,gh2,gh3,gh4,
     &                coef(spn),ga )
                 call overlap( 2,ga,ga2 )
c
                 do m=-2,2	! m-loop
                    if ( ( m.ne.0 ).and.( iabs(m).le.iabs(l) ) ) then
                       call cvecinit( nn,g )
                       call calg2( l,m,spo,r0,mt,mu0,coef(spn),
     &                      ga,aa,ga2,gdr,g( isp(spn) ) )
                       if( mod(l,100).eq.0) then
                          if ( (m.eq.-2).or.(m.eq.-l) ) then
                             call dclisb0( a,nn,1,lda,g,eps,dr,z,ier)
                          else
                             call dcsbsub0( a,nn,1,lda,g,eps,dr,z,ier)
                          endif
                          do jj=1,nn !sum up c of the same l
                             tmpr(jj) = tmpr(jj) + cdabs(g(jj))
                          enddo
                       else
                          if ( (m.eq.-2).or.(m.eq.-l) ) then
                             call dclisb( a(1,kc),nn-kc+1,1,lda,ns-kc+1
     &                            ,g(kc),eps,dr,z,ier)
                          else
                             call dcsbsub( a(1,kc),nn-kc+1,1,lda,ns-kc+1
     &                            ,g(kc),eps,dr,z,ier)
                          endif
                       endif
c     
                       if( mod(l,100).eq.0) then
                          call calcutd(nzone,nlayer,tmpr,ratc,nn,ra,kc)
                       endif
c     
                       call calamp(g(nn),l,lsuf,maxamp,ismall,ratl)
                       do ir=1,nr
                          call calu( g(nn),lsq,bvec(1,m,ir),u(1,ir) )
                       enddo
                    endif
                 enddo          ! m-loop
              enddo             ! l-loop
           endif
c ************************** Files Handling **************************
           do ir=1,nr
              open( unit=11,file=output(ir),
     &             position='append',status='old')
              write(11,*) i,dble(u(1,ir)),dimag(u(1,ir))
              write(11,*) dble(u(2,ir)),dimag(u(2,ir))
              write(11,*) dble(u(3,ir)),dimag(u(3,ir))
              close(11)
           enddo
           if(ilog.eq.1) then
              open(unit=11,file='llog.log',
     &         position='append',status='old')
              write(11,*) i,llog,nnlayer
              close(11)
           endif
        enddo                   ! omega-loop
c     
c
	end
