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
	parameter ( maxnzone = 35 )
	parameter ( maxnr = 10000 )
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
	complex*16 :: u(3,maxnr)
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
       real*8 t1i,t1f
c
	data lda/ 2 /
	data eps/ -1.d0 /
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **************MPI***********************************
	include 'mpif.h'
	integer :: mpii
	integer :: petot, my_rank, ierr
	integer :: filenum, mpios
	integer :: outputmemory ! MB
	integer :: outputinterval
	real(8) :: memoryperomega ! MB
	integer :: outputindex
	integer, allocatable, dimension (:) :: outputi
	complex*16, allocatable, dimension(:,:,:) :: outputu
c       when the values to be output use memory over outputmemory MB,
c       they are written in output files. The interval is outputinterval
c       memoryperomega is the quantity of memory used for one omega step
	integer, allocatable, dimension (:) :: mpimin, mpimax
	character *2 :: char_rank
	double precision :: ark, angel
	data outputmemory /300/
	call mpi_init (ierr)
	call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
	call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
c	write(*,*) "myrank",my_rank
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	allocate(mpimin(PETOT), mpimax(PETOT))
cccccccccccccccccccccccccccccccccccccccccccccccc
	filenum = 10+my_rank
	write (char_rank,'(I1)') my_rank 
c	memoryperomega = 3*16 * maxnr*0.000001
c	outputinterval = outputmemory/memoryperomega != integer * nr
c	allocate (outputu(outputinterval*3))
c	allocate (outputi(outputinterval))
c	call mpi_finalize(ierr)
c	stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *************** Inputting and computing the parameters ***************
c --- inputting parameter ---
	if (my_rank .eq.0) then
	   call pinput2( maxnlay,maxnzone,maxnr,re,ratc,ratl,
	1	tlen,np,omegai,imin,imax,
	2	nzone,vrmin,vrmax,rrho,vsv,vsh,qmu,
	3	r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)
	endif
	call MPI_BCAST (re, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (ratc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (ratl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (tlen, 1, MPI_DOUBLE_PRECISION, 0,
	1    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (np, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (omegai, 1, MPI_DOUBLE_PRECISION, 0,
	1    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (imin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (imax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (nzone, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (vrmin, maxnzone, MPI_DOUBLE_PRECISION, 0,
	1    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (vrmax, maxnzone, MPI_DOUBLE_PRECISION, 0,
	1    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (rrho, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,
	2    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (vsv, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,
	3    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (vsh, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,
	4    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (qmu, maxnzone, MPI_DOUBLE_PRECISION, 0,
	1    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (r0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,
	2    ierr)
	call MPI_BCAST (eqlat, 1, MPI_DOUBLE_PRECISION, 0,
	3    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (eqlon, 1, MPI_DOUBLE_PRECISION, 0,
	4    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (mt, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,
	5    ierr)
	call MPI_BCAST (nr, 1, MPI_INTEGER, 0,
	1    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (theta, maxnr, MPI_DOUBLE_PRECISION, 0,
	2    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (phi, maxnr, MPI_DOUBLE_PRECISION, 0,
	3    MPI_COMM_WORLD, ierr)
	call MPI_BCAST (lat(1), maxnr, MPI_DOUBLE_PRECISION,
	1    0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (lon(1), maxnr, MPI_DOUBLE_PRECISION,
	2    0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST (output, 80*maxnr, MPI_CHARACTER,
	3    0, MPI_COMM_WORLD, ierr)
ccccccccccccccccccccccccccccccccccccccccccccc
	memoryperomega = 1*3*16 *nr*0.000001
	outputinterval = outputmemory/memoryperomega != integer * nr
c	write (*,*) outputinterval
c	call mpi_finalize(ierr)
c	stop
	allocate (outputi(outputinterval))
	allocate (outputu(3,nr,outputinterval))
c	do i=1,3
c	   do j = 1,nr
c	      do ii = 1, outputinterval
c		 outputu(i,j,ii)=cmplx(1,0)
c	      enddo
c	   enddo
c	enddo
	
ccccccccccccccccccccccccccccccccccccccccccccccc
	
c --- computing the required parameters ---
c computing and checking the parameters
	rmin = vrmin(1)
	rmax = vrmax(nzone)
	ndc = nzone - 1
c	write(*,*) my_rank  ! debug
c	call mpi_finalize (ierr)
c	stop
	do ir=1,nr
           theta(ir)= theta(ir) / 1.8d2 * pi
           phi(ir)= phi(ir)   / 1.8d2 * pi
        enddo
	if ( (r0.lt.rmin) .or. (r0.gt.rmax) )
     &       pause 'Location of the source is improper.'
c ************************** Files Handling **************************
	if (my_rank.eq.0)then
	   do ir=1,nr
c	   open( unit=filenum,file=trim(output(ir))//"."//char_rank,
c	1	status='unknown')
c	   close(filenum)
	      open( 10,file=trim(output(ir)),
     $      status='unknown',form='binary',
     $         convert='big_endian')
	      write(10) tlen
	      write(10) np,1,3
	      write(10) omegai,lat(ir),lon(ir)
c       write(11,*) theta(ir)*1.8d2/pi,phi(ir)*1.8d2/pi
	      write(10) eqlat,eqlon,r0
	      close (10)
	   enddo
	endif



c	if(ilog.eq.1) then
c	   open(unit=filenum,file='llog.log.'//char_rank,
c	1	status='unknown')
c	   close(filenum)
c	endif
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

ccccccccc MPIccccccccccccccccccccccccccc


c	write(*,*) imin, "  ", imax, "  ",iamari
	
c	call simplesplit (imin, imax, PETOT, mpimin, mpimax)
	call trianglesplit (imin, imax, PETOT, mpimin, mpimax)

	outputindex = 1
c	write (*,*) my_rank
c	write (*,*) mpimin(my_rank+1), mpimax(my_rank+1)
        t1i = MPI_Wtime()
	do i= mpimin(my_rank+1), mpimax(my_rank+1)! omega-loop
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
c              t1i = MPI_Wtime()
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
c              t1f = MPI_Wtime()
c              write (*,*) "One l-loop took",t1f-t1i,"s",i
           endif
c ************************** Files Handling **************************
	   
c	   filenum = 11+my_rank
c	   write (char_rank,'(I1)') my_rank 
c       outputinterval
	   outputi(outputindex)=i
           do ir=1,nr
c	      call mpi_finalize(ierr)
c	      stop
c              open( unit=filenum,file=trim(output(ir))//"."//char_rank,
c     &             position='append',status='old')
c              write(filenum,*) i,dble(u(1,ir)),dimag(u(1,ir))
c              write(filenum,*) dble(u(2,ir)),dimag(u(2,ir))
c              write(filenum,*) dble(u(3,ir)),dimag(u(3,ir))
c              close(filenum)
c	      outputi (outputindex)=i
	      outputu (1,ir,outputindex) =u (1,ir)
	      outputu (2, ir, outputindex)= u (2,ir)	 
c	      outputu (2, ir, outputindex)=
c	1	   cmplx(dble(u (2,ir)), dimag(u(2,ir)))
	      outputu (3, ir, outputindex)= u(3,ir)
c	      write (*,*) outputu (1,ir,outputindex)
c	      write(*,*) dimag(outputu (1,ir,outputindex)
c	      write(*,*)outputu (2,ir,outputindex) 
c	      write(*,*)u (2,ir)
c	      write (*,*) dimag(outputu(2,ir,outputindex))
c	      write(*,*) dimag(u(2,ir))
           enddo

	   if (outputindex .ge. outputinterval .or.
	1	i .eq. mpimax(my_rank+1)) then
c	      write(*,*) "kakikomimasu"
	      do ir = 1 ,nr
 120		continue       
                open(unit=filenum,file=trim(output(ir)),
     &     position='append',status='unknown',
     &          form='binary', convert='big_endian',
     &                 err=100, share='denyrw' )
        goto 555
 100    continue       
c       write(*,*) my_rank, "waiting for 100s", output(ir)
        call system("sleep 100")
        goto 120
 555                   continue
		 do mpii= 1, outputindex
		    write(filenum) outputi(mpii),
	1		 dble(outputu(1,ir,mpii)),
	1		 dimag(outputu(1,ir,mpii))
		    write(filenum) dble(outputu(2,ir,mpii)),
	1		 dimag(outputu(2,ir,mpii))
		    write(filenum) dble(outputu(3,ir,mpii)),
	2		 dimag(outputu(3,ir,mpii))
		 enddo
		 close(filenum)
	      enddo
	      outputindex =0
	   endif
	   outputindex = outputindex + 1
           if(ilog.eq.1) then
              open(unit=filenum,file='llog.log.'//char_rank,
     &         position='append',status='old')
              write(filenum,*) i,llog,nnlayer
              close(filenum)
           endif
	   

c	   call mpi_finalize(ierr)
c	   stop
	

        enddo                   ! omega-loop
        t1f = MPI_Wtime()
        write(*,*) "Run finished in ",t1f-t1i,"s"
c     
c
c	call mpi_finalize(ierr)
c	stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	write(*,*) my_rank, "Ivalice looks to the horizon!"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	call mpi_finalize (ierr)
	stop
	end
