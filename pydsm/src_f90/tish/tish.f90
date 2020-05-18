subroutine tish(re,ratc,ratl,tlen,np,omegai,imin,imax, &
   nzone,vrmin,vrmax,rrho,vsv,vsh,qmu, &
   r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output,write_to_file, &
   outputu)
!------------------------------------------------------------------------
!  ************** tish.f ****************
! Computation of SH synthetic seismograms
! in transversely isotropic media for anisotropic PREM
! using modified DSM operators & modified source representation.
! Synthetics for shallow events can be computed.
!
!                                                 2002.10 K.Kawai
!------------------------------------------------------------------------
      use parameters
      implicit none

! ----------------------------<<variables>>----------------------------
! write to files?
      logical, intent(in) :: write_to_file
! variable for the trial function
      integer nnlayer,nlayer(maxnzone)
      integer l,m
      real(dp) ra(maxnlay+maxnzone+1),gra(3),plm(3,0:3,maxnr)
      complex(dp) bvec(3,-2:2,maxnr)
! variable for the structure
      integer nzone
      integer ndc,vnp
      real(dp) rmin,rmax
      real(dp) vrmin(maxnzone),vrmax(maxnzone)
      real(dp) rrho(4,maxnzone),vsv(4,maxnzone),vsh(4,maxnzone)
      real(dp) qmu(maxnzone)
      real(dp) vra(maxnlay+2*maxnzone+1)
      real(dp) rho(maxnlay+2*maxnzone+1)
      real(dp) ecL(maxnlay+2*maxnzone+1)
      real(dp) ecN(maxnlay+2*maxnzone+1)
      real(dp) gvra(3),grho(3),gecL(3),gecN(3)
      complex(dp) coef(maxnzone)
! variable for the periodic range
      integer np,imin,imax
      real(dp) tlen,omega,omegai
      complex(dp) u(3,maxnr)
! variable for the source
      integer spn,ns
      real(dp) r0,mt(3,3),spo,mu0,eqlat,eqlon
! variable for the station
      integer nr,ir
      real(dp) theta(maxnr),phi(maxnr)
      real(dp) lat(maxnr),lon(maxnr)
! variable for the matrix elements
      complex(dp) a0( 2,maxnlay+1 ), a2( 2,maxnlay+1 )
      complex(dp)  a( 2,maxnlay+1 )
      real(dp) t( 4*maxnlay )
      real(dp) h1( 4*maxnlay ),h2( 4*maxnlay )
            real(dp) h3( 4*maxnlay ),h4( 4*maxnlay )
      real(dp) gt(8),gh1(8),gh2(8),gh3(8),gh4(8)
      complex(dp) aa(4),ga(8),ga2(2,3),gdr(3)
      complex(dp) g( maxnlay+1 )
! variable for the file
      character*80 :: output(maxnr)
! variable for grid spacing
      real(dp) tmpr(maxnlay+1)
      real(dp) gridpar(maxnzone),dzpar(maxnzone),vmin(maxnzone)
      real(dp) re,ratc,ratl,maxamp
      integer kc,lsuf,ismall,llog
! variable for the stack point
      integer isp(maxnzone),jsp(maxnzone),ins
! other variables
      integer i,j,ii,jj,nn,lda,ier
      real(dp) eps,work( 4*maxnlay ),lsq
      complex(dp) dr(maxnlay+1),z(maxnlay+1)
      complex(dp) cwork( 4*maxnlay )
      integer ltmp(2),iimax

      data lda/ 2 /
      data eps/ -1.d0 /

      integer :: mpii
      complex(dp), intent(out) :: outputu(3,nr,imin:imax)

! --- computing the required parameters ---
! computing and checking the parameters
      rmin = vrmin(1)
      rmax = vrmax(nzone)
      ndc = nzone - 1
      do ir=1,nr
         theta(ir)= theta(ir) / 1.8d2 * pi
         phi(ir)= phi(ir) / 1.8d2 * pi
      enddo
      if ( (r0 < rmin) .or. (r0 > rmax) ) then
         write(*,*) 'Location of the source is improper.'
         return
      endif
! ************************** Files Handling **************************
      if (write_to_file) then
         do ir=1,nr
            open(unit=11,file=output(ir),status='replace', &
               access='stream', form='unformatted', convert='big_endian')
            write(11) tlen
            write(11) np,1,3
            write(11) omegai,lat(ir),lon(ir)
   !        write(11,*) theta(ir)*1.8d2/pi,phi(ir)*1.8d2/pi
            write(11) eqlat,eqlon,r0
            close(11)
         enddo
      endif

      if(ilog.eq.1) then
         open(unit=11,file='llog.log',status='unknown')
         close(11)
      endif

      iimax = imax
      if( (rmax-r0) < shallowdepth) then ! option for shallow events
! computing of the number and the location of grid points
         iimax = int(tlen * 2.d0)
         call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax, &
            iimax,1,tlen,vmin,gridpar,dzpar )
         call calra ( maxnlay,maxnzone, &
            nnlayer,gridpar,dzpar,nzone,vrmin,vrmax, &
            rmin,rmax,nlayer,ra,re )
! --- checking the parameter
         if ( nnlayer.gt.maxnlay ) then
            write(*,*) 'The number of grid points is too large.'
            return
         endif
! computing the stack points
         call calsp( ndc,nlayer,isp,jsp )
! computing the source location
         call calspo( ndc,vrmax,nnlayer,r0,rmin,rmax,ra,isp,spo,spn )
! computing grids for source computations
         call calgra( isp,ra,r0,spn,spo,gra )
! ******************* Computing the matrix elements *******************
! computing the structure grid points
         call calstg( nzone,rrho,vsv,vsh, nnlayer,nlayer,ra,rmax, &
            vnp,vra,rho,ecL,ecN)
         call calgstg( spn,rrho,vsv,vsh,gra,gvra,rmax,grho,gecL,gecN,r0,mu0 )
         do i=1,ndc+1
	         call calmatc( nlayer(i),vnp,vra,rho,2,0,0, &
               ra( isp(i) ),t( jsp(i) ),work( jsp(i) ) )
	         call calmatc( nlayer(i),vnp,vra,ecL ,2,1,1, &
               ra( isp(i) ),h1( jsp(i) ),work( jsp(i) ) )
	         call calmatc( nlayer(i),vnp,vra,ecL ,1,1,0, &
               ra( isp(i) ),h2( jsp(i) ),work( jsp(i) ) )
	         call calmatc( nlayer(i),vnp,vra,ecL ,0,0,0, &
               ra( isp(i) ),h3( jsp(i) ),work( jsp(i) ) )
	         call calmatc( nlayer(i),vnp,vra,ecN ,0,0,0, &
               ra( isp(i) ),h4( jsp(i) ),work( jsp(i) ) )
	         call caltl( nlayer(i),vnp,vra,rho, &
               ra( isp(i) ),work( jsp(i) ) )
	         call calt( nlayer(i),  t( jsp(i) ),  work( jsp(i) ), &
               t( jsp(i) ) )
	         call calhl( nlayer(i),vnp,vra,ecL, &
               ra( isp(i) ),work( jsp(i) ) )
	         call calt( nlayer(i), h3( jsp(i) ), work( jsp(i) ), &
               h3( jsp(i) ) )
	         call calhl( nlayer(i),vnp,vra,ecN, &
               ra( isp(i) ),work( jsp(i) ) )
	         call calt( nlayer(i), h4( jsp(i) ), work( jsp(i) ), &
               h4( jsp(i) ) )
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
!
         nn = nnlayer + 1
         ns = isp(spn) + dint(spo)
         ins = 4 * ns - 3
!
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
!
            call cmatinit( lda,nn,a0 )
            call cmatinit( lda,nn,a2 )
            do j=1,ndc+1
               call cala0( nlayer(j),omega,omegai, &
                     t(jsp(j)), h1(jsp(j)),h2(jsp(j)), h3(jsp(j)), &
                     h4(jsp(j)),coef(j), cwork(jsp(j)) )
               call overlap( nlayer(j),cwork(jsp(j)),a0( 1,isp(j) ) )
               call cala2( nlayer(j),h4(jsp(j)),coef(j), cwork(jsp(j)) )
               call overlap( nlayer(j),cwork(jsp(j)),a2( 1,isp(j) ) )
            enddo
!
            kc = 1
            ismall = 0
            maxamp = -1.d0
            ltmp(ii) = maxlmax
            do l=0,maxlmax	! l-loop
               if( ismall.gt.20 ) then
                  if(ltmp(ii).gt.l) ltmp(ii) = l
                  exit
               endif
!
               do jj=1,maxnlay+1 ! initialize
                  tmpr(jj) = 0.d0
               enddo
               lsq = dsqrt( dble(l)*dble(l+1) )
! computing the coefficient matrix elements
! --- Initializing the matrix elements
               call cmatinit( lda,nn,a )
               call cmatinit( lda,3,ga2 )
               call cala( nn,l,lda,a0,a2,a )
               call calga( 1,omega,omegai,l, &
                     t(ins),h1(ins),h2(ins),h3(ins),h4(ins), &
                     coef(spn),aa )
               call calga( 2,omega,omegai,l,gt,gh1,gh2,gh3,gh4, &
                     coef(spn),ga )
               call overlap( 2,ga,ga2 )
!
               do m=-2,2	! m-loop
                  if ( ( m.ne.0 ).and.( iabs(m).le.iabs(l) ) ) then
                     call cvecinit( nn,g )
                     call calg2( l,m,spo,r0,mt,mu0,coef(spn), &
                           ga,aa,ga2,gdr,g( isp(spn) ) )
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
                           call dclisb( a(1,kc),nn-kc+1,1,lda,ns-kc+1 &
                                 ,g(kc),eps,dr,z,ier)
                        else
                           call dcsbsub( a(1,kc),nn-kc+1,1,lda,ns-kc+1 &
                                 ,g(kc),eps,dr,z,ier)
                        endif
                     endif
!
                     if( mod(l,100).eq.0) then
                        call calcutd(nzone,nlayer,tmpr,ratc,nn,ra,kc)
                     endif
!
                     call calamp(g(nn),l,lsuf,maxamp,ismall,ratl)
                  endif
               enddo          ! m-loop
            enddo             ! l-loop
         enddo		! omega-loop
         iimax = dble(max(ltmp(1),ltmp(2))) * tlen / lmaxdivf
      endif			! option for shallow events
!
! computing of the number and the location of grid points
      call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax, &
            iimax,1,tlen,vmin,gridpar,dzpar )
      call calra ( maxnlay,maxnzone, &
            nnlayer,gridpar,dzpar,nzone,vrmin,vrmax, &
            rmin,rmax,nlayer,ra,re )
! --- checking the parameter
      if ( nnlayer.gt.maxnlay ) then
         write(*,*) 'The number of grid points is too large.'
         return
      endif
! computing the stack points
      call calsp( ndc,nlayer,isp,jsp )
! computing the source location
      call calspo( ndc,vrmax,nnlayer,r0,rmin,rmax,ra,isp,spo,spn )
! computing grids for source computations
      call calgra( isp,ra,r0,spn,spo,gra )
! ******************* Computing the matrix elements *******************
! computing the structure grid points
      call calstg( nzone,rrho,vsv,vsh, &
            nnlayer,nlayer,ra,rmax,vnp,vra,rho,ecL,ecN)
      call calgstg( spn,rrho,vsv,vsh, &
            gra,gvra,rmax,grho,gecL,gecN,r0,mu0 )
      do i=1,ndc+1
           call calmatc( nlayer(i),vnp,vra,rho,2,0,0, &
               ra( isp(i) ),t( jsp(i) ),work( jsp(i) ) )
           call calmatc( nlayer(i),vnp,vra,ecL ,2,1,1, &
               ra( isp(i) ),h1( jsp(i) ),work( jsp(i) ) )
           call calmatc( nlayer(i),vnp,vra,ecL ,1,1,0, &
               ra( isp(i) ),h2( jsp(i) ),work( jsp(i) ) )
           call calmatc( nlayer(i),vnp,vra,ecL ,0,0,0, &
               ra( isp(i) ),h3( jsp(i) ),work( jsp(i) ) )
           call calmatc( nlayer(i),vnp,vra,ecN ,0,0,0, &
               ra( isp(i) ),h4( jsp(i) ),work( jsp(i) ) )
           call caltl( nlayer(i),vnp,vra,rho, &
               ra( isp(i) ),work( jsp(i) ) )
           call calt( nlayer(i),  t( jsp(i) ),  work( jsp(i) ), &
               t( jsp(i) ) )
           call calhl( nlayer(i),vnp,vra,ecL, &
               ra( isp(i) ),work( jsp(i) ) )
           call calt( nlayer(i), h3( jsp(i) ), work( jsp(i) ), &
               h3( jsp(i) ) )
           call calhl( nlayer(i),vnp,vra,ecN, &
               ra( isp(i) ),work( jsp(i) ) )
           call calt( nlayer(i), h4( jsp(i) ), work( jsp(i) ), &
               h4( jsp(i) ) )
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
!
! ******************** Computing the displacement *********************
      nn = nnlayer + 1
      ns = isp(spn) + dint(spo)
      ins = 4 * ns - 3
!
      llog = 0
      do i=imin,imax		! omega-loop
         call cmatinit( 3,nr,u )
         if ( i.ne.0 ) then
            omega = 2.d0 * pi * dble(i) / tlen
            call callsuf(omega,nzone,vrmax,vsv,lsuf)
            do ir=1,nr
               call matinit( 3,4,plm(1,0,ir) )
            enddo
!	    if ( lmin(ipband).gt.0 ) then
!	      do 160 l=0,lmin(ipband)-1
!	        do 150 ir=1,nr
!	          call calbvec( l,theta(ir),phi(ir),
!     &	                        plm(1,0,ir),bvec(1,-2,ir) )
!  150	        continue
!  160	      continue
!	    endif
            call calcoef( nzone,omega,qmu,coef )
!
            call cmatinit( lda,nn,a0 )
            call cmatinit( lda,nn,a2 )
            do j=1,ndc+1
               call cala0( nlayer(j),omega,omegai, &
                     t(jsp(j)), h1(jsp(j)),h2(jsp(j)), h3(jsp(j)), &
                     h4(jsp(j)),coef(j), cwork(jsp(j)) )
               call overlap( nlayer(j),cwork(jsp(j)),a0( 1,isp(j) ) )
               call cala2( nlayer(j),h4(jsp(j)),coef(j), cwork(jsp(j)) )
               call overlap( nlayer(j),cwork(jsp(j)),a2( 1,isp(j) ) )
            enddo
!
            kc = 1
            ismall = 0
            maxamp = -1.d0
            llog = maxlmax
            do l=0,maxlmax	! l-loop
               if( ismall.gt.20 ) then
                  if(llog.gt.l) llog = l
                  cycle
               endif
!
               do jj=1,maxnlay+1 ! initialize
                  tmpr(jj) = 0.d0
               enddo
               lsq = dsqrt( dble(l)*dble(l+1) )
! ***** Computing the trial function *****
               do ir=1,nr
                  call calbvec( l,theta(ir),phi(ir),plm(1,0,ir),bvec(1,-2,ir) )
               enddo
! computing the coefficient matrix elements
! --- Initializing the matrix elements
               call cmatinit( lda,nn,a )
               call cmatinit( lda,3,ga2 )
               call cala( nn,l,lda,a0,a2,a )
               call calga( 1,omega,omegai,l, &
                  t(ins),h1(ins),h2(ins),h3(ins),h4(ins),coef(spn),aa )
!	      call calga2( 2,omega,omegai,l,gt,gh1,gh2,gh3,
!     &                    coef(spn),ga )
               call calga( 2,omega,omegai,l,gt,gh1,gh2,gh3,gh4,coef(spn),ga )
               call overlap( 2,ga,ga2 )
!
               do m=-2,2	! m-loop
                  if ( ( m.ne.0 ).and.( iabs(m).le.iabs(l) ) ) then
                     call cvecinit( nn,g )
                     call calg2( l,m,spo,r0,mt,mu0,coef(spn), &
                        ga,aa,ga2,gdr,g( isp(spn) ) )
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
                           call dclisb( a(1,kc),nn-kc+1,1,lda,ns-kc+1 &
                            ,g(kc),eps,dr,z,ier)
                        else
                           call dcsbsub( a(1,kc),nn-kc+1,1,lda,ns-kc+1 &
                                 ,g(kc),eps,dr,z,ier)
                        endif
                     endif
!
                     if( mod(l,100).eq.0) then
                        call calcutd(nzone,nlayer,tmpr,ratc,nn,ra,kc)
                     endif
!
                     call calamp(g(nn),l,lsuf,maxamp,ismall,ratl)
                     do ir=1,nr
                        call calu( g(nn),lsq,bvec(1,m,ir),u(1,ir) )
                     enddo

                  endif
               enddo          ! m-loop
            enddo             ! l-loop
         endif
! ************************** Files Handling **************************
         do ir=1,nr
            outputu(1,ir,i) = u(1,ir)
            outputu(2,ir,i) = u(2,ir)
	         outputu(3,ir,i) = u(3,ir)

         enddo
         if(ilog.eq.1) then
            open(unit=11,file='llog.log',position='append',status='old')
            write(11,*) i,llog,nnlayer
            close(11)
         endif
!
   	   if (write_to_file .and. i .eq. imax) then
	         write(*,*) "kakikomimasu"
	         do ir = 1 ,nr
               open(unit=10,file=output(ir),position='append',status='old', &
                  access='stream', form='unformatted', convert='big_endian')
               do mpii= imin, imax
                  write(10) mpii,dble(outputu(1,ir,mpii)), &
                     dimag(outputu(1,ir,mpii))
		            write(10) dble(outputu(2,ir,mpii)),dimag(outputu(2,ir,mpii))
		            write(10) dble(outputu(3,ir,mpii)),dimag(outputu(3,ir,mpii))
               enddo
               close(10)
            enddo
         endif
      enddo                   ! omega-loop
!
!
      !write(*,*) "Ivalice looks to the horizon"

      return

   end subroutine tish

program main

   use parameters
   implicit none

   character(len=160) :: parameter_file

   integer np
   integer imin
   integer imax
   integer nzone
   integer nr
   real(dp) :: tlen
   real(dp) :: omegai
   real(dp) :: re
   real(dp) :: ratc
   real(dp) :: ratl
   real(dp) :: vrmin(maxnzone)
   real(dp) :: vrmax(maxnzone)
   real(dp) :: qmu(maxnzone)
   real(dp) :: rho(4,maxnzone)
   real(dp) :: vsv(4,maxnzone)
   real(dp) :: vsh(4,maxnzone)
   real(dp) :: r0
   real(dp) :: mt(3,3)
   real(dp) :: lat(maxnr)
   real(dp) :: lon(maxnr)
   real(dp) :: theta(maxnr)
   real(dp) :: phi(maxnr)
   real(dp) :: eqlat
   real(dp) :: eqlon
   character(len=80) :: output(maxnr)
   complex(dp), allocatable, dimension(:,:,:) :: outputu
   logical :: write_to_file = .true.

   ! read input parameters
   call get_command_argument(1, parameter_file)
   call pinput_fromfile(parameter_file, &
      re,ratc,ratl,tlen,np,omegai,imin,imax, &
      nzone,vrmin,vrmax,rho,vsv,vsh,qmu, &
      r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)

   allocate(outputu(3,nr,imin:imax))

   ! main loop
   write(*,*) 'Enter main loop'
   call tish(re,ratc,ratl,tlen,np,omegai,imin,imax, &
      nzone,vrmin,vrmax,rho,vsv,vsh,qmu, &
      r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output,write_to_file,outputu)
   write(*,*) 'Done!'

end program main
