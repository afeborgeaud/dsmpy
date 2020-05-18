subroutine tipsv(re,ratc,ratl,tlen,np,omegai,imin,imax, &
   nzone,vrmin,vrmax,rrho,vpv,vph,vsv,vsh,eta,qmu,qkappa, &
   r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output,write_to_file,&
   outputu)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  ************** tipsv.f ****************
!c Computation of PSV synthetic seismograms
!c in transversely isotropic for anisotropic PREM
!c using modified DSM operators & modified source representation.
!c Synthetics for shallow events can be computed.
!c                                                 2002.12  K.Kawai
!c -> f90 Kensuke Konishi 2020 May 5
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use parameters
    implicit none
!c ---------------------------<< constants >>---------------------------
! write to files?
    logical, intent(in) :: write_to_file
!c ---------------------------<< variables >>---------------------------
!c variable for the trial function
    integer nlayer(maxnzone)
    integer nslay,nllay
    integer inlayer,jnlayer,jnslay,jnllay
    integer l,m
    double precision ra(maxnlay+maxnzone+1),plm(3,0:3,maxnr)
    complex(dp) bvec(3,-2:2,maxnr)
!c variable for the structure
    integer nzone,isl,ill,nsl,nll
    integer iphase(maxnzone),ndc,vnp
    double precision rmin,rmax
    double precision,intent(in):: vrmin(maxnzone),vrmax(maxnzone)
    double precision,intent(in):: rrho(4,maxnzone), vpv(4,maxnzone),vph(4,maxnzone)
    double precision,intent(in):: vsv(4,maxnzone),vsh(4,maxnzone),eta(4,maxnzone)
    double precision,intent(in):: qmu(maxnzone),qkappa(maxnzone)
    double precision vra(maxnlay+2*maxnzone+1)
    double precision rho(maxnlay+2*maxnzone+1)
    double precision kappa(maxnlay+2*maxnzone+1)
    double precision ecKx(maxnlay+2*maxnzone+1) !3*Kx=3A-4N
    double precision ecKy(maxnlay+2*maxnzone+1) !3*Ky=3F+2N
    double precision ecKz(maxnlay+2*maxnzone+1) !3*Kz=2F+C
    double precision mu(maxnlay+2*maxnzone+1)
    double precision ecL(maxnlay+2*maxnzone+1)
    double precision ecN(maxnlay+2*maxnzone+1)
    double precision rhoinv(maxnlay+2*maxnzone+1)
    double precision kappainv(maxnlay+2*maxnzone+1)
    complex(dp) coef1(maxnzone),coef2(maxnzone),coef(maxnzone)
!c variable for the periodic range
    integer np,imin,imax
    double precision tlen,omega,omegai
    complex(dp) u(3,maxnr)
!c variable for the source
    integer spn,ns
    double precision,intent(in):: r0,mt(3,3),eqlat,eqlon
    double precision:: spo,ecC0,ecF0,ecL0
    complex(dp) ya(4),yb(4),yc(4),yd(4)
!c variable for the station
    integer,intent(in):: nr
    double precision:: theta(maxnr),phi(maxnr)
    double precision,intent(in):: lat(maxnr),lon(maxnr)
!c variable for the matrix elements
    complex(dp) a0(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
    complex(dp) a1(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
    complex(dp) a2(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
    complex(dp) a(4,2*(maxnslay+1) + (maxnllay+1) )
    complex(dp) c(2, (maxnslay+1) + (maxnllay+1) )
    double precision t(8*maxnslay )
    double precision h1x(8*maxnslay ),h1y( 8*maxnslay ),h1z( 8*maxnslay )
    double precision h2L(8*maxnslay ),h2N( 8*maxnslay )
    double precision h3ax(8*maxnslay )
    double precision h3ay(8*maxnslay ),h3az( 8*maxnslay )
    double precision h4aL(8*maxnslay ),h4aN( 8*maxnslay )
    double precision h5ax(8*maxnslay ),h5ay( 8*maxnslay ),h5az( 8*maxnslay )
    double precision h6aL(8*maxnslay ),h6aN( 8*maxnslay )
    double precision h3x(8*maxnslay ),h3y( 8*maxnslay ),h3z( 8*maxnslay )
    double precision h4L(8*maxnslay ),h4N( 8*maxnslay )
    double precision h5x(8*maxnslay ),h5y( 8*maxnslay ),h5z( 8*maxnslay )
    double precision h6L(8*maxnslay ),h6N( 8*maxnslay )
    double precision h7x(8*maxnslay ),h7y( 8*maxnslay ),h7z( 8*maxnslay )
    double precision h8L(8*maxnslay ),h8N( 8*maxnslay )
    double precision h3mx(-2:1,2*(maxnslay+maxnzone) )
    double precision h3my(-2:1,2*(maxnslay+maxnzone) )
    double precision h3mz(-2:1,2*(maxnslay+maxnzone) )
    double precision h5mx(-1:2,2*(maxnslay+maxnzone) )
    double precision h5my(-1:2,2*(maxnslay+maxnzone) )
    double precision h5mz(-1:2,2*(maxnslay+maxnzone) )
    double precision h4m1L(-1:2,2*(maxnslay+maxnzone) )
    double precision h4m1N(-1:2,2*(maxnslay+maxnzone) )
    double precision h4m2L(-2:1,2*(maxnslay+maxnzone) )
    double precision h4m2N(-2:1,2*(maxnslay+maxnzone) )
    double precision h6m1L(-1:2,2*(maxnslay+maxnzone) )
    double precision h6m1N(-1:2,2*(maxnslay+maxnzone) )
    double precision h6m2L(-2:1,2*(maxnslay+maxnzone) )
    double precision h6m2N(-2:1,2*(maxnslay+maxnzone) )
    double precision p1(8*maxnllay ),p2( 8*maxnllay ),p3( 8*maxnllay )
    complex(dp) g(2*(maxnslay+1) + (maxnllay+1) )
    complex(dp) d((maxnslay+1) + (maxnllay+1) )
    !c variable for the file
    character*80 output(maxnr)
    !c variable for the stack point
    integer isp(maxnzone)
    integer issp(maxnzone)
    integer ilsp(maxnzone),jssp(maxnzone)
    integer jsp(maxnzone)
    integer ksp(maxnzone),lsp(maxnzone)
    integer isdr,jsdr,ildr
    !c variables for the gridding
    integer jjdr(maxnzone),kkdr(maxnzone)
    integer jdr,kdr
    double precision vmin(maxnzone),gridpar(maxnzone),dzpar(maxnzone)
    !c variables for l cut off
    complex(dp) tmpc(2*(maxnslay+1) + (maxnllay+1))
    integer sufzone,ismall,kc,lsuf,llog
    double precision maxamp,ratc,ratl,re
    !c variables for the numerical integration
    complex(dp) anum(4,4,10),bnum(4,4,10)
    !c other variables
    integer i,j,ii,nn,lda,ir,ier,itmp,jtmp,mtmp,kkdr0,nn0
    integer ll(12),lli(12),llj(12)
    double precision eps,work(8*maxnslay),l2,lsq
    complex(dp) z( 2*(maxnslay+1) + (maxnllay+1) )
    complex(dp) w( 2*(maxnslay+1) + (maxnllay+1) )
    complex(dp) cwork( 2*(16*maxnslay + 4*maxnllay) )
    integer ltmp(2),iimax

    data lda/ 4 /
    data eps/ -1.d0 /
    complex(dp), intent(out) :: outputu(3,nr,imin:imax)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c *************** Inputting and computing the parameters ***************
!c --- inputting parameter ---
    !c --- computing the required parameters ---
    !c counting of the nsl and nll
    call calnl( nzone,vsv,iphase,nsl,nll )
    !c computing and checking the parameters
    rmin = vrmin(1)
    rmax = vrmax(nzone)
    ndc = nzone - 1

    theta=theta/180*pi
    phi=phi/180*pi

    if ( (r0<rmin) .or. (r0>rmax) ) stop 'Location of the source is improper.'

    iimax = imax

    if( (rmax-r0)<shallowdepth) then ! option for shallow events
        !c computing of the number and the location of grid points
        iimax = int(tlen * 2.d0)
        call calgrid( nzone,vrmin,vrmax,vpv,vsv,rmin,rmax,iimax,1,tlen, vmin,gridpar,dzpar )
        call calra(inlayer,jnlayer,jnslay,jnllay,&
                dzpar,nzone,vrmin,vrmax,iphase,rmin,nslay,nllay,nlayer,ra,re )
        !c --- checking the parameter
        if ( inlayer > maxnlay ) stop 'The number of grid points is too large.'
        if ( nslay > maxnslay ) stop 'The number of the grid points in the solid is too large.'
        if ( nllay > maxnllay ) stop 'The number of the grid points in the liquid is too large.'
        if ( jnlayer > 2*maxnlay ) stop 'The number of the total grid points is too large.'
        if ( ( jnslay > 2*maxnslay ).or.( jnllay>2*maxnllay ) ) stop 'The number of the total grid points is too large.'
        !c computing the stack points
        call calsp( maxnzone,ndc,nsl,nll,iphase,nlayer,nllay,&
                isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr,ildr,jdr,kdr )
        !c computing the source location
        call calspo( maxnlay,maxnzone,vrmax,iphase,inlayer,r0,rmin,rmax,ra,isp,spo,spn )
        !c ******************* Computing the matrix elements *******************

    !c computing the structure grid points
    call calstg( maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,&
                vnp,vra,rho,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN,r0,spn,ecC0,ecF0,ecL0 )

    rhoinv(1:vnp)=1.d0/rho(1:vnp)
    kappainv(1:vnp)=1.d0/kappa(1:vnp)
    isl = 0
    ill = 0
    do i=1,ndc+1
        if ( iphase(i)==1 ) then
            isl = isl + 1
            itmp = isdr+issp(isl)
            call calmatc( nlayer(i),vnp,vra,rho,2,0,0,ra(isp(i)), t(itmp) )
            call caltl( nlayer(i),vnp,vra,rho,ra(isp(i)),work(itmp) )
            !call calt( nlayer(i), t(itmp),work(itmp), t(itmp) )
      t(itmp:itmp+4*nlayer(i)-1)=(t(itmp:itmp+4*nlayer(i)-1)+work(itmp:itmp+4*nlayer(i)-1))/2d0
            call calmatc( nlayer(i),vnp,vra,ecKx,0,0,0,ra(isp(i)),h1x(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKx,ra(isp(i)),work(itmp) )
      !      call calt( nlayer(i), h1x(itmp),work(itmp),h1x(itmp) )
      h1x(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h1x(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecKy,0,0,0,ra(isp(i)),h1y(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKy,ra(isp(i)),work(itmp) )
        !    call calt( nlayer(i), h1y(itmp),work(itmp),h1y(itmp) )
            h1y(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h1y(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecKz,0,0,0,ra(isp(i)),h1z(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKz,ra(isp(i)),work(itmp) )
         !   call calt( nlayer(i), h1z(itmp),work(itmp),h1z(itmp) )
            h1z(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h1z(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecL,0,0,0,ra(isp(i)),h2L(itmp) )
            call calhl( nlayer(i),vnp,vra,ecL,ra(isp(i)),work(itmp) )
          !  call calt( nlayer(i), h2L(itmp),work(itmp),h2L(itmp) )
            h2L(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h2L(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h2N(itmp) )
            call calhl( nlayer(i),vnp,vra,ecN,ra(isp(i)),work(itmp) )
           ! call calt( nlayer(i), h2N(itmp),work(itmp),h2N(itmp) )
            h2N(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h2N(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecKx,1,0,1,ra(isp(i)),h5ax(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKy,1,0,1,ra(isp(i)),h5ay(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKz,1,0,1,ra(isp(i)),h5az(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecL,1,0,1,ra(isp(i)),h6aL(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecN,1,0,1,ra(isp(i)),h6aN(itmp) )
            call mtrnp( nlayer(i),h5ax(itmp),h3ax(itmp) )
            call mtrnp( nlayer(i),h5ay(itmp),h3ay(itmp) )
            call mtrnp( nlayer(i),h5az(itmp),h3az(itmp) )
            call mtrnp( nlayer(i),h6aL(itmp),h4aL(itmp) )
            call mtrnp( nlayer(i),h6aN(itmp),h4aN(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKx,2,1,1,ra(isp(i)), h7x(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKy,2,1,1,ra(isp(i)), h7y(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKz,2,1,1,ra(isp(i)), h7z(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecL,2,1,1,ra(isp(i)), h8L(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecN,2,1,1,ra(isp(i)), h8N(itmp) )
        else
            ill = ill + 1
            itmp = ildr+ilsp(ill)
            call calmatc( nlayer(i),vnp,vra,rhoinv,2,1,1,ra(isp(i)),p1(itmp) )
            call calmatc( nlayer(i),vnp,vra,rhoinv,0,0,0,ra(isp(i)),p2(itmp) )
            call calhl( nlayer(i),vnp,vra,rhoinv,ra(isp(i)),work(itmp) )
            !call calt( nlayer(i),p2(itmp),work(itmp),p2(itmp) )
            p2(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+p2(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,kappainv,2,0,0,ra(isp(i)),p3(itmp) )
            call caltl( nlayer(i),vnp,vra,kappainv,ra(isp(i)),work(itmp) )
            p3(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+p3(itmp:itmp+4*nlayer(i)-1))/2
        !call calt( nlayer(i),p3(itmp),work(itmp),p3(itmp) )
        endif
    enddo
 !       c Computing the modified operator of the 1st derivative
    call caltstg( maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,vra,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN )
    isl = 0
    do i=1,ndc+1
        if ( iphase(i)==1 ) then
            isl = isl + 1
            itmp = isdr+issp(isl)
            jtmp = isp(i)+i-1
            call calh5( nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),work(itmp) )
            h5x(itmp:itmp+4*nlayer(i)-1)= h5ax(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
            call calh5( nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),work(itmp) )
            h5y(itmp:itmp+4*nlayer(i)-1)= h5ay(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
            call calh5( nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),work(itmp) )
            h5z(itmp:itmp+4*nlayer(i)-1)= h5az(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
            call calh5( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),work(itmp) )
            h6L(itmp:itmp+4*nlayer(i)-1)= h6aL(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
            call calh5( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),work(itmp) )
            h6N(itmp:itmp+4*nlayer(i)-1)= h6aN(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
            call mtrnp( nlayer(i),h5x(itmp),h3x(itmp) )
            call mtrnp( nlayer(i),h5y(itmp),h3y(itmp) )
            call mtrnp( nlayer(i),h5z(itmp),h3z(itmp) )
            call mtrnp( nlayer(i),h6L(itmp),h4L(itmp) )
            call mtrnp( nlayer(i),h6N(itmp),h4N(itmp) )
            itmp = jsdr+jssp(isl)
            call calhm1( nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),h5mx(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),h5my(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),h5mz(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),h6m1L(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),h6m1N(-1,itmp) )
            call calhm2( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),h6m2L(-2,itmp) )
            call calhm2( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),h6m2N(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5mx(-1,itmp),h3mx(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5my(-1,itmp),h3my(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5mz(-1,itmp),h3mz(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h6m1L(-1,itmp),h4m2L(-2,itmp))
            call mtrnp2( nlayer(i),1,2,h6m1N(-1,itmp),h4m2N(-2,itmp))
            call mtrnp2( nlayer(i),2,1,h6m2L(-2,itmp),h4m1L(-1,itmp))
            call mtrnp2( nlayer(i),2,1,h6m2N(-2,itmp),h4m1N(-1,itmp))
        endif
    enddo

    do ii=1,2              ! omega-loop
        if(ii==1) then
            if(imin==0) then
                i=1
            else
                i=imin
            endif
        endif
        if(ii==2) i=imax
        if ( i/=0 ) then
            omega = 2.d0 * pi * dble(i) / tlen
            call callsuf(omega,nzone,vrmax,vsv,lsuf)
            call calcoef( nzone,omega,qmu,qkappa,coef1,coef2,coef )
            mtmp = isp(spn) + int(spo)
            if ( spo==int(spo) ) mtmp = mtmp - 1
                call calabnum( omega,omegai,rmax,rrho(1,spn),vpv(1,spn),vph(1,spn),&
                              vsv(1,spn),vsh(1,spn),eta(1,spn),ra(mtmp),r0,coef1(spn),coef2(spn),anum(1,1,1),bnum(1,1,1) )
                !c computing the matrix elements independent of l
                isl = 0
                ill = 0
                do j=1,ndc+1
                    if ( iphase(j)==1 ) then
                        isl = isl + 1
                        itmp = isdr+issp(isl)
                        jtmp = jdr+jsp(j)
                        mtmp = kdr+ksp(j)
                        call cala0( nlayer(j),omega,omegai,t(itmp),h1x(itmp), &
                                            h2L(itmp), h2N(itmp),h3ay(itmp),&
                                            h4aL(itmp),h4aN(itmp),h5ay(itmp),&
                                            h6aL(itmp),h6aN(itmp),h7y(itmp),h7z(itmp),&
                                            h8L(itmp), h8N(itmp),coef1(j),coef2(j),cwork(jtmp) )
                        call overlapa( nlayer(j),cwork(jtmp),a0(1,mtmp) )
                        call cala1( nlayer(j),h1x(itmp),&
                                            h2L(itmp),h2N(itmp),h3y(itmp),&
                                            h4L(itmp),h4N(itmp),h5y(itmp),&
                                            h6L(itmp),h6N(itmp),coef1(j),coef2(j),cwork(jtmp) )
                        call overlapa( nlayer(j),cwork(jtmp),a1(1,mtmp) )
                        call cala2(nlayer(j),h1x(itmp),h2L(itmp),h2N(itmp),coef1(j),coef2(j),cwork(jtmp))
                        call overlapa( nlayer(j),cwork(jtmp),a2(1,mtmp) )
                        jtmp = jsdr+jssp(isl)
                        call calhml( nlayer(j),coef1(j),coef2(j),&
                             h3my(-2,jtmp),h5my(-1,jtmp),h4m1L(-1,jtmp),h4m2N(-2,jtmp),&
                             h6m1N(-1,jtmp),h6m2L(-2,jtmp),a1(1,mtmp) )
                    else
                        ill = ill + 1
                        itmp = ildr+ilsp(ill)
                        jtmp = jdr+jsp(j)
                        mtmp = kdr+ksp(j)
                        call calb0( nlayer(j),omega,omegai,p1(itmp),p3(itmp),coef(j),cwork(jtmp))
                        call overlapb( nlayer(j),cwork(jtmp),a0(1,mtmp) )
                        call calb2( nlayer(j),omega,omegai,p2(itmp),cwork(jtmp))
                        call overlapb( nlayer(j),cwork(jtmp),a2(1,mtmp) )
                endif
            enddo

            maxamp = -1.0
            ismall = 0
            ltmp(ii) = maxlmax
            kc = 1
            do l=0,maxlmax   ! l-loop
                ! initialize
                tmpc(1:maxnlay+1) = dcmplx(0.d0)
                if ( ismall>20 ) then
                    if(ltmp(ii)>l) ltmp(ii) = l
                    exit
                endif
                l2 = dble(l)*dble(l+1)
                lsq = dsqrt( l2 )
                !    c computing the coefficient matrix elements
                 !   c --- renewing  mdr
                if (  mod(l,50)==0  ) then
                    call calmdr( omega,l,nzone,vrmax,vmin,rmax,sufzone )
                    call calspdr( maxnzone,nzone,iphase,nlayer,&
                                jjdr,kkdr )
                    nn = kkdr(nzone) + 2 * nlayer(nzone) + 1
                endif
                !c computing the matrix elements
                call cala( ndc,iphase,nlayer,kkdr,kdr,ksp,&
                                    l2,lsq,nn,a0,a1,a2,a )
                !c computing the boundary condition elements
                call calbc( ndc,vrmax,iphase,kkdr,a )

                jtmp = kkdr(spn) + 2 * int(spo)
                mtmp = isp(spn) + int(spo)
                if ( spo==int(spo) ) then
                    jtmp = jtmp - 2
                    mtmp = mtmp - 1
                endif
                call calya( anum(1,1,1),bnum(1,1,1),l2,ra(mtmp),r0,ya,yb,yc,yd )
                do m=-2,2     ! m-loop
                    if ( iabs(m)<=iabs(l) ) then
                !        c computing the excitation vector
                        g(1:nn)=0
                        call calg( l,m,coef1(spn),coef2(spn),lsq,&
                                               ecC0,ecF0,ecL0,ya,yb,yc,yd,&
                                               ra(mtmp),r0,mt,g(jtmp) )
                        if ( l==0 ) then
                            !c rearranging the matrix elements for l=0
                            call rea2( nn,a,g,c,d,nzone,iphase,kkdr,spn,kkdr0,nn0 )
                            !c computing the expansion coefficient
                            itmp=1
                            if ( rmin==0.d0 ) itmp=2
                            ns = kkdr0 + ( nint(spo) - 1 )
                            call dcsymbdl( c(1,itmp),1,nn0-itmp+1,1,eps,z(itmp),w(itmp),ll,lli,llj,ier )
                            call dcsbdlv( c(1,itmp),d(itmp),1,nn0-itmp+1,ns-itmp+1,eps,z(itmp),ier )
                        else
                            !c computing the expansion coefficient
                            itmp=1
                            if ( rmin==0.d0 ) itmp=3
                            ns = kkdr(spn) + 2 * ( nint(spo) - 1 )
                            if( mod(l,100)==0) then
                                if ( ( m==-2 ).or.( m==-l ) )&
                                 call dcsymbdl0( a(1,itmp),3,nn-itmp+1,6,eps, z(itmp),w(itmp),ll,lli,llj,ier )

                                call dcsbdlv0( a(1,itmp),g(itmp),3,nn-itmp+1,eps,z(itmp),ier )
                                ! sum up c of the same l
                                tmpc(1:nn) = tmpc(1:nn) + g(1:nn)
                            else
                                if ( kc<3 ) kc = 3
                                itmp = kc
                                if ( ( m==-2 ).or.( m==-l ) )&
                                 call dcsymbdl( a(1,itmp),3,nn-itmp+1,6,eps,z(itmp),w(itmp),ll,lli,llj,ier )

                                call dcsbdlv( a(1,itmp),g(itmp),3,nn-itmp+1,ns-itmp+1,eps,z(itmp),ier )
                            endif
                            !c computing ampratio
                            call calamp( g(nn-1),l,lsuf,maxamp,ismall,ratl )

                            if( mod(l,100)==0) call calcutd(nzone,nlayer,tmpc,ratc,nn,iphase,ra,kkdr,kc)

!                                c computing the displacement
                        endif
                    endif
                enddo         ! m-loop
            enddo            ! l-loop
        endif
    enddo                  ! omega-loop
        iimax = dble(max(ltmp(1),ltmp(2))) * tlen / lmaxdivf
    endif                     ! option for shallow events

    call calgrid( nzone,vrmin,vrmax,vpv,vsv,rmin,rmax,iimax,1,tlen, vmin,gridpar,dzpar )
    call calra(inlayer,jnlayer,jnslay,jnllay,&
         dzpar,nzone,vrmin,vrmax,iphase,rmin,nslay,nllay,nlayer,ra,re )
    !c --- checking the parameter
    if ( inlayer>maxnlay ) stop 'The number of grid points is too large.'
    if ( nslay>maxnslay ) stop 'The number of the grid points in the solid is too large.'
    if ( nllay>maxnllay ) stop 'The number of the grid points in the liquid is too large.'
    if ( jnlayer>2*maxnlay ) stop 'The number of the total grid points is too large.'
    if ( ( jnslay>2*maxnslay ).or.( jnllay>2*maxnllay ) ) stop 'The number of the total grid points is too large.'
    !c computing the stack points
    call calsp( maxnzone,ndc,nsl,nll,iphase,nlayer,nllay,&
         isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr,ildr,jdr,kdr )
!    c computing the source location
    call calspo(maxnlay,maxnzone,vrmax,iphase,inlayer,r0,rmin,rmax,ra,isp,spo,spn )

    !c ******************* Computing the matrix elements *******************
   !    c computing the structure grid points
    call calstg(maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,&
         vnp,vra,rho,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN,r0,spn,ecC0,ecF0,ecL0 )
    !call calinv( vnp,rho,kappa,rhoinv,kappainv )
    rhoinv(1:vnp) = 1.d0 / rho(1:vnp)
    kappainv(1:vnp) = 1.d0 / kappa(1:vnp)
    isl = 0
    ill = 0
    do i=1,ndc+1
        if ( iphase(i)==1 ) then
            isl = isl + 1
            itmp = isdr+issp(isl)
            call calmatc( nlayer(i),vnp,vra,rho,2,0,0,ra(isp(i)), t(itmp) )
            call caltl( nlayer(i),vnp,vra,rho,ra(isp(i)),work(itmp) )
            t(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+t(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecKx,0,0,0,ra(isp(i)),h1x(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKx,ra(isp(i)),work(itmp) )
            h1x(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h1x(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecKy,0,0,0,ra(isp(i)),h1y(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKy,ra(isp(i)),work(itmp) )
            h1y(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h1y(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecKz,0,0,0,ra(isp(i)),h1z(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKz,ra(isp(i)),work(itmp) )
            h1z(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h1z(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecL,0,0,0,ra(isp(i)),h2L(itmp) )
            call calhl( nlayer(i),vnp,vra,ecL,ra(isp(i)),work(itmp) )
            h2L(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h2L(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h2N(itmp) )
            call calhl( nlayer(i),vnp,vra,ecN,ra(isp(i)),work(itmp) )
            h2N(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+h2N(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,ecKx,1,0,1,ra(isp(i)),h5ax(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKy,1,0,1,ra(isp(i)),h5ay(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKz,1,0,1,ra(isp(i)),h5az(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecL,1,0,1,ra(isp(i)),h6aL(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecN,1,0,1,ra(isp(i)),h6aN(itmp) )
            call mtrnp( nlayer(i),h5ax(itmp),h3ax(itmp) )
            call mtrnp( nlayer(i),h5ay(itmp),h3ay(itmp) )
            call mtrnp( nlayer(i),h5az(itmp),h3az(itmp) )
            call mtrnp( nlayer(i),h6aL(itmp),h4aL(itmp) )
            call mtrnp( nlayer(i),h6aN(itmp),h4aN(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKx,2,1,1,ra(isp(i)), h7x(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKy,2,1,1,ra(isp(i)), h7y(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKz,2,1,1,ra(isp(i)), h7z(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecL,2,1,1,ra(isp(i)), h8L(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecN,2,1,1,ra(isp(i)), h8N(itmp) )
        else
            ill = ill + 1
            itmp = ildr+ilsp(ill)
            call calmatc( nlayer(i),vnp,vra,rhoinv,2,1,1,ra(isp(i)),p1(itmp) )
            call calmatc( nlayer(i),vnp,vra,rhoinv,0,0,0,ra(isp(i)),p2(itmp) )
            call calhl( nlayer(i),vnp,vra,rhoinv,ra(isp(i)),work(itmp) )
            p2(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+p2(itmp:itmp+4*nlayer(i)-1))/2
            call calmatc( nlayer(i),vnp,vra,kappainv,2,0,0,ra(isp(i)),p3(itmp) )
            call caltl( nlayer(i),vnp,vra,kappainv,ra(isp(i)),work(itmp) )
            p3(itmp:itmp+4*nlayer(i)-1)=(work(itmp:itmp+4*nlayer(i)-1)+p3(itmp:itmp+4*nlayer(i)-1))/2
        endif
    enddo
    !c Computing the modified operator of the 1st derivative
    call caltstg( maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,vra,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN )
    isl = 0
    do i=1,ndc+1
        if ( iphase(i)==1 ) then
            isl = isl + 1
!           !call submat( nlayer(i),h6aN(itmp),work(itmp),h6N(itmp) )

            itmp = isdr+issp(isl)
            jtmp = isp(i)+i-1
            call calh5( nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),work(itmp) )
 !           call submat( nlayer(i),h5ax(itmp),work(itmp),h5x(itmp) )
       h5x(itmp:itmp+4*nlayer(i)-1)= h5ax(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
            call calh5( nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),work(itmp) )
!            call submat( nlayer(i),h5ay(itmp),work(itmp),h5y(itmp) )
       h5y(itmp:itmp+4*nlayer(i)-1)= h5ay(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
            call calh5( nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),work(itmp) )
!            call submat( nlayer(i),h5az(itmp),work(itmp),h5z(itmp) )
        h5z(itmp:itmp+4*nlayer(i)-1)= h5az(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
            call calh5( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),work(itmp) )
!            call submat( nlayer(i),h6aL(itmp),work(itmp),h6L(itmp) )
             h6L(itmp:itmp+4*nlayer(i)-1)= h6aL(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)

            call calh5( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),work(itmp) )
!            call submat( nlayer(i),h6aN(itmp),work(itmp),h6N(itmp) )
            h6N(itmp:itmp+4*nlayer(i)-1)= h6aN(itmp:itmp+4*nlayer(i)-1)-work(itmp:itmp+4*nlayer(i)-1)
           call mtrnp( nlayer(i),h5x(itmp),h3x(itmp) )
            call mtrnp( nlayer(i),h5y(itmp),h3y(itmp) )
            call mtrnp( nlayer(i),h5z(itmp),h3z(itmp) )
            call mtrnp( nlayer(i),h6L(itmp),h4L(itmp) )
            call mtrnp( nlayer(i),h6N(itmp),h4N(itmp) )
            itmp = jsdr+jssp(isl)
            call calhm1( nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),h5mx(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),h5my(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),h5mz(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),h6m1L(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),h6m1N(-1,itmp) )
            call calhm2( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),h6m2L(-2,itmp) )
            call calhm2( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),h6m2N(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5mx(-1,itmp),h3mx(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5my(-1,itmp),h3my(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5mz(-1,itmp),h3mz(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h6m1L(-1,itmp),h4m2L(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h6m1N(-1,itmp),h4m2N(-2,itmp) )
            call mtrnp2( nlayer(i),2,1,h6m2L(-2,itmp),h4m1L(-1,itmp) )
            call mtrnp2( nlayer(i),2,1,h6m2N(-2,itmp),h4m1N(-1,itmp) )
        endif
    enddo
!    c ******************** Computing the displacement *********************

    llog = 0
    do i=imin,imax            ! omega-loop
        u(1:3,1:nr)=0
        if ( i/=0 ) then
            omega = 2.d0 * pi * i / tlen
            call callsuf(omega,nzone,vrmax,vsv,lsuf)

            plm(1:3,0:3,1:nr)=0
            call calcoef( nzone,omega,qmu,qkappa,coef1,coef2,coef )
            mtmp = isp(spn) + int(spo)
            if ( spo==int(spo) ) mtmp = mtmp - 1
            call calabnum( omega,omegai,rmax,rrho(1,spn),vpv(1,spn),vph(1,spn),vsv(1,spn),vsh(1,spn),&
             eta(1,spn),ra(mtmp),r0,coef1(spn),coef2(spn),anum(1,1,1),bnum(1,1,1) )
!            c computing the matrix elements independent of l
            isl = 0
            ill = 0
            do j=1,ndc+1
                if ( iphase(j)==1 ) then
                    isl = isl + 1
                    itmp = isdr+issp(isl)
                    jtmp = jdr+jsp(j)
                    mtmp = kdr+ksp(j)
                    call cala0( nlayer(j),omega,omegai,t(itmp),h1x(itmp),&
                     h2L(itmp),h2N(itmp),h3ay(itmp),&
                     h4aL(itmp),h4aN(itmp),h5ay(itmp),&
                     h6aL(itmp),h6aN(itmp),h7y(itmp),h7z(itmp),&
                     h8L(itmp),h8N(itmp),coef1(j),coef2(j),cwork(jtmp) )
                    call overlapa( nlayer(j),cwork(jtmp),a0(1,mtmp))
                    call cala1( nlayer(j),h1x(itmp),&
                     h2L(itmp),h2N(itmp),h3y(itmp),&
                     h4L(itmp),h4N(itmp),h5y(itmp),&
                     h6L(itmp),h6N(itmp),coef1(j),coef2(j),cwork(jtmp) )
                    call overlapa( nlayer(j),cwork(jtmp),a1(1,mtmp) )
                    call cala2( nlayer(j),h1x(itmp),h2L(itmp),h2N(itmp),coef1(j),coef2(j),cwork(jtmp))
                    call overlapa( nlayer(j),cwork(jtmp),a2(1,mtmp) )
                    jtmp = jsdr+jssp(isl)
                    call calhml( nlayer(j),coef1(j),coef2(j),&
                     h3my(-2,jtmp),h5my(-1,jtmp),&
                     h4m1L(-1,jtmp),h4m2N(-2,jtmp),h6m1N(-1,jtmp),h6m2L(-2,jtmp),a1(1,mtmp) )
                else
                    ill = ill + 1
                    itmp = ildr+ilsp(ill)
                    jtmp = jdr+jsp(j)
                    mtmp = kdr+ksp(j)
                    call calb0( nlayer(j),omega,omegai,p1(itmp),p3(itmp),coef(j),cwork(jtmp) )
                    call overlapb( nlayer(j),cwork(jtmp),a0(1,mtmp) )
                    call calb2( nlayer(j),omega,omegai,p2(itmp),cwork(jtmp))
                    call overlapb( nlayer(j),cwork(jtmp),a2(1,mtmp) )
                endif
            enddo
            maxamp = -1
            ismall = 0
            llog = maxlmax
            kc = 1
            do l=0,maxlmax      ! l-loop
                tmpc(1:maxnlay+1) = 0  ! initialize
                if ( ismall>20 ) then
                    if(llog>l) llog = l
                    exit
                endif
                l2 = dble(l*(l+1))
                lsq = dsqrt( l2 )
                !c computing the coefficient matrix elements
                !c --- renewing  mdr
                if (  mod(l,50)==0  ) then
                    call calmdr( omega,l,nzone,vrmax,vmin,rmax,sufzone )
                    call calspdr( maxnzone,nzone,iphase,nlayer,jjdr,kkdr )
                    nn = kkdr(nzone) + 2 * nlayer(nzone) + 1
                endif
                !c ***** Computing the trial function *****
                do ir=1,nr
                    call calbvec( l,theta(ir),phi(ir),plm(1,0,ir),bvec(1,-2,ir) )
                enddo
!              c computing the matrix elements
                call cala( ndc,iphase,nlayer,kkdr,kdr,ksp,l2,lsq,nn,a0,a1,a2,a )
!                c computing the boundary condition elements
                call calbc( ndc,vrmax,iphase,kkdr,a )

                jtmp = kkdr(spn) + 2 * int(spo)
                mtmp = isp(spn) + int(spo)
                if ( spo==int(spo) ) then
                    jtmp = jtmp - 2
                    mtmp = mtmp - 1
                endif
                call calya( anum(1,1,1),bnum(1,1,1),l2,ra(mtmp),r0,ya,yb,yc,yd )
                do m=-2,2        ! m-loop
                    if ( iabs(m)<=iabs(l) ) then
                        !c computing the excitation vector
                        g(1:nn)=0
                        call calg( l,m,coef1(spn),coef2(spn),lsq,ecC0,ecF0,ecL0,ya,yb,yc,yd,ra(mtmp),r0,mt,g(jtmp) )
                        if ( l==0 ) then
!                            c rearranging the matrix elements for l=0
                            call rea2( nn,a,g,c,d,nzone,iphase,kkdr,spn,kkdr0,nn0 )
                            !c computing the expansion coefficient
                            itmp=1
                            if ( rmin==0.d0 ) itmp=2
                            ns = kkdr0 + ( nint(spo) - 1 )
                            call dcsymbdl( c(1,itmp),1,nn0-itmp+1,1,eps,z(itmp),w(itmp),ll,lli,llj,ier )
                            call dcsbdlv( c(1,itmp),d(itmp),1,nn0-itmp+1,ns-itmp+1,eps,z(itmp),ier )
                            !c computing the displacement
                            u(:,1:nr)=d(nn0)*bvec(:,m,1:nr)
                        else
                            !c computing the expansion coefficient
                            itmp=1
                            if ( rmin==0.d0 ) itmp=3
                            ns = kkdr(spn) + 2 * ( nint(spo) - 1 )
                            if( mod(l,100)==0) then
                                if ( ( m==-2 ).or.( m==-l ) ) then
                                    call dcsymbdl0( a(1,itmp),3,nn-itmp+1,6,eps,z(itmp),w(itmp),ll,lli,llj,ier )
                                endif
                                call dcsbdlv0( a(1,itmp),g(itmp),3,nn-itmp+1,eps,z(itmp),ier )
                                 ! sum up c of the same l
                                tmpc(1:nn) = tmpc(1:nn) + g(1:nn)
                            else
                                if ( kc<3 ) kc = 3
                                itmp = kc
                                if ( ( m==-2 ).or.( m==-l ) )&
                                 call dcsymbdl( a(1,itmp),3,nn-itmp+1,6,eps,z(itmp),w(itmp),ll,lli,llj,ier )
                                call dcsbdlv( a(1,itmp),g(itmp),3,nn-itmp+1,ns-itmp+1,eps,z(itmp),ier )
                            endif
                            !c computing ampratio
                            call calamp( g(nn-1),l,lsuf,maxamp,ismall,ratl )
                            if( mod(l,100)==0)&
                             call calcutd(nzone,nlayer,tmpc,ratc,nn,iphase,ra,kkdr,kc)

                            !c computing the displacement calu
                            u(1,1:nr) = u(1,1:nr) + g(nn-1) * bvec(1,m,1:nr)
                            u(2:3,1:nr) = u(2:3,1:nr) + g(nn) * bvec(2:3,m,1:nr) / lsq
                        endif
                    endif
                enddo            ! m-loop
            enddo               ! l-loop
        endif
        outputu (1:3,1:nr,i)=u(1:3,1:nr)
    enddo                   ! omega-loop

!c ************************** Files Handling **************************
    if (write_to_file) then
        write(*,*) "kakikomimasu"
        do ir=1,nr
            open(unit=10,file=trim(output(ir)),status='unknown',form='unformatted',access='stream',convert='big_endian')
            write(10) tlen,np,1,3,omegai,lat(ir),lon(ir),eqlat,eqlon,r0

            do i= imin, imax
                write(10) i, dble(outputu(1,ir,i)), dimag(outputu(1,ir,i))
                write(10) dble(outputu(2,ir,i)), dimag(outputu(2,ir,i))
                write(10) dble(outputu(3,ir,i)), dimag(outputu(3,ir,i))
            enddo
            close(10)
        enddo
        write(*,*) "Ivalice looks to the horizon"
    endif

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    return
 end subroutine tipsv


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
   real(dp) :: qmu(maxnzone),qkappa(maxnzone)
   real(dp) :: rho(4,maxnzone)
   real(dp) :: vpv(4,maxnzone),vph(4,maxnzone),vsv(4,maxnzone),vsh(4,maxnzone),eta(4,maxnzone)
   real(dp) :: r0
   real(dp) :: mt(3,3)
   real(dp) :: lat(maxnr)
   real(dp) :: lon(maxnr)
   real(dp) :: theta(maxnr)
   real(dp) :: phi(maxnr)
   real(dp) :: eqlat
   real(dp) :: eqlon
   character(len=80) :: output(maxnr)
   complex(dp), allocatable :: outputu(:,:,:)
   logical :: write_to_file = .true.

   ! read input parameters
   call get_command_argument(1, parameter_file)
   print*,parameter_file

   call pinput_fromfile(parameter_file, &
      re,ratc,ratl,tlen,np,omegai,imin,imax, &
      nzone,vrmin,vrmax,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa, &
      r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)

   allocate(outputu(3,nr,imin:imax))

   ! main loop
   write(*,*) 'Enter main loop'

   call tipsv(re,ratc,ratl,tlen,np,omegai,imin,imax, &
      nzone,vrmin,vrmax,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,&
      r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output,write_to_file,&
      outputu)
   write(*,*) 'Done!'

end program main


















