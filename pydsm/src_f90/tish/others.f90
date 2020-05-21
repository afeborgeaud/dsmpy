subroutine pinput_tish(parameter_file, &
    re,ratc,ratl,tlen,np,omegai,imin,imax, &
    nzone,vrmin,vrmax,rho,vsv,vsh,qmu, &
    r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)
!----------------------------------------------------------
! Parameter Input
!----------------------------------------------------------
    use parameters
    implicit none
    integer, parameter :: maxnline=2*maxnr+6*maxnzone+50
    character*160, intent(in) :: parameter_file
    integer, intent(out):: np,imin,imax,nzone,nr
    real(dp),intent(out) :: tlen,omegai,re,ratc,ratl
    real(dp),dimension(maxnzone), intent(out):: vrmin,vrmax,qmu
    real(dp),dimension(4,maxnzone), intent(out):: rho,vsv,vsh
    real(dp),dimension(maxnr), intent(out) :: theta,phi,lat,lon
    real(dp),intent(out) :: eqlat,eqlon,r0,mt(3,3)
    character*80,dimension(maxnr),intent(out) :: output
    real(dp) :: stlat,stlon,eqlattmp
    integer i,linenum,io
    logical:: file_exists
    character*80::buffer
    character*80,dimension(maxnline) :: lines

    inquire(file=parameter_file,exist=file_exists)
    if (.not. file_exists) stop 'parameter file does not exist.'

    linenum=0
    open(unit=1,file=parameter_file,status='old',action='read')
    do
        read(1, '(a)', iostat=io) buffer
        buffer = adjustl(buffer)
        if(buffer(1:1)=='c'.or.buffer(1:1)=='c'.or.buffer(1:1)=='!') cycle
        if(io/=0) exit
        linenum=linenum+1
        lines(linenum) = buffer
    enddo
    close(1)

    read(lines(1),*) tlen,np
    read(lines(2),*) re ! relative error (vertical grid)
    read(lines(3),*) ratc ! ampratio (vertical grid cut-off)
    read(lines(4),*) ratl ! ampratio (for l-cutoff)
    read(lines(5),*) omegai ! omegai
    omegai = -dlog(omegai)/tlen
    read(lines(6),*) imin,imax
    read(lines(7),*) nzone
    if (nzone > maxnzone) stop 'nzone is too large. (pinput)'

! structure
    do i=1,nzone
        read(lines(7+3*(i-1)+1),*) vrmin(i),vrmax(i),rho(1:4,i)
        read(lines(7+3*(i-1)+2),*) vsv(1:4,i)
        read(lines(7+3*(i-1)+3),*) vsh(1:4,i),qmu(i)
    enddo
! source parameter
    read(lines(3*nzone+8),*) r0,eqlat,eqlon
    eqlattmp=eqlat
    call translat(eqlattmp,eqlattmp)
    read(lines(3*nzone+9),*) mt(1,1:3), mt(2,2:3), mt(3,3)
    read(lines(3*nzone+10),*) nr
! station
    if (nr > maxnr) stop 'nr is too large. (pinput)'
    do i=1,nr
        read(lines(3*nzone+10+i),*) lat(i),lon(i)
        stlat = lat(i)
        stlon = lon(i)
        call translat(stlat,stlat)
        call calthetaphi(eqlattmp,eqlon,stlat,stlon,theta(i),phi(i))
    enddo
    theta(1:nr) = theta(1:nr) / 1.8d2 * pi
    phi(1:nr) = phi(1:nr) / 1.8d2 * pi
    do i=1,nr
        read(lines(3*nzone+10+nr+i),'(a)') output(i)
        output(i)=trim(output(i))
    enddo
    return
    end

!----------------------------------------------------------
subroutine calthetaphi(ievla,ievlo,istla,istlo,theta,phi)
!----------------------------------------------------------
    use parameters
    implicit none
    double precision:: ievla,ievlo,istla,istlo
    double precision:: evla,evlo,stla,stlo
    double precision:: theta,phi
    double precision:: gcarc,az
    double precision:: tc,ts

    ! transformation to spherical coordinates

    evla = 90.d0 - ievla
    stla = 90.d0 - istla

    evla = evla / 1.8d2 * pi
    evlo = ievlo / 1.8d2 * pi
    stla = stla / 1.8d2 * pi
    stlo = istlo / 1.8d2 * pi

    gcarc = dacos( dcos(evla) * dcos(stla)&
        + dsin(evla) * dsin(stla) * dcos(evlo - stlo) )

    tc = ( dcos(stla) * dsin(evla) &
        - dsin(stla) * dcos(evla) * dcos(stlo - evlo) )&
        / dsin(gcarc)
    ts = dsin(stla) * dsin(stlo - evlo) / dsin(gcarc)

    az = dacos(tc)
    if( ts < 0.d0 ) az = -1.d0 * az

    az = az * 1.8d2 / pi

    gcarc = gcarc * 1.8d2 / pi

    theta = gcarc
    phi   = 180.d0 - az
    return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine translat(geodetic,geocentric)
!----------------------------------------------------------
    use parameters
    implicit none
    double precision:: geocentric, geodetic

    integer:: flag

    flag = 0
    if(geodetic > 90.d0) then
        geodetic = 180 - geodetic
        flag = 1
    endif

    geodetic = geodetic / 180 * pi
    geocentric = datan( (1 - flattening) * (1 - flattening)&
        * dtan(geodetic) )
    geocentric = geocentric * 180 / pi
    if(flag == 1) geocentric = 180 - geocentric

    return
end

!----------------------------------------------------------
subroutine calgrid( nzone,vrmin,vrmax,vs,rmin,rmax,&
    imax,lmin,tlen,vmin,gridpar,dzpar )
!----------------------------------------------------------
    use parameters
    implicit none
    integer:: nzone,imax,lmin
    double precision:: vrmin(*),vrmax(*),vs(4,*)
    double precision:: rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
    integer:: izone,j
    double precision:: coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp

    do izone=1,nzone
        ! computing the S-velocity at each zone
        v(:) = vs(:,izone)
        vs1 = 0.d0
        vs2 = 0.d0
        do j=1,4
            if ( j==1 ) then
                coef1 = 1.d0
            else
                coef1 = coef1 * ( vrmin(izone) / rmax )
            endif
            if ( j==1 ) then
                coef2 = 1.d0
            else
                coef2 = coef2 * ( vrmax(izone) / rmax )
            endif
            vs1 = vs1 + v(j) * coef1
            vs2 = vs2 + v(j) * coef2
        enddo
        ! computing rh
        rh = vrmax(izone) - vrmin(izone)
        ! computing omega,amax
        omega = 2.d0 * pi * dble(imax) / tlen
        if ( vs1>=vs2 ) then
            vmin(izone) = vs2
        else
            vmin(izone) = vs1
        endif
        amax = vrmax(izone)
        gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) ) &
            - ( (dble(lmin)+0.5d0) * (dble(lmin)+0.5d0) )&
            / ( amax * amax )
        if ( gtmp>0.d0 ) then
            dzpar(izone)   = dsqrt( 1.d0/gtmp )
            gridpar(izone) = rh / dzpar(izone)
        else
            dzpar(izone)   = 0.d0
            gridpar(izone) = 0.d0
        endif
    enddo
    ! rearangement of gridpar
    gtmp = sum(gridpar(1:nzone))
    do izone=1,nzone
        if ( gridpar(izone)>0.d0 ) then
            gridpar(izone) = gridpar(izone) / gtmp
        else
            rh = vrmax(izone) - vrmin(izone)
            gridpar(izone) = rh / ( rmax - rmin ) * 0.1d0
        endif
    enddo
    ! re-rearangement of gridpar
    gtmp = sum(gridpar(1:nzone))
    gridpar(1:nzone)=gridpar(1:nzone)/gtmp

    return
end

!----------------------------------------------------------
subroutine calra(nlayer,gridpar,dzpar,nzone,vrmin,vrmax,&
    rmin,rmax,nnl,ra,re )
!----------------------------------------------------------
!c Computing the number and the location of grid points.
!----------------------------------------------------------
    use parameters
    implicit none
    integer:: nlayer
    integer:: nzone,nnl(maxnzone)
    double precision:: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax
    double precision:: ra(maxnlay+maxnzone+1)
    integer:: izone,itmp,i,ntmp
    double precision:: rh,re

    ! Initializing the data
    ra(1:maxnlay+maxnzone+1) = 0.d0
    nnl(1:nzone) = 0

    ! computing the number and the location of the grid points
    ra(1) = rmin
    itmp = 1
    do  izone=1,nzone
        rh = vrmax(izone) - vrmin(izone)
        if(dzpar(izone)==0.d0) then
            ntmp = 1
        else
            ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) &
                / 2.d0 / pi  / 7.d-1 + 1 )
        endif
        ! ntmp (see Geller & Takeuchi 1995 6.2)
        nnl(izone) = ntmp
        if ( nnl(izone)<5 ) nnl(izone)=5
        do i=1,nnl(izone)
            itmp = itmp + 1
            ra(itmp) = vrmin(izone)&
                + rh * dble(i) / dble( nnl(izone) )
        enddo
    enddo

    ! recouting the total number of grid points
    nlayer =sum(nnl(1:nzone))

    return
end

!----------------------------------------------------------
subroutine calra0( nlayer,nzone,vrmin,vrmax,nnl,ra )
!----------------------------------------------------------
!c Computing the number and the location of grid points.
!----------------------------------------------------------
    integer:: nlayer,nzone,nnl(nzone)
    double precision:: vrmin(*),vrmax(*),ra(*)
    integer:: izone,itmp,i
    double precision:: rmin,rmax,rh

    ! computing the number and the location of the grid points
    rmin = vrmin(1)
    rmax = vrmax(nzone)
    ra(1) = rmin
    itmp = 1
    do izone=1,nzone
        rh = vrmax(izone) - vrmin(izone)
        nnl(izone) = int( dble(nlayer) * rh / ( rmax - rmin ) ) + 1
        do i=1,nnl(izone)
            itmp = itmp + 1
            ra(itmp) = vrmin(izone) + rh * dble(i) / dble( nnl(izone) )
        enddo
    enddo
    ! recouting the total number of grid points
    nlayer = sum(nnl)

    return
end

!----------------------------------------------------------
subroutine calgra( isp,ra,r0,spn,spo,gra )
!----------------------------------------------------------
    integer:: isp(*),spn,itmp
    double precision,intent(in):: ra(*),r0,spo
    double precision,intent(out):: gra(3)

    itmp = isp(spn) + dint( spo )
    gra(1) = ra(itmp)
    gra(2) = r0
    gra(3) = ra(itmp+1)
    return
end

!----------------------------------------------------------
subroutine calsp( ndc,nlayer,isp,jsp )
!----------------------------------------------------------
! Computing the stack points.
!----------------------------------------------------------
    integer:: ndc,nlayer(*)
    integer:: isp(*),jsp(*)
    integer:: i

    ! computation of isp,jsp,ksp,lsp
    isp(1) = 1
    jsp(1) = 1

    do i=1,ndc
        isp(i+1) = isp(i) + nlayer(i)
        jsp(i+1) = jsp(i) + 4 * nlayer(i)
    enddo

    return
end

!----------------------------------------------------------
subroutine calspo( ndc,rdc,nlayer,r0,rmin,rmax,ra,&
    isp,spo,spn )
!----------------------------------------------------------
! Computing the source location.

    integer:: ndc,nlayer,isp(*),spn
    double precision:: rdc(*),r0,rmin,rmax,ra(*),spo
    integer:: itmp

    ! checking the parameter
    if ( r0<rmin .or. r0>rmax ) stop 'The source location is improper.(calspo)'
    ! computing 'spo'
    if ( r0==rmax ) then
        spo = dble(nlayer) - 0.01d0
        r0 = ra(nlayer) + (spo-dble(nlayer-1)) * ( ra(nlayer+1)-ra(nlayer) )
    else
        itmp = 2
        do
            if(r0<ra(itmp)) exit
            itmp=itmp+1
        enddo

        spo = dble(itmp-2) + ( r0-ra(itmp-1) ) / ( ra(itmp)-ra(itmp-1) )
        !c temporal handling
        if ( (spo-dble(itmp-2))<0.01d0 ) then
            spo = dble(itmp-2) + 0.01d0
            r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
        endif
        if ( (spo-dble(itmp-2))>0.99d0 ) then
            spo = dble(itmp-2) + 0.99d0
            r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
        endif

    endif
! computing 'spn'
    spn = 0
    itmp = 1

    do
        spn=spn+1
        if(r0<=rdc(itmp))exit
        itmp=itmp+1
    enddo

    ! changing 'spo'
    spo = spo - dble( isp(spn) - 1 )

    return
end


!----------------------------------------------------------
subroutine calstg( nzone,rrho,vsv,vsh,&
    nlayer,nnl,ra,rmax,vnp,vra,rho,&
    ecL,ecN )
!----------------------------------------------------------
!c Computing the structure grid points.
!----------------------------------------------------------
    implicit none
    integer:: nzone,nlayer,nnl(*),vnp
    double precision:: rrho(4,*),vsv(4,*),vsh(4,*),ra(*),rmax
    double precision:: vra(*),rho(*),ecL(*),ecN(*)
    double precision:: trho,tvsv,tvsh,coef
    integer:: izone,i,j,itmp,jtmp

    ! initializing the data
    vra(1:nlayer+nzone+1) = 0.d0
    rho(1:nlayer+nzone+1) = 0.d0
    ecL(1:nlayer+nzone+1) = 0.d0
    ecN(1:nlayer+nzone+1) = 0.d0

    !c computing the structure grid points
    itmp = 0
    jtmp = 0
    do izone=1,nzone
        do i=1,nnl(izone)+1
            itmp = itmp + 1
            jtmp = jtmp + 1
            vra(itmp) = ra(jtmp)
    !c --- evaluating the density and elastic constants at this point
            trho = 0.d0
            tvsv = 0.d0
            tvsh = 0.d0
            do j=1,4
                if ( j==1 ) then
                    coef = 1.d0
                else
                    coef = coef * ( vra(itmp) / rmax )
                endif
                trho = trho + rrho(j,izone) * coef
                tvsv = tvsv  + vsv(j,izone)   * coef
                tvsh = tvsh  + vsh(j,izone)   * coef
            enddo
            rho(itmp) = trho
            ecL(itmp) = rho(itmp) * tvsv * tvsv
            ecN(itmp) = rho(itmp) * tvsh * tvsh
        enddo
        jtmp = jtmp - 1
    enddo
    vnp = itmp
    return
end

!----------------------------------------------------------
subroutine calgstg( spn,rrho,vsv,vsh,&
    ra,vra,rmax,rho,ecL,ecN,r0,mu0 )
!----------------------------------------------------------
!c Computing the structure grid points.
    implicit none
    integer:: spn
    double precision:: rrho(4,*),vsv(4,*),vsh(4,*)
    double precision:: ra(*),rmax
    double precision:: vra(*),rho(*),ecL(*),ecN(*),r0,mu0
    double precision:: trho,tvsv,tvsh,coef
    integer:: i,j

    ! initializing the data
    vra(1:3) = 0.d0
    rho(1:3) = 0.d0
    ecL(1:3) = 0.d0
    ecN(1:3) = 0.d0

    !c computing the structure grid points
    do i=1,3
        vra(i) = ra(i)
        !c --- evaluating the density and elastic constants at this point
        trho = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        do j=1,4
            if ( j==1 ) then
                coef = 1.d0
            else
                coef = coef * ( vra(i) / rmax )
            endif
            trho = trho + rrho(j,spn) * coef
            tvsv = tvsv + vsv(j,spn) * coef
            tvsh = tvsh + vsh(j,spn) * coef
        enddo
        rho(i) = trho
        ecL(i) = rho(i) * tvsv * tvsv
        ecN(i) = rho(i) * tvsh * tvsh
    enddo

    mu0 = ecL(2)

    return
end

!----------------------------------------------------------
subroutine calcoef( nzone,omega,q,coef )
!----------------------------------------------------------
    use parameters
    implicit none

    integer:: izone
    integer,intent(in):: nzone
    double precision,intent(in):: omega,q(*)
    complex(dp),intent(out):: coef(*)
    double precision:: aa,bb

    do izone=1,nzone
        if ( omega==0.d0 ) then
            aa = 1.d0
        else
            aa = 1.d0 + dlog( omega / ( 2.d0 * pi ) ) / ( pi * Q(izone) )
        endif
        bb = 1.d0 / ( 2.d0 * Q(izone) )
        coef(izone) = dcmplx( aa, bb ) * dcmplx( aa, bb )
    enddo

    return
end

!----------------------------------------------------------
subroutine calcutd(nzone,nnl,tmpr,rat,nn,ra,kc)
!----------------------------------------------------------
    use parameters
    implicit none
    integer,intent(in):: nzone,nn,nnl(*)
    integer,intent(out)::kc
    complex(dp),intent(in):: tmpr(*)
    double precision:: rat,ra(*)
    integer:: nc

    double precision:: cU(nn),rc
    double precision:: maxamp,amp(nn)
    integer:: jz,jj,i,ml(nzone),tzone

    cU(1:nn) = tmpr(1:nn)

    maxamp = -1.d0
    do i=1,nn
        amp(i) = cU(i)
        maxamp = max(amp(i),maxamp)
    enddo

    maxamp = maxamp * rat ! threshold value
    if(maxamp==0.d0) then
        kc = 1
        return
    endif

    do i=1,nn
        if(amp(i)>maxamp) then
            nc = i
            exit
        endif
    enddo

    i = 1
    do jj=1,nzone
        i = i + nnl(jj)
        ml(jj) = i
    enddo

    do jj=nzone,1,-1
        if(ml(jj)>nc) tzone = jj
    enddo

    rc = ra(nc)
    kc = nc
    return
end

!----------------------------------------------------------
subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)
!----------------------------------------------------------
    implicit none
    integer,intent(in):: nzone
    integer,intent(out)::lsuf
    double precision,intent(in):: omega,vrmax(*),vsv(4,*)
    double precision:: tvs
    integer:: i

    tvs = 0.d0

    do i=1,4
        tvs = tvs + ( vsv(i,nzone) )
    enddo

    lsuf = int(omega * vrmax(nzone) / tvs - 0.5d0) + 1
    return
end

!------------------------------------------------------------------------
subroutine calamp(g,l,lsuf,maxamp,ismall,ratl)
!------------------------------------------------------------------------
    use parameters
    implicit none
    integer:: l,lsuf,ismall
    double precision:: maxamp,ratl
    complex(dp),intent(in):: g

    double precision:: amp,ampratio

    ampratio = 0.d0
    amp = cdabs(g)
    if( maxamp<amp ) maxamp = amp
    if ( amp /= 0.d0 .and. maxamp /= 0.d0 ) then
        ampratio = amp / maxamp
    endif
    if( ampratio<ratl .and. lsuf<l ) then
        ismall = ismall + 1
    else
        ismall = 0
    endif
    return
end
