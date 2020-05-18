!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine pinput_fromfile(parameter_file,&
    re,ratc,ratl,tlen,np,omegai,imin,imax,&
    nzone,vrmin,vrmax,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,&
    r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)
    use parameters
    implicit none
    character*160, intent(in):: parameter_file
    integer,intent(out):: np,imin,imax,nzone,nr
    double precision,intent(out):: tlen,omegai,re,ratc,ratl
    double precision,dimension(maxnzone),intent(out):: vrmin,vrmax,qmu,qkappa
    double precision,dimension(4,maxnzone),intent(out)::rho,vpv,vph,vsv,vsh,eta
    double precision,dimension(maxnr),intent(out)::theta,phi,lat,lon
    double precision,intent(out):: eqlat,eqlon,r0,mt(3,3)
    character*80,dimension(maxnr),intent(out):: output
    integer:: i,linenum,io
    logical:: file_exists
    double precision:: eqlattmp,stlat,stlon
    character*80::buffer
    character*80,dimension(1000):: lines

    inquire(file=parameter_file, exist=file_exists)

    if (.not. file_exists) stop 'parameter file does not exist.'

    linenum=0
    open(1, file=parameter_file,status='old',action='read')
    do
        read (1,'(a)', iostat=io) buffer
        buffer = adjustl(buffer)
        if(buffer(1:1) =='#'.or.buffer(1:1) =='c' .or. buffer(1:1) == '!') cycle
        if (io/=0) exit
        linenum = linenum+1
        lines(linenum) = buffer
    enddo
    close(1)

! reading the parameter
    read(lines(1),*) tlen,np
    read(lines(2),*) re! relative error (vertical grid)
    read(lines(3),*) ratc! ampratio (vertical grid cut-off)
    read(lines(4),*) ratl! ampratio (for l-cutoff)
    read(lines(5),*) omegai! artificial damping
    omegai = - dlog(omegai) / tlen
    read(lines(6),*) imin,imax
    read(lines(7),*) nzone
    if ( nzone>maxnzone ) stop 'nzone is too large. (pinput)'
    do i=1, nzone
        read(lines(7+6*(i-1)+1),*) vrmin(i),vrmax(i),rho(1:4,i)
        read(lines(7+6*(i-1)+2),*) vpv(1:4,i)
        read(lines(7+6*(i-1)+3),*) vph(1:4,i)
        read(lines(7+6*(i-1)+4),*) vsv(1:4,i)
        read(lines(7+6*(i-1)+5),*) vsh(1:4,i)
        read(lines(7+6*(i-1)+6),*) eta(1:4,i),qmu(i),qkappa(i)
    enddo
    read(lines(8+6*nzone),*) r0,eqlat,eqlon
    eqlattmp = eqlat
    call translat(eqlattmp,eqlattmp)
    read(lines(9+6*nzone),*) mt(1,1:3),mt(2,2:3),mt(3,3)
    read(lines(10+6*nzone),*) nr
    if ( nr>maxnr ) stop 'nr is too large. (pinput)'
    do  i=1,nr
        read(lines(10+6*nzone+i),*) lat(i),lon(i)
        stlat = lat(i)
        stlon = lon(i)
        call translat(stlat,stlat)
        call calthetaphi(eqlattmp,eqlon,stlat,stlon,theta(i),phi(i))
    enddo
    do i=1,nr
        read(lines(10+6*nzone+nr+i),'(a)') output(i)
        output(i)=trim(output(i))
    enddo
    return
   end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calthetaphi(ievla,ievlo,istla,istlo,theta,phi)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use parameters
    implicit none
    double precision,intent(in):: ievla,ievlo,istla,istlo
    double precision,intent(out):: theta,phi
    double precision:: evla,evlo,stla,stlo
    double precision:: gcarc,az,tc,ts

    ! transformation to spherical coordinates
    evla = 90 - ievla
    stla = 90 - istla

    evla = evla / 180 * pi
    evlo = ievlo / 180 * pi
    stla = stla / 180 * pi
    stlo = istlo / 180 * pi

    gcarc = dacos(dcos(evla) * dcos(stla)+dsin(evla)*dsin(stla)*dcos(evlo-stlo))

    tc = (dcos(stla)*dsin(evla)-dsin(stla)*dcos(evla)*dcos(stlo-evlo))/dsin(gcarc)
    ts = dsin(stla) * dsin(stlo - evlo) / dsin(gcarc)

    az = dacos(tc)
    if( ts < 0 ) az = - az

    az = az * 180 / pi

    gcarc = gcarc * 180 / pi

    theta = gcarc
    phi   = 180 - az
    return
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine translat(geodetic,geocentric)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use parameters
    implicit none
    double precision,intent(in):: geodetic
    double precision,intent(out):: geocentric
    double precision:: tmp_geodetic
    integer:: flag

    flag = 0
    tmp_geodetic=geodetic
    if(90<geodetic ) then
        tmp_geodetic = 180 - geodetic
        flag = 1
    endif

    tmp_geodetic = tmp_geodetic / 180 * pi
    geocentric = datan((1-flattening)*(1-flattening)*dtan(tmp_geodetic))
    geocentric = geocentric * 180 / pi
    !      if(geocentric < 0.d0 ) geocentric = 1.8d2 + geocentric
    if(flag == 1) geocentric = 180 - geocentric

    return
    end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calnl( nzone,vs,iphase,nsl,nll )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! counting of nsl and nll.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: nzone,iphase(*),nsl,nll
    double precision:: vs(4,*)
    integer:: i

    nsl = 0
    nll = 0
    do   i=1,nzone
        if ( ( vs(1,i)==0.d0 ).and.( vs(2,i)==0.d0 ).and.&
            ( vs(3,i)==0.d0 ).and.( vs(4,i)==0.d0 ) ) then
            nll = nll + 1
            iphase(i) = 2
        else
            nsl = nsl + 1
            iphase(i) = 1
        endif
    enddo
    return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calgrid( nzone,vrmin,vrmax,vp,vs,rmin,rmax,&
    imax,lmin,tlen,vmin,gridpar,dzpar )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use parameters
    implicit none
    integer:: nzone,imax,lmin
    double precision:: vrmin(*),vrmax(*),vp(4,*),vs(4,*)
    double precision:: rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
    integer:: izone,i,j
    double precision:: coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp

    do izone=1,nzone
        ! computing the S-velocity at each zone
        if ( vs(1,izone)==0.d0 ) then
            v(1:4) = vp(1:4,izone)
        else
            v(1:4) = vs(1:4,izone)
        endif
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
    gtmp = 0.d0
    do izone=1,nzone
        gtmp = gtmp + gridpar(izone)
    enddo
    do izone=1,nzone
        if ( gridpar(izone)>0.d0 ) then
            gridpar(izone) = gridpar(izone) / gtmp
        else
            rh = vrmax(izone) - vrmin(izone)
            gridpar(izone) = rh / ( rmax - rmin ) * 0.1d0
        endif
    enddo
    ! re-rearangement of gridpar
    gtmp = 0.d0
    do izone=1,nzone
        gtmp = gtmp + gridpar(izone)
    enddo
    do izone=1,nzone
        gridpar(izone) = gridpar(izone) / gtmp
    enddo
    !
    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calra(inlayer,jnlayer,jnslay,jnllay,&
    dzpar,nzone,vrmin,vrmax,iphase,&
    rmin,nslay,nllay,nnl,ra,re )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the number and the location of grid points.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use parameters
    implicit none
    integer:: inlayer,jnlayer,jnslay,jnllay
    integer:: nzone,iphase(*),nslay,nllay,nnl(maxnzone)
    double precision:: dzpar(*),vrmin(*),vrmax(*),rmin
    double precision:: ra(maxnlay+maxnzone+1)
    integer:: izone,itmp,i,ntmp
    double precision:: rh,re

    ! Initializing the data
    nslay = 0
    nllay = 0
    inlayer = 0
    ra(1:maxnlay+maxnzone+1)=0

    nnl(1:nzone)=0

    jnlayer = 0
    jnslay = 0
    jnllay = 0
    !
    !	tnlayer = nlayer / (2**(idr-1))
    ! computing the number and the location of the grid points
    ra(1) = rmin
    itmp = 1
    do izone=1,nzone
        rh = vrmax(izone) - vrmin(izone)
        if(dzpar(izone)==0.d0) then
            ntmp = 1
        else
            ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone)&
                / 2.d0 / pi  / 7.d-1 + 1 )
        endif
        ! ntmp (see Geller & Takeuchi 1995 6.2)
        !	  nnl(izone) = dint( dble(nlayer) * gridpar(izone) ) + 1
        nnl(izone) = ntmp
        if ( nnl(izone)<5 ) nnl(izone)=5
        if ( iphase(izone)==1 ) nslay = nslay + nnl(izone)
        if ( nslay>maxnslay ) stop  'nslay is too large. (calra)'
        if ( iphase(izone)==2 ) nllay = nllay + nnl(izone)
        if ( nllay>maxnllay ) stop  'nllay is too large. (calra)'
        do  i=1,nnl(izone)
            itmp = itmp + 1
            if ( itmp>maxnlay ) stop  'nlay is too large. (calra)'
            ra(itmp) = vrmin(izone) + rh * dble(i) / dble( nnl(izone) )
        enddo
    enddo
    !
    ! recouting the total number of grid points
    inlayer = 0
    do  izone=1,nzone
        inlayer = inlayer + nnl(izone)
    enddo
    jnlayer = jnlayer + inlayer
    jnslay  = jnslay  + nslay
    jnllay  = jnllay  + nllay
    !
    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calsp( maxnzone,ndc,nsl,nll,&
    iphase,nlayer,nllay,&
    isp,jsp,ksp,issp,ilsp,lsp,jssp,&
    isdr,jsdr,ildr,jdr,kdr )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the stack points.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: maxnzone
    integer:: ndc,nsl,nll,iphase(*),nlayer(maxnzone)
    integer:: nllay
    integer:: isp(maxnzone),jsp(maxnzone),ksp(maxnzone)
    integer:: issp(maxnzone),ilsp(maxnzone)
    integer:: lsp(maxnzone),jssp(maxnzone)
    integer:: isdr,jsdr,ildr,jdr,kdr
    integer:: i,isl,ill
    !
    ! Initialization of the data
    isp = 0
    jsp = 0
    ksp = 0
    issp = 0
    ilsp = 0
    lsp  = 0
    jssp = 0

    isdr = 0
    jsdr = 0
    ildr = 0
    jdr = 0
    kdr = 0
    ! computation of isp,jsp,ksp,issp,ilsp,lsp
    isp(1)  = 1
    jsp(1)  = 1
    ksp(1)  = 1
    issp(1) = 1
    ilsp(1) = 1
    lsp(1)  = 1
    jssp(1) = 1
    isl = 0
    ill = 0
    do  i=1,ndc
        isp(i+1) = isp(i) + nlayer(i)
        if ( iphase(i)==1 ) then
            jsp(i+1) = jsp(i) + 16 * nlayer(i)
            ksp(i+1) = ksp(i) + 2 * ( nlayer(i) + 1 )
            lsp(i+1) = lsp(i) + 4 * nlayer(i)
            isl = isl + 1
            if ( isl/=nsl ) then
                issp(isl+1) = issp(isl) + 4 * nlayer(i)
                jssp(isl+1) = jssp(isl) + nlayer(i) + 1
            endif
        else
            jsp(i+1) = jsp(i) + 4 * nlayer(i)
            ksp(i+1) = ksp(i) + ( nlayer(i) + 1 )
            lsp(i+1) = lsp(i) + 2 * nlayer(i)
            ill = ill + 1
            if ( ill/=nll )  ilsp(ill+1) = ilsp(ill) + 4 * nlayer(i)
        endif
    enddo
    isdr = 0
    jsdr = 0
    ildr = 0
    jdr = 0
    isdr = isdr + issp(nsl)-1 + 4 * nlayer(ndc+1)
    jsdr = jsdr + jssp(nsl)-1 + nlayer(ndc+1) + 1
    ildr = ildr + 4 * nllay
    jdr =  jdr  + jsp(ndc+1)-1 + 16 * nlayer(ndc+1)
    kdr =  kdr + ksp(ndc+1)-1 + 2 * ( nlayer(ndc+1)+1 )
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calspo( maxnlay,maxnzone,&
    rdc,iphase,inlayer,&
    r0,rmin,rmax,ra,isp,spo,spn )
    implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the source location.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: maxnlay,maxnzone,iphase(*)
    integer:: inlayer,isp(maxnzone),spn
    double precision:: rdc(*),r0,rmin,rmax,ra(maxnlay+maxnzone+1),spo
    integer:: itmp
    !
    ! checking the parameter
    if ( (r0<rmin).or.(r0>rmax) ) stop 'The source location is improper.(calspo)'
    spo = 0
    ! computing 'spo'
    if ( r0==rmax ) then
        spo = dble( inlayer ) - 0.01d0
        r0 = ra(inlayer) + (spo-dble(inlayer-1)) * ( ra(inlayer+1) -ra(inlayer) )
    else
        itmp = 2
110 continue
        if ( r0<ra(itmp) ) then
            continue
        else
           itmp = itmp + 1
            goto 110
        endif
        spo = dble(itmp-2)+ ( r0-ra(itmp-1) )   / ( ra(itmp)-ra(itmp-1) )
    ! temporal handling
        if ( (spo-dble(itmp-2))<0.01d0 ) then
            spo = dble(itmp-2) + 0.01d0
            r0 = ra(itmp-1)+ (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
        endif
        if ( (spo-dble(itmp-2))>0.99d0 ) then
            spo = dble(itmp-2) + 0.99d0
            r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
        endif
    endif
! computing 'spn'
    spn = 0
    itmp = 1
130 continue
    if ( iphase(itmp)==1 ) then
        spn = spn + 1
        if ( r0<=rdc(itmp) ) then
          continue
        else
            itmp = itmp + 1
            goto 130
        endif
    else
        spn = spn + 1
        if ( r0<=rdc(itmp) ) stop 'The source is in the liquid layer.(calspo)'
        itmp = itmp + 1
        goto 130
    endif
! changing 'spo'
    spo = spo - dble( isp(spn) - 1 )

    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calstg( maxnlay,maxnzone,nzone,rrho,&
    vpv,vph,vsv,vsh,eta,nnl,ra,rmax,&
    vnp,vra,rho,kappa,ecKx,ecKy,ecKz,&
    mu,ecL,ecN, r0,spn,ecC0,ecF0,ecL0 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the structure grid points.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: maxnlay,maxnzone,nzone,nnl(*),vnp,spn
    double precision:: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
    double precision:: ra(*),rmax
    double precision:: vra(*),rho(*),kappa(*),ecKx(*),ecKy(*),ecKz(*)
    double precision:: mu(*),ecL(*),ecN(*)
    double precision:: ecA,ecC,ecF
    double precision:: r0,ecA0,ecC0,ecF0,ecL0
    double precision:: trho,tvpv,tvph,tvsv,tvsh,teta,coef
    integer:: izone,i,j,itmp,jtmp
    !
    ! initializing the data
    vra(1: maxnlay+2*maxnzone+1)=0
    rho(1: maxnlay+2*maxnzone+1)=0
    kappa(1: maxnlay+2*maxnzone+1)=0
    ecKx(1: maxnlay+2*maxnzone+1)=0
    ecKy(1: maxnlay+2*maxnzone+1 )=0
    ecKz(1: maxnlay+2*maxnzone+1 )=0
    mu(1: maxnlay+2*maxnzone+1 )=0
    ecL(1: maxnlay+2*maxnzone+1 )=0
    ecN(1: maxnlay+2*maxnzone+1 )=0
    ! computing the structure grid points
    itmp = 0
    jtmp = 0
    do  izone=1,nzone
        do  i=1,nnl(izone)+1
            itmp = itmp + 1
            jtmp = jtmp + 1
            vra(itmp) = ra(jtmp)
            ! --- evaluating the density and elastic constants at this point
            trho = 0.d0
            tvpv = 0.d0
            tvph = 0.d0
            tvsv = 0.d0
            tvsh = 0.d0
            teta = 0.d0
            do j=1,4
                if ( j==1 ) then
                    coef = 1.d0
                else
                    coef = coef * ( vra(itmp) / rmax )
                endif
                trho  = trho  + rrho(j,izone)  * coef
                tvpv  = tvpv  + vpv(j,izone)   * coef
                tvph  = tvph  + vph(j,izone)   * coef
                tvsv  = tvsv  + vsv(j,izone)   * coef
                tvsh  = tvsh  + vsh(j,izone)   * coef
                teta  = teta  + eta(j,izone)   * coef
            enddo
            rho(itmp) = trho
            ecL(itmp)  = rho(itmp) * tvsv * tvsv
            ecN(itmp)  = rho(itmp) * tvsh * tvsh
            ecA = trho * tvph * tvph
            ecC = trho * tvpv * tvpv
            ecF = teta * ( ecA - 2.d0 * ecL(itmp) )
            kappa(itmp) = ( 4.d0 * ecA + ecC     + 4.d0 * ecF - 4.d0 * ecN(itmp) ) / 9.d0
            ecKx(itmp) = ecA - 4.d0 / 3.d0 * ecN(itmp)
            ecKy(itmp) = ecF + 2.d0 / 3.d0 * ecN(itmp)
            ecKz(itmp) = ( ecC + 2.d0 * ecF ) / 3.d0
        enddo
        jtmp = jtmp - 1
    enddo
    vnp = itmp
    !
    trho = 0.d0
    tvpv = 0.d0
    tvph = 0.d0
    tvsv = 0.d0
    tvsh = 0.d0
    teta = 0.d0
    do j=1,4
        if ( j==1 ) then
            coef = 1.d0
        else
            coef = coef * ( r0 / rmax )
        endif
        trho  = trho  + rrho(j,spn)  * coef
        tvpv  = tvpv  + vpv(j,spn)   * coef
        tvph  = tvph  + vph(j,spn)   * coef
        tvsv  = tvsv  + vsv(j,spn)   * coef
        tvsh  = tvsh  + vsh(j,spn)   * coef
        teta  = teta  + eta(j,spn)   * coef
    enddo
    ecL0 = trho * tvsv * tvsv
    ecA0 = trho * tvph * tvph
    ecC0 = trho * tvpv * tvpv
    ecF0 = teta * ( ecA0 - 2.d0 * ecL0 )
    !
    return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine caltstg( maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,nnl,ra,rmax,tvra,tkappa,tecKx,tecKy,tecKz,tmu,tecL,tecN)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the structure grid points.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: maxnlay,maxnzone,nzone,nnl(*)
    double precision:: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
    double precision:: ra(*),rmax
    double precision:: tvra(*),tkappa(*),tmu(*)
    double precision:: tecKx(*),tecKy(*),tecKz(*),tecL(*),tecN(*)
    double precision:: trho,tvpv,tvph,tvsv,tvsh,teta,coef
    double precision:: ecA,ecC,ecF
    integer:: izone,i,j,itmp,jtmp
    !
    tvra(1:maxnlay+2*maxnzone+1)=0
    tkappa(1:maxnlay+2*maxnzone+1)=0
    tecKx(1:maxnlay+2*maxnzone+1)=0
    tecKy(1:maxnlay+2*maxnzone+1)=0
    tecKz(1:maxnlay+2*maxnzone+1)=0
    tmu(1:maxnlay+2*maxnzone+1)=0
    tecL(1:maxnlay+2*maxnzone+1)=0
    tecN(1:maxnlay+2*maxnzone+1)=0
    ! computing the structure grid points
    itmp = 0
    jtmp = 0
    do izone=1,nzone
        do i=1,nnl(izone)+1
            itmp = itmp + 1
            jtmp = jtmp + 1
            tvra(itmp) = ra(jtmp)
            ! --- evaluating the density and elastic constants at this point
            trho = 0.d0
            tvpv = 0.d0
            tvph = 0.d0
            tvsv = 0.d0
            tvsh = 0.d0
            teta = 0.d0
            do j=1,4
                if ( j==1 ) then
                    coef = 1.d0
                else
                    coef = coef * ( tvra(itmp) / rmax )
                endif
                trho = trho + rrho(j,izone) * coef
                tvpv  = tvpv  + vpv(j,izone)   * coef
                tvph  = tvph  + vph(j,izone)   * coef
                tvsv  = tvsv  + vsv(j,izone)   * coef
                tvsh  = tvsh  + vsh(j,izone)   * coef
                teta  = teta  + eta(j,izone)   * coef
            enddo
            tecL(itmp)  = trho * tvsv * tvsv
            tecN(itmp)  = trho * tvsh * tvsh
            ecA = trho * tvph * tvph
            ecC = trho * tvpv * tvpv
            ecF = teta * ( ecA - 2.d0 * tecL(itmp) )
            tkappa(itmp) = ( 4.d0 * ecA + ecC  + 4.d0 * ecF - 4.d0 * tecN(itmp) )/ 9.d0
            tecKx(itmp) = ecA - 4.d0 / 3.d0 * tecN(itmp)
            tecKy(itmp) = ecF + 2.d0 / 3.d0 * tecN(itmp)
            tecKz(itmp) = ( ecC + 2.d0 * ecF ) / 3.d0
        enddo
        jtmp = jtmp - 1
    enddo
    return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calspdr( maxnzone,nzone,iphase,nlayer,jjdr,kkdr )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: maxnzone,nzone,iphase(*)
    integer:: nlayer(maxnzone),jjdr(*),kkdr(*)
    integer:: izone
    jjdr(1) = 1
    kkdr(1) = 1
    do  izone=1,nzone-1
        if ( iphase(izone)==1 ) then
            jjdr(izone+1) = jjdr(izone) + 16 * nlayer(izone)
            if ( iphase(izone+1)==1 ) then
                kkdr(izone+1)  = kkdr(izone) + 2 * nlayer(izone)
            else
                kkdr(izone+1) = kkdr(izone) + 2 * ( nlayer(izone)+1 )
            endif
        else
            jjdr(izone+1) = jjdr(izone) + 4 * nlayer(izone)
            if ( iphase(izone+1)==1 ) then
                kkdr(izone+1)   = kkdr(izone) + ( nlayer(izone)+1 )
            else
                kkdr(izone+1)  = kkdr(izone) + nlayer(izone)
            endif
        endif
    enddo
    return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calmdr( omega,l,nzone,vrmax,vmin,rmax,sufzone )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use parameters
    implicit none
    integer:: l,nzone,sufzone
    double precision:: omega,vrmax(*),vmin(*),rmax
    integer:: izone
    double precision:: gtmp,tdzpar
    !
    sufzone = 0
    do  izone=1,nzone
        gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) )&
            - ( (dble(l)+0.5d0) * (dble(l)+0.5d0) )&
            / ( vrmax(izone) * vrmax(izone) )
        if ( gtmp>0.d0 ) then
            tdzpar = dsqrt( 1.d0/gtmp )
        else
            if ( vrmax(izone)>rmax*(1-2.d0*pi/(dble(l)+0.50)) ) then
                tdzpar = 0.d0
            else
                sufzone = izone
                tdzpar = 0.d0
            endif
        endif
    enddo
    !
    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calamp( g,l,lsuf,maxamp,ismall,ratl )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use parameters
    implicit none
    integer:: ismall,l,lsuf
    double precision:: maxamp,ratl
    complex(dp) g(2)
    double precision:: amp,ampratio

    ampratio = 0.d0
    amp = dsqrt( cdabs( g(1) )**2 + cdabs( g(2) )**2 )
    if ( amp>maxamp ) maxamp = amp
    if ( (amp/=0.d0).and.(maxamp/=0.d0) )   ampratio = amp / maxamp
    if ( ( ampratio<ratl ).and.( l>=lsuf ) ) then
        ismall = ismall + 1
    else
        ismall = 0
    endif
    return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calcutd(nzone,nnl,tmpc,rat,nn,iphase,ra,kkdr,kc)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use parameters
    implicit none
    integer:: nzone,nn,kkdr(*),kc,iphase(*),nnl(*)
    complex(dp) ::tmpc(*)
    double precision:: rat,ra(*)
    integer:: nc
    double precision:: cU(nn),cV(nn),rc
    double precision:: maxamp,amp(nn)
    integer:: iz,jz,jj,i,ml(nzone),tzone

    cU = 0
    cV = 0

    iz = 2
    jz = 1
    do  jj=1,nn
        if(jj==kkdr(iz)) then
            if(iphase(iz)/=iphase(iz-1)) jz = jz - 1
            iz = iz + 1
        endif
        if(iphase(iz-1)==1) then
            if(mod((jj-kkdr(iz-1)),2)==1) then ! U
                cU(jz) = cdabs(tmpc(jj))
                jz = jz + 1
                ! V nothing
            endif
        else ! U in fluid
            cU(jz) = cdabs(tmpc(jj))
            jz = jz + 1
        endif
    enddo
    !
    maxamp = -1.d0
    do i=1,jz-1
        amp(i) = cU(i)
        if(maxamp<amp(i)) maxamp = amp(i)
    enddo
    !
    maxamp = maxamp * rat ! threshold value
    !
    nc = 1
    do i=1,jz-1
        if(amp(i)>maxamp) then
            nc = i
            goto 140
        endif
    enddo
140 continue

    i = 1
    do jj=1,nzone
        i = i + nnl(jj)
        ml(jj) = i
    enddo
    !
    do jj=nzone,1,-1
        if(ml(jj)>nc) tzone = jj
    enddo

    rc = ra(nc)

    do i=1,jz-1
        if( (ra(i)<=rc).and.(rc<ra(i+1)) ) then
            nc = i
            if(tzone==1) then ! case(tzone is innermost zone)
                if(iphase(tzone)==1) kc = 1 + 2 * nc
                if(iphase(tzone)==2) kc = 1 + nc
            else
                if(iphase(tzone)==1) then
                    kc = kkdr(tzone) + 2 * (nc - ml(tzone-1))
                endif
                if(iphase(tzone)==2) then
                    kc = kkdr(tzone) + nc - ml(tzone-1)
                endif
            endif
        endif
    enddo

    return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nzone,lsuf
    double precision:: omega,vrmax(*),vsv(4,*)

    double precision:: tvs,coef
    integer:: i

    tvs = 0.d0
    do i=1,4
        if(i==1) then
            coef = 1.d0
        else
            coef = coef
        endif
        tvs = tvs + ( vsv(i,nzone) ) * coef
    enddo

    lsuf = int(omega * vrmax(nzone) / tvs - 0.5d0) + 1
    return
end
!
