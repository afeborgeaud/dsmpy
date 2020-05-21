! calmat.f for wcalprem.f
!----------------------------------------------------------
subroutine calmatc(nlayer,vnp,vra,con,rpow,w1dn,w2dn,ra,m,work)
!----------------------------------------------------------
! Computing \int r^rpow con W_p^(w1dn) W_q^(w2dn) dr.
!----------------------------------------------------------
    implicit none
    integer:: maxrpow
    parameter ( maxrpow = 2 )
    integer,intent(in):: nlayer,vnp,rpow,w1dn,w2dn
    double precision:: vra(vnp),con(vnp),ra(*),m(*),work(*)
    integer:: i,j,k,l,nn,snp
    double precision:: a(2,2),b(2,2),c(5),rh
    ! parameter check
    if (maxrpow<rpow) stop 'Invalid arguments.(calmatc)'
    ! computing of the matrix elements
    snp = 1
    do i=1,nlayer
        rh = ra(i+1) - ra(i)
        if ( w1dn == 0 ) then
            a(2,1) = -1.d0 / rh
            a(1,1) = ra(i+1) / rh
            a(2,2) = 1.d0 / rh
            a(1,2) = -ra(i) / rh
        elseif ( w1dn == 1 ) then
                a(2,1) = 0.d0
                a(1,1) = -1.d0 / rh
                a(2,2) = 0.d0
                a(1,2) = 1.d0 / rh
        else
            stop 'Invalid arguments.(calmatc)'
        endif

        if ( w2dn == 0 ) then
            b(2,1) = -1.d0 / rh
            b(1,1) = ra(i+1) / rh
            b(2,2) = 1.d0 / rh
            b(1,2) = -ra(i) / rh
        elseif ( w2dn == 1 ) then
                b(2,1) = 0.d0
                b(1,1) = -1.d0 / rh
                b(2,2) = 0.d0
                b(1,2) = 1.d0 / rh
        else
            stop 'Invalid arguments.(calmatc)'
        endif

        do j=1,2
            do k=1,2
                c(1:5) = 0.d0

                call pmulti( 2,a(1,j),2,b(1,k),3,c )
                do l=3,1,-1
                    c(l+rpow) = c(l)
                    if ( rpow>0 ) c(l)=0.d0
                enddo
                nn = 4 * (i-1) + 2 * (j-1) + k
                call pinteg(snp,5,c,ra(i),ra(i+1),vnp,vra,con,work(nn))
            enddo
        enddo
    enddo
    if ( w1dn /= w2dn ) then
        m(1:nlayer*4-3:4)=2*work(1:nlayer*4-3:4)
        m(2:nlayer*4-2:4)=work(2:nlayer*4-2:4)+work(3:nlayer*4-1:4)
        m(3:nlayer*4-1:4)=work(2:nlayer*4-2:4)+work(3:nlayer*4-1:4)
        m(4:nlayer*4:4)=2*work(4:nlayer*4:4)
    else
        m(1:4*nlayer) = work(1:4*nlayer)
    endif
    return
end

!----------------------------------------------------------
subroutine caltl( nlayer,vnp,vra,rho,ra,tl)
!----------------------------------------------------------
! Computing of lumped mass matrix.
!----------------------------------------------------------
    implicit none
    integer:: nlayer,vnp
    double precision:: vra(vnp),rho(vnp),ra(*),tl(*)
    integer:: i,nn,snp
    double precision:: c(3),from,to

    snp = 1
    do i=1,nlayer
        c(1:2) = 0.d0
        c(3) = 1.d0

        from = ra(i)
        to = ( ra(i) + ra(i+1) ) / 2.d0
        nn = 4 * (i-1)
        call pinteg( snp,3,c,from,to,vnp,vra,rho,tl(nn+1) )

        tl(nn+2:nn+3) = 0.d0

        from = to
        to = ra(i+1)
        call pinteg( snp,3,c,from,to,vnp,vra,rho,tl(nn+4) )

    enddo
    return
end

!----------------------------------------------------------
subroutine calhl( nlayer,vnp,vra,mu,ra,hl)
!----------------------------------------------------------
! Computing of lumped rigidity matrix.
!----------------------------------------------------------
    implicit none
    integer:: nlayer,vnp
    double precision:: vra(vnp),mu(vnp),ra(*),hl(*)
    integer:: i,nn,snp
    double precision:: c(1),from,to

    snp = 1
    do  i=1,nlayer
        c(1) = 1.d0

        from = ra(i)
        to = ( ra(i) + ra(i+1) ) / 2.d0
        nn = 4 * (i-1)
        call pinteg( snp,1,c,from,to,vnp,vra,mu,hl(nn+1) )

        hl(nn+2:nn+3) = 0.d0

        from = to
        to = ra(i+1)
        call pinteg( snp,1,c,from,to,vnp,vra,mu,hl(nn+4) )

    enddo
    return
end

!----------------------------------------------------------
subroutine pmulti(n,a,m,b,l,c)
!----------------------------------------------------------
! Computing the (l-1) degrees polynomial c(n) which is the product of
! (n-1) degrees polynomial a(n) and (m-1) degrees polynomial b(n).
!----------------------------------------------------------
    integer:: n,m,l,i,j
    double precision:: a(n),b(m),c(l)

    if (n+m-1/=l) stop 'Invalid arguments.(pmulti)'

    c(1:l)=0.d0

    do i=1,n
        do j=1,m
            c(i+j-1) = c(i+j-1) + a(i)*b(j)
        enddo
    enddo

    return
end

!----------------------------------------------------------
subroutine pinteg(snp,n,p,from,to,vnp,vra,con,pint)
!----------------------------------------------------------
! Evaluating the integrated value pint from 'from' to 'to' of p(n)*con
! which is the product of (n-1) degrees polynomial 'p(n)' and 'con'.
!----------------------------------------------------------
    implicit none
    integer:: maxn
    parameter ( maxn = 5 )

    integer:: snp,n,vnp
    double precision:: from,to,p(n),vra(vnp),con(vnp),pint
    double precision:: x1,x2,q(2),pq(maxn+1),psint

    if ( maxn<n )stop 'Degrees of a polynomial is too large.(pinteg)'
    if ( snp>=vnp ) snp = 1

    pint = 0.d0
    x1 = from
    do
        if (  vra(snp)<=x1 .and. vra(snp+1)>x1  ) then
            x2 = dmin1( to, vra(snp+1) )
        else
            snp = snp + 1
            if ( snp==vnp ) snp = 1
            cycle
        endif
        ! evaluating the integrated value
        q(2) = ( con(snp+1)-con(snp) ) / ( vra(snp+1)-vra(snp) )
        q(1) = - q(2) * vra(snp) + con(snp)
        call pmulti( n,p,2,q,n+1,pq )
        call polint( n+1,pq,x1,x2,psint )
        pint = pint + psint

        if ( x2==to ) exit
        x1 = x2
    enddo
    return
end

!----------------------------------------------------------
subroutine polint( n,p,x1,x2,pint )
!----------------------------------------------------------
!c Evaluating the integrated value 'pint' from 'x1' to 'x2'
!c of (n-1) degrees polynomial 'p(n)'.
!----------------------------------------------------------
    implicit none
    integer,parameter:: maxn=6
    integer:: n
    double precision:: p(*),x1,x2,pint
    integer:: i,j
    double precision:: a(maxn),b(maxn),dx,xx

    if (maxn<n) stop 'Degrees of a polynomial is too large.(polint)'

    a(1) = 1.d0
    b(1) = 1.d0
    if ( n>=2 ) then
        do i=2,n
            a(i) = a(i-1) * x1
            b(i) = b(i-1) * x2
        enddo
    endif
    dx = x2 - x1

    pint = 0.d0
    do i=1,n
        xx = 0.d0
        do j=1,i
            xx = xx + a(j) * b(i-j+1) / dble(i)
        enddo
        pint = pint + p(i) * dx * xx
    enddo
    return
end

!----------------------------------------------------------
subroutine overlap( nlayer,a,a2 )
!----------------------------------------------------------
! Overlapping the coefficient matrix elements in the solid part.
!----------------------------------------------------------
    use parameters
    implicit none
    integer,intent(in):: nlayer
    complex(dp):: a(*),a2(2,*)
    integer:: i,j,k,mu,m,i1,i2,k1,k2,nsize

    nsize = nlayer+1
    mu = 1
    m = mu + 1
    do j=1,nsize
        i1 = max0(1,j-mu)
        i2 = j
        do i=i1,i2
            k = i - j + m
            if ( i==j ) then
                if ( i==1 ) then
                    a2(k,j) = a2(k,j) + a(1)
                else
                    if ( i==nsize ) then
                        a2(k,j) = a2(k,j) + a( 4*nlayer )
                    else
                        k1 = 4 * i - 4
                        k2 = k1 + 1
                        a2(k,j) = a2(k,j) + a(k1) + a(k2)
                    endif
                endif
            endif
            if (i+1==j) then
                k1 = 4 * i - 2
                a2(k,j) = a2(k,j) + a(k1)
            endif
        enddo
    enddo
    return
end

!----------------------------------------------------------
subroutine calg2( l,m,spo,r0,mt,mu0,coef,ga,a,ga2,dr,g2 )
!----------------------------------------------------------
    use parameters
    implicit none

    integer:: l,m
    double precision:: spo,r0,mt(3,3),mu0
    complex(dp) coef,ga(*),a(*),ga2(2,*),g2(*)
    integer:: i,itmp
    double precision:: b1,sgn,eps,ier
    complex(dp) dd,cg2(3),dr(3),z(3)

    data eps/ -1.d0 /

    ! computing of particular solution with free surface boundary conditions
    cg2=0
    dd = 0
    if ( m>=0 ) then
        sgn = 1.d0
    else
        sgn = - 1.d0
    endif
    if ( iabs(m)==1 ) then
        b1 = dsqrt( dble( 2*l+1 ) / ( 16.d0 * pi ) )
        !c
        dd = dcmplx( b1 ) * ( dcmplx( sgn * mt(1,3), mt(1,2) ) )/ ( dcmplx( r0 * r0 * mu0 ) * coef )
        itmp = 4
        cg2=0
        do  i=2,3
            cg2(i) = - dd * ( ga(itmp+1) + ga(itmp+2) )
            itmp = itmp + 2
        enddo
    else
        if ( iabs(m)==2 ) then
            b1 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2)/(64.d0*pi) )
            cg2(2) = dcmplx( b1 / r0 ) * dcmplx( 2.d0*mt(2,3), sgn*( mt(2,2)-mt(3,3) ) )
        endif
    endif
    if ( m==-2 .or. m==-l ) then
        call dclisb0( ga2,3,1,2,cg2,eps,dr,z,ier)
    else
        call dcsbsub0( ga2,3,1,2,cg2,eps,dr,z,ier)
    endif
    cg2(3) = cg2(3) + dd
    ! computation of the excitation vector
    itmp = dint(spo)
    g2(itmp+1) = a(1) * cg2(1) + a(2) * cg2(3)
    g2(itmp+2) = a(3) * cg2(1) + a(4) * cg2(3)

    return
end
