! calmat.f for wcalprem.f
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calmatc( nlayer,vnp,vra,con,rpow,w1dn,w2dn,ra,m )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing \int r^rpow con W_p^(w1dn) W_q^(w2dn) dr.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer,parameter:: maxrpow= 2
    integer:: nlayer,vnp,rpow,w1dn,w2dn
    double precision:: vra(vnp),con(vnp),ra(*),m(*)
    integer:: i,j,k,kk,l,nn,snp
    double precision:: a(2,2),b(2,2),c(5),rh
    ! parameter check
    if (maxrpow<rpow ) stop 'Invalid arguments.(calmatc)'
    ! computing of the matrix elements
    snp = 1
    do i=1,nlayer
        rh = ra(i+1) - ra(i)
        if ( w1dn == 0 ) then
            a(2,1) = -1.d0 / rh
            a(1,1) = ra(i+1) / rh
            a(2,2) = 1.d0 / rh
            a(1,2) = -ra(i) / rh
        else
            if ( w1dn == 1 ) then
                a(2,1) = 0.d0
                a(1,1) = -1.d0 / rh
                a(2,2) = 0.d0
                a(1,2) = 1.d0 / rh
            else
                stop 'Invalid arguments.(calmatc)'
            endif
        endif
        if ( w2dn == 0 ) then
            b(2,1) = -1.d0 / rh
            b(1,1) = ra(i+1) / rh
            b(2,2) = 1.d0 / rh
            b(1,2) = -ra(i) / rh
        else
            if ( w2dn == 1 ) then
                b(2,1) = 0.d0
                b(1,1) = -1.d0 / rh
                b(2,2) = 0.d0
                b(1,2) = 1.d0 / rh
            else
                stop 'Invalid arguments.(calmatc)'
            endif
        endif
        do j=1,2
            do k=1,2
                do kk=1,5
                    c(kk) = 0.d0
                enddo
                call pmulti( 2,a(1,j),2,b(1,k),3,c )
                do l=3,1,-1
                    c(l+rpow) = c(l)
                    if ( rpow>0 ) c(l)=0.d0
                enddo
                nn = 4 * (i-1) + 2 * (j-1) + k
                call pinteg( snp,5,c,ra(i),ra(i+1),&
                    vnp,vra,con,m(nn) )
            enddo
        enddo
    enddo
    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine pmulti(n,a,m,b,l,c)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the (l-1) degrees polynomial c(n) which is the product of
! (n-1) degrees polynomial a(n) and (m-1) degrees polynomial b(n).
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: n,m,l,i,j
    double precision:: a(n),b(m),c(l)

    if (n+m-1/=l) stop 'Invalid arguments.(pmulti)'
    c=0
    do i=1,n
        do j=1,m
            c(i+j-1) = c(i+j-1) + a(i)*b(j)
        enddo
    enddo

    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine pinteg(snp,n,p,from,to,vnp,vra,con,pint)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Evaluating the integrated value pint from 'from' to 'to' of p(n)*con
! which is the product of (n-1) degrees polynomial 'p(n)' and 'con'.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer,parameter::maxn = 5
    integer:: snp,n,vnp
    double precision:: from,to,p(n),vra(vnp),con(vnp),pint
    double precision:: x1,x2,q(2),pq(maxn+1),psint

    if ( maxn<n ) stop 'Degrees of a polynomial is too large.(pinteg)'
    if ( vnp <=snp ) snp = 1

    pint = 0.d0
    x1 = from
100 continue
    if ( (vra(snp)<=x1).and.(vra(snp+1)>x1) ) then
        x2 = dmin1( to, vra(snp+1) )
    else
        snp = snp + 1
        if ( snp==vnp ) snp = 1
        goto 100
    endif
    ! evaluating the integrated value
    q(2) = ( con(snp+1)-con(snp) ) / ( vra(snp+1)-vra(snp) )
    q(1) = - q(2) * vra(snp) + con(snp)
    call pmulti( n,p,2,q,n+1,pq )
    call polint( n+1,pq,x1,x2,psint )
    pint = pint + psint

    if ( x2==to ) then
        continue
    else
        x1 = x2
        goto 100
    endif
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine polint( n,p,x1,x2,pint )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Evaluating the integrated value 'pint' from 'x1' to 'x2'
! of (n-1) degrees polynomial 'p(n)'.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer,parameter :: maxn = 6
    integer:: i,j,n
    double precision:: p(*),x1,x2,pint
    double precision:: a(maxn),b(maxn),dx,xx

    if (maxn<n ) stop 'Degree of a polynomial is too large.(polint)'

    a(1) = 1.d0
    b(1) = 1.d0
    if ( 2<=n ) then
        do  i=2,n
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine caltl( nlayer,vnp,vra,rho,ra,tl)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing of lumped mass matrix.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer,vnp
    double precision:: vra(vnp),rho(vnp),ra(*),tl(*)
    integer:: i,nn,snp
    double precision:: c(3),from,to
    !
    snp = 1
    do i=1,nlayer
        c(1) = 0.d0
        c(2) = 0.d0
        c(3) = 1.d0

        from = ra(i)
        to = ( ra(i) + ra(i+1) ) / 2.d0
        nn = 4 * (i-1)
        call pinteg( snp,3,c,from,to,vnp,vra,rho,tl(nn+1) )

        tl(nn+2) = 0.d0
        tl(nn+3) = 0.d0

        from = to
        to = ra(i+1)
        call pinteg( snp,3,c,from,to,vnp,vra,rho,tl(nn+4) )
    enddo
    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calhl( nlayer,vnp,vra,mu,ra,hl)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing of lumped rigidity matrix.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer,vnp
    double precision:: vra(vnp),mu(vnp),ra(*),hl(*)
    integer:: i,nn,snp
    double precision:: c(1),from,to

    snp = 1
    do i=1,nlayer
        c(1) = 1.d0
        from = ra(i)
        to = ( ra(i) + ra(i+1) ) / 2.d0
        nn = 4 * (i-1)
        call pinteg( snp,1,c,from,to,vnp,vra,mu,hl(nn+1) )

        hl(nn+2) = 0.d0
        hl(nn+3) = 0.d0

        from = to
        to = ra(i+1)
        call pinteg( snp,1,c,from,to,vnp,vra,mu,hl(nn+4) )
    enddo
    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calt( nlayer, tl, tc, t )
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: t(*), tl(*), tc(*)
    integer:: i

    do i=1,4*nlayer
        t(i) = ( tl(i) + tc(i) ) / 2.d0
    enddo

    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calh5( nlayer,vra,con,ra,h5 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computation of the submatrix `h5'.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: vra(*),con(*),ra(*),h5(*)
    integer:: itmp,i
    !
    ! Data initialization
    h5(1:4*nlayer)=0

    itmp = 0
110 continue
    itmp = itmp + 1
    if ( vra(itmp)==ra(1) ) then
        if ( vra(itmp+1)==vra(itmp) ) itmp = itmp + 1
    else
        goto 110
    endif

    itmp = itmp - 1

    do i=1,nlayer
        itmp = itmp + 1
        h5(4*i-3) = - 3.d0 / 8.d0 * con(itmp) * ra(i) - 1.d0 / 8.d0 * con(itmp+1) * ra(i+1)
        h5(4*i-2) = - h5(4*i-3)
        h5(4*i-1) = - 1.d0 / 8.d0 * con(itmp) * ra(i) - 3.d0 / 8.d0 * con(itmp+1) * ra(i+1)
        h5(4*i)   = - h5(4*i-1)
    enddo

    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calhm1( nlayer,vra,con,ra,hm1 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the modified submatrix `hm1'.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: vra(*),con(*),ra(*),hm1(-1:2,*)
    integer:: itmp,i
    !
    itmp = 0
100 continue
    itmp = itmp + 1
    if ( vra(itmp)==ra(1) ) then
        if ( vra(itmp+1)==vra(itmp) ) itmp = itmp + 1
    else
        goto 100
    endif

    hm1(0,1) = hm1(0,1) - 7.d0 / 12.d0 * con(itmp) * ra(1)
    hm1(1,1) = hm1(1,1) + 8.d0 / 12.d0 * con(itmp) * ra(1)
    hm1(2,1) = hm1(2,1) - 1.d0 / 12.d0 * con(itmp) * ra(1)
    do  i=2,nlayer-1
        itmp = itmp + 1
        hm1(-1,i) = hm1(-1,i) - 5.d0 / 12.d0 * con(itmp) * ra(i)
        hm1( 0,i) = hm1( 0,i) - 3.d0 / 12.d0 * con(itmp) * ra(i)
        hm1( 1,i) = hm1( 1,i) + 9.d0 / 12.d0 * con(itmp) * ra(i)
        hm1( 2,i) = hm1( 2,i) - 1.d0 / 12.d0 * con(itmp) * ra(i)
    enddo
    itmp = itmp + 1
    hm1(-1,nlayer)   = hm1(-1,nlayer)      - 5.d0 / 12.d0 * con(itmp) * ra(nlayer)
    hm1( 0,nlayer)   = hm1( 0,nlayer)   - 3.d0 / 12.d0 * con(itmp) * ra(nlayer)
    hm1( 1,nlayer)   = hm1( 1,nlayer)     + 8.d0 / 12.d0 * con(itmp) * ra(nlayer)
    itmp = itmp + 1
    hm1(-1,nlayer+1) = hm1(-1,nlayer+1) - 5.d0 / 12.d0 * con(itmp) * ra(nlayer+1)
    hm1( 0,nlayer+1) = hm1( 0,nlayer+1)  + 5.d0 / 12.d0 * con(itmp) * ra(nlayer+1)
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calhm2( nlayer,vra,con,ra,hm2 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the modified submatrix `hm2'.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer ::nlayer
    double precision:: vra(*),con(*),ra(*),hm2(-2:1,*)
    integer ::itmp,i

    itmp = 0
100 continue
    itmp = itmp + 1
    if ( vra(itmp)==ra(1) ) then
        if ( vra(itmp+1)==vra(itmp) ) itmp = itmp + 1
    else
        goto 100
    endif
    !
    hm2(0,1) = hm2(0,1) - 5.d0 / 12.d0 * con(itmp) * ra(1)
    hm2(1,1) = hm2(1,1) + 5.d0 / 12.d0 * con(itmp) * ra(1)
    itmp = itmp + 1
    hm2(-1,2) = hm2(-1,2) - 8.d0 / 12.d0 * con(itmp) * ra(2)
    hm2( 0,2) = hm2( 0,2) + 3.d0 / 12.d0 * con(itmp) * ra(2)
    hm2( 1,2) = hm2( 1,2) + 5.d0 / 12.d0 * con(itmp) * ra(2)
    do  i=3,nlayer
        itmp = itmp + 1
        hm2(-2,i) = hm2(-2,i) + 1.d0 / 12.d0 * con(itmp) * ra(i)
        hm2(-1,i) = hm2(-1,i) - 9.d0 / 12.d0 * con(itmp) * ra(i)
        hm2( 0,i) = hm2( 0,i) + 3.d0 / 12.d0 * con(itmp) * ra(i)
        hm2( 1,i) = hm2( 1,i) + 5.d0 / 12.d0 * con(itmp) * ra(i)
    enddo
    itmp = itmp + 1
    hm2(-2,nlayer+1) = hm2(-2,nlayer+1)   + 1.d0 / 12.d0 * con(itmp) * ra(nlayer+1)
    hm2(-1,nlayer+1) = hm2(-1,nlayer+1)   - 8.d0 / 12.d0 * con(itmp) * ra(nlayer+1)
    hm2( 0,nlayer+1) = hm2( 0,nlayer+1)   + 7.d0 / 12.d0 * con(itmp) * ra(nlayer+1)
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mtrnp(nlayer,a1,a2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing 'a2' which is the transpose of 'a1'.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: a1(*),a2(*)
    integer:: i,nn

    nn = 4 * nlayer
    do i=1,nn-3,4
        a2(i) = a1(i)
    enddo
    do  i=2,nn-2,4
        a2(i) = a1(i+1)
    enddo
    do i=3,nn-1,4
        a2(i) = a1(i-1)
    enddo
    do i=4,nn,4
        a2(i) = a1(i)
    enddo

    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mtrnp2( nlayer,p,q,hm1,hm2 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer::nlayer,p,q
    double precision:: hm1(-p:q,*),hm2(-q:p,*)
    integer::i,j,m,n

    do i=-q,p
        do j=1,nlayer+1
            m = -i
            n = i + j
            if ( ( n>=1 ).and.( n<=nlayer+1 ) ) hm2(i,j) = hm1(m,n)
        enddo
    enddo
    !
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calcoef( nzone,omega,qmu,qkappa,coef1,coef2,coef )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nzone
    double precision:: omega,qmu(*),qkappa(*)
    complex(kind(0d0)):: coef1(*),coef2(*),coef(*)
    integer:: i
    double precision:: aa,bb
    double precision,parameter:: pi=3.1415926535897932d0

    ! evaluating the effect of anelasticity
    do  i=1,nzone
        if ( qmu(i)<=0.d0 ) then
            coef1(i) = dcmplx( 1.d0 )
        else
            if ( omega==0.d0 ) then
                aa = 1.d0
            else
                aa = 1.d0  + dlog( omega / ( 2.d0 * pi ) ) / ( pi * qmu(i) )
            endif
            bb = 1.d0 / ( 2.d0 * qmu(i) )
            coef1(i) = dcmplx( aa, bb ) * dcmplx( aa, bb )
        endif
        if ( qkappa(i)<=0.d0 ) then
            coef2(i) = dcmplx( 1.d0 )
            coef(i) = dcmplx( 1.d0 ) / coef2(i)
        else
            if ( omega==0.d0 ) then
                aa = 1.d0
            else
                aa = 1.d0  + dlog( omega / ( 2.d0 * pi ) ) / ( pi * qkappa(i) )
            endif
            bb = 1.d0 / ( 2.d0 * qkappa(i) )
            coef2(i) = dcmplx( aa, bb ) * dcmplx( aa, bb )
            coef(i) = dcmplx( 1.d0 ) / coef2(i)
        endif
    enddo
    return
    end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine cala0( nlayer,omega,omegai,t,h1x,h1y,h1z,h2L,h2N,&
    h3ax,h3ay,h3az,h4aL,h4aN,h5ax,h5ay,h5az,h6aL,h6aN,&
    h7x,h7y,h7z,h8L,h8N,coef1,coef2,a0 )
    implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the coefficient matrix 'a' in the solid part.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: nlayer
    double precision:: omega,omegai
    double precision:: t(*)
    double precision:: h1x(*),h1y(*),h1z(*),h2L(*),h2N(*)
    double precision:: h3ax(*),h3ay(*),h3az(*),h4aL(*),h4aN(*)
    double precision:: h5ax(*),h5ay(*),h5az(*),h6aL(*),h6aN(*)
    double precision:: h7x(*),h7y(*),h7z(*),h8L(*),h8N(*)
    complex(kind(0d0))::coef1,coef2,a0(*)
    integer:: i,nnn
    double precision:: tt
    complex(kind(0d0))::hh0,comega2
    !
    a0( 1:16*nlayer)=0
    comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
    do i=1,16*nlayer-3,4
        nnn = ( i + 3 ) / 4
        if ( mod(nnn,4)==3 ) cycle
        tt = t(nnn)
        hh0 = coef2 * dcmplx( 4.d0 * h1x(nnn) )& !4A
            + coef1 * dcmplx( 16.d0/3.d0 * h2N(nnn) )& !4A
            + coef1 * dcmplx( -4.d0 * h2N(nnn) )& !-4N
            + coef2 * dcmplx( 2.d0 * ( h3ay(nnn) + h5ay(nnn) ) )& !2F
            - coef1 * dcmplx( 4.d0/3.d0 * ( h4aN(nnn) + h6aN(nnn) ) )&
            + coef2 * dcmplx( 3.d0 *h7z(nnn) -2.d0*h7y(nnn) )& !C
            + coef1 * dcmplx( 4.d0/3.d0 * h8N(nnn) ) !C
        a0(i) = comega2 * dcmplx( tt ) - hh0
    enddo
    do i=4,16*nlayer,4
        nnn = i / 4
        if ( mod(nnn,4)==3 ) cycle
        tt = t(nnn)
        hh0 = coef1 * dcmplx(   h2L(nnn) ) &!L
            + coef1 * dcmplx( -2.d0 * h2N(nnn) )& !-2N
            - coef1 * dcmplx( h4aL(nnn) + h6aL(nnn) )& !-L
            + coef1 * dcmplx( h8L(nnn) ) !L
        a0(i) = comega2 * dcmplx( tt ) - hh0
    enddo
    return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine cala1( nlayer,h1x,h1y,h1z,h2L,h2N,h3x,h3y,h3z,&
    h4L,h4N,h5x,h5y,h5z,h6L,h6N,coef1,coef2,a1 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the coefficient matrix 'a' in the solid part.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: h1x(*),h1y(*),h1z(*),h2L(*),h2N(*)
    double precision:: h3x(*),h3y(*),h3z(*),h4L(*),h4N(*)
    double precision:: h5x(*),h5y(*),h5z(*),h6L(*),h6N(*)
    complex(kind(0d0)):: coef1,coef2,a1(*)
    integer:: i,nnn
    complex(kind(0d0)):: hh1

    a1( 1:16*nlayer )=0
    do i=2,16*nlayer-2,4
        nnn = ( i + 2 ) / 4
        if ( mod(nnn,4)==3 ) cycle
        hh1 = - coef2 * dcmplx( 2.d0 * h1x(nnn) )&  !2A
            - coef1 * dcmplx( 8.d0/3.d0 * h2N(nnn) ) &!2A
            - coef1 * dcmplx( h2L(nnn) )& !L
            - coef1 * dcmplx( -2.d0 * h2N(nnn) ) &!-2N
            - coef2 * dcmplx( h3y(nnn) )& !F
            - coef1 * dcmplx( -2.d0/3.d0 * h4N(nnn) ) &!F
            + coef1 * dcmplx( h6L(nnn) ) !-L
        a1(i) = - hh1
    enddo
    do i=3,16*nlayer-1,4
        nnn = ( i + 1 ) / 4
        if ( mod(nnn,4)/=2 ) cycle
        hh1 = - coef2 * dcmplx( 2.d0 * h1x(nnn) )& !2A
            - coef1 * dcmplx( 8.d0/3.d0 * h2N(nnn) ) & !2A
            - coef1 * dcmplx( -2.d0 * h2N(nnn) )& !-2N
            - coef1 * dcmplx(  h2L(nnn) ) & !L
            + coef1 * dcmplx( h4L(nnn) ) & !-L
            - coef2 * dcmplx( h5y(nnn) )& !F
            + coef1 * dcmplx( 2.d0/3.d0 * h6N(nnn) ) !F
        a1(i) = - hh1
    enddo
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine cala2( nlayer,h1x,h1y,h1z,h2L,h2N,coef1,coef2,a2 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the coefficient matrix 'a' in the solid part.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: h1x(*),h1y(*),h1z(*),h2L(*),h2N(*)
    complex(kind(0d0)):: coef1,coef2,a2(*)
    integer:: i,nnn
    complex(kind(0d0)):: hh2

    a2( 1:16*nlayer)=0
    do i=1,16*nlayer-3,4
        nnn = ( i + 3 ) / 4
        if ( mod(nnn,4)==3 ) cycle
        hh2 = coef1 * dcmplx( h2L(nnn) )
        a2(i) = - hh2
    enddo

    do i=4,16*nlayer,4
        nnn = i / 4
        if ( mod(nnn,4)==3 ) cycle
        hh2 = coef2 * dcmplx( h1x(nnn) )  + coef1 * dcmplx( 4.d0/3.d0 * h2N(nnn) )
        a2(i) = - hh2
    enddo
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calb0( nlayer,omega,omegai,p1,p3,coef,b0 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the coefficient matrix 'b' in the solid part.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: omega,omegai,p1(*),p3(*)
    complex(kind(0d0)):: coef,b0(*)
    integer:: i
    complex(kind(0d0)):: comega2

    b0( 1:4*nlayer )=0
    comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
    do i=1,4*nlayer-3,4
        b0(i) = - dcmplx( p1(i) ) / comega2  + coef * dcmplx( p3(i) )
    enddo
    do i=2,4*nlayer-2,4
        b0(i) = - dcmplx( p1(i) ) / comega2  + coef * dcmplx( p3(i) )
    enddo
    do i=4,4*nlayer,4
        b0(i) = - dcmplx( p1(i) ) / comega2     + coef * dcmplx( p3(i) )
    enddo
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calb2( nlayer,omega,omegai,p2,coef,b2 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the coefficient matrix 'b' in the solid part.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: omega,omegai,p2(*)
    complex(kind(0d0)):: coef,b2(*)
    integer:: i
    complex(kind(0d0)):: comega2

    b2( 1:4*nlayer )=0
    comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
    do  i=1,4*nlayer-3,4
        b2(i) = - dcmplx( p2(i) ) / comega2
    enddo
    do  i=2,4*nlayer-2,4
        b2(i) = - dcmplx( p2(i) ) / comega2
    enddo
    do  i=4,4*nlayer,4
        b2(i) = - dcmplx( p2(i) ) / comega2
    enddo
    !
    return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calhml( nlayer,coef1,coef2, h3mx,h3my,h3mz,h5mx,h5my,h5mz,&
    h4m1L,h4m1N,h4m2L,h4m2N,h6m1L,h6m1N,h6m2L,h6m2N,a1 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Summing up the modified first derivatives matrix operators.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    double precision:: h3mx(-2:1,*),h3my(-2:1,*),h3mz(-2:1,*)
    double precision:: h5mx(-1:2,*),h5my(-1:2,*),h5mz(-1:2,*)
    double precision:: h4m1L(-1:2,*),h4m1N(-1:2,*),h4m2L(-2:1,*),h4m2N(-2:1,*)
    double precision:: h6m1L(-1:2,*),h6m1N(-1:2,*),h6m2L(-2:1,*),h6m2N(-2:1,*)
    complex(kind(0d0)):: coef1,coef2,a1(4,*)
    integer:: i,j,m,n
    !
    do  j=2,2*(nlayer+1),2
        do 100 i=1,3,2
            if ( (i==1).and.(j==2) ) goto 100
            m = (-i+3)/2
            n = (i+j-3)/2
            a1(i,j) = a1(i,j)&
                - ( - coef2 * dcmplx( h3my(m,n) )& !-F
                + coef1 * dcmplx( 2.d0/3.d0 * h4m2N(m,n) )&
                + coef1 * dcmplx( h6m2L(m,n) ) )!L
100     continue
        enddo
        do   j=3,2*nlayer+1,2
            do 120 i=1,3,2
                if ( (i==1).and.(j==3) ) goto 120
                m = (-i+5)/2
                n = (i+j-4)/2
                a1(i,j) = a1(i,j)&
                    - (   coef1 * dcmplx( h4m1L(m,n) ) &!L
                    - coef2 * dcmplx( h5my(m,n) ) &!-F
                    + coef1 * dcmplx( 2.d0/3.d0 * h6m1N(m,n) ) )
120         continue
            enddo
            !
            return
        end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine overlapa(nlayer,a,a2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Overlapping the coefficient matrix elements in the solid part.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nlayer
    complex(kind(0d0)):: a(*),a2(4,*)
    integer:: i,nsize,msize,n1,n2
!
    nsize = 2 * ( nlayer + 1 )
    msize = 16 * nlayer
!
    a2(1:4,1:nsize)=0

    a2(4,1) = a(1)
    a2(3,2)  = a(2)
    a2(4,2)  = a(4)
    a2(4,nsize-1) = a(msize-3)
    a2(3,nsize)   = a(msize-2)
    a2(4,nsize)   = a(msize)
    do i=3,nsize-1,2
        n1 = 16 * (i-3)/2 + 5
        a2(2,i) = a(n1)
        n1 = n1 + 2
        a2(3,i) = a(n1)
        if ( i/=nsize-1 ) then
            n1 = n1 + 6
            n2 = n1 + 4
            a2(4,i) = a(n1) + a(n2)
        endif
        enddo
        do i=4,nsize,2
            n1 = 16 * (i-4)/2 + 6
            a2(1,i) = a(n1)
            n1 = n1 + 2
            a2(2,i) = a(n1)
            if ( i/=nsize ) then
                n1 = n1 + 6
                n2 = n1 + 4
                a2(3,i) = a(n1) + a(n2)
                n1 = n1 + 2
                n2 = n1 + 4
                a2(4,i) = a(n1) + a(n2)
            endif
        enddo
!
        return
    end
        !
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine overlapb(nlayer,b,b2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Overlapping the coefficient matrix elements in the liquid part.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: nlayer
    complex(kind(0d0))::b(*),b2(4,*)
    integer:: i,n1,n2

    b2( 1:4,1:nlayer+1 )=0

    b2(4,1) = b(1)
    do i=2,nlayer
        n1 = 4 * ( i - 1 )
        n2 = n1 + 1
        b2(4,i) = b(n1) + b(n2)
    enddo
    b2(4,nlayer+1) = b(4*nlayer)

    do i=2,nlayer+1
        n1 = 4 * ( i - 2 ) + 2
        b2(3,i) = b(n1)
    enddo
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine cala( maxnzone,ndc,iphase,nlayer,kkdr,kdr,ksp,&
        l2,lsq,nn,a0,a1,a2,a )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: maxnzone,ndc,iphase(*)
    integer:: nlayer(maxnzone)
    integer:: kkdr(*),kdr,ksp(maxnzone),nn
    double precision:: l2,lsq
    complex(kind(0d0)):: a0(4,*),a1(4,*),a2(4,*),a(4,*)
    integer:: j,isl,ill,ii,jj,mtmp,i1,i2

    ! --- Initializing the matrix elements
    !call cmatinit( 4,nn,a )
    a(1:4,1:nn)=0

    isl = 0
    ill = 0
    do j=1,ndc+1
        if ( iphase(j)==1 ) then
            isl = isl + 1
            i1 = kkdr(j)
            i2 = kkdr(j)+2*nlayer(j)+1
            do jj=i1,i2
                mtmp = kdr+ksp(j)+jj-kkdr(j)
                do  ii=2,4,2
                    a(ii,jj) = a(ii,jj) +  dcmplx(  dble(a0(ii,mtmp)) + l2 *  dble(a2(ii,mtmp)),&
                            dimag(a0(ii,mtmp)) + l2 * dimag(a2(ii,mtmp))  )
                enddo
                do ii=1,3,2
                    a(ii,jj) = a(ii,jj)  + dcmplx( lsq * dble(a1(ii,mtmp)),&
                    lsq * dimag(a1(ii,mtmp)) )
                enddo
            enddo
        else
            ill = ill + 1
            i1 = kkdr(j)
            i2 = kkdr(j)+nlayer(j)
            do jj=i1,i2
                mtmp = kdr+ksp(j)+jj-kkdr(j)
                do ii=3,4
                    a(ii,jj) = a(ii,jj)+dcmplx( dble(a0(ii,mtmp)) + l2 * dble(a2(ii,mtmp)),&
                    dimag(a0(ii,mtmp)) + l2 * dimag(a2(ii,mtmp)) )
                enddo
            enddo
        endif
    enddo
!
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calbc( maxnzone,ndc,rdc,iphase,kkdr,a )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the boundary condition elements.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: maxnzone,ndc,iphase(*),kkdr(*)
    double precision:: rdc(*)
    complex(kind(0d0)):: a(4,*)
    integer:: i

    do i=1,ndc
        if ( (iphase(i)==1).and.(iphase(i+1)==2) )  a( 2,kkdr(i+1) ) =   dcmplx( rdc(i) * rdc(i) )
        if ( (iphase(i)==2).and.(iphase(i+1)==1) )  a( 3,kkdr(i+1) ) = - dcmplx( rdc(i) * rdc(i) )
    enddo
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calabnum( omega,omegai,rmax,rrho, vpv,vph,vsv,vsh,eta,&
    ra,r0,coef1,coef2,anum,bnum )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    double precision:: omega,omegai,rmax
    double precision:: rrho(4),vpv(4),vph(4),vsv(4),vsh(4),eta(4)
    double precision:: ra(2),r0
    complex(kind(0d0)):: coef1,coef2,anum(4,4,10),bnum(4,4,10)
    integer:: i,j,k
    double precision:: trho,tmu,tecA,tecC,tecL,tecN,tAkappa,tCkappa
    double precision:: tvpv,tvph,tvsv,tvsh,teta,coef,r,r2
    complex(kind(0d0)):: xrho,xecA,xecC,xecF,xecL,xecN
    complex(kind(0d0)):: xmu,xFdC,xAmF2dC
    complex(kind(0d0)):: xAkappa,xCkappa,comega2
    !
    do k=1,10
        do  j=1,4
            do  i=1,4
                anum(i,j,k) = dcmplx( 0.d0 )
                bnum(i,j,k) = dcmplx( 0.d0 )
            enddo
        enddo
    enddo            ! computation of mu,lam and rho at r
    do i=1,10
        if ( i<=5 ) then
            r = ra(1) + dble(i-1) / 4.d0 * ( r0 - ra(1) )
        else
            r = ra(2) - dble(10-i) / 4.d0 * ( ra(2) - r0 )
        endif
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
            coef = coef * ( r / rmax )
            endif
            trho  = trho  + rrho(j)  * coef
            tvpv  = tvpv  + vpv(j)   * coef
            tvph  = tvph  + vph(j)   * coef
            tvsv  = tvsv  + vsv(j)   * coef
            tvsh  = tvsh  + vsh(j)   * coef
            teta  = teta  + eta(j)   * coef
        enddo
        tecA = trho * tvph * tvph
        tecC = trho * tvpv * tvpv
        tecL = trho * tvsv * tvsv
        tecN = trho * tvsh * tvsh
        tmu  = trho * tvsv * tvsv
!
        tAkappa = tecA - tmu * 4.d0 / 3.d0
        tCkappa = tecC - tmu * 4.d0 / 3.d0
        xAkappa = dcmplx(tAkappa) * coef2
        xCkappa = dcmplx(tCkappa) * coef2
        xrho = dcmplx(trho)
        xmu  = dcmplx(tmu) *coef1
        xecL = dcmplx(tecL) * coef1
        xecN = dcmplx(tecN) * coef1
        xecA = xAkappa + xmu * dcmplx( 4.d0 / 3.d0 )
        xecC = xCkappa + xmu * dcmplx( 4.d0 / 3.d0 )
        xecF = dcmplx(teta) * ( xecA - dcmplx(2.d0) * xecL )
        xFdC = xecF / xecC
        xAmF2dC = xecA - xecF * xFdC
        !
        r2 = r * r
        comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
                !
        anum(1,1,i) = - dcmplx( 2.d0 / r ) * xFdC
        anum(1,2,i) =   dcmplx( 1.d0 ) / xecC
        anum(2,1,i) = - xrho * comega2   + dcmplx( 4.d0 / r2 ) * ( xAmF2dC - xecN )
        anum(2,2,i) =   dcmplx( 2.d0 / r ) * ( xFdC - 1.d0 )
        anum(3,1,i) = - dcmplx( 1.d0 / r )
        anum(3,3,i) =   dcmplx( 1.d0 / r )
        anum(3,4,i) =   dcmplx( 1.d0 ) / xecL
        anum(4,1,i) = - dcmplx( 2.d0 / r2 ) * ( xAmF2dC - xecN )
        anum(4,2,i) = - dcmplx( 1.d0 / r ) * xFdC
        anum(4,3,i) = - xrho * comega2 - dcmplx( 2.d0 / r2 ) * xecN
        anum(4,4,i) = - dcmplx( 3.d0 / r )
        bnum(1,3,i) =   dcmplx( 1.d0 / r ) * xFdC
        bnum(2,3,i) = - dcmplx( 2.d0 / r2 ) * ( xAmF2dC - xecN )
        bnum(2,4,i) =   dcmplx( 1.d0 / r )
        bnum(4,3,i) =   xAmF2dC / dcmplx(r2)
    enddo
    !
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calya( aa,bb,l2,ra,r0,ya,yb,yc,yd )
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the excitation vector using Geller and Hatori(1995).
!                                1995.7     N.Takeuchi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    common a,b,cl2,itmp
    external eqmotion1
! ---------------------------<< constants >>---------------------------
    double precision,parameter ::pi=3.1415926535897932d0
! ---------------------------<< variables >>---------------------------
! variables for the source
    double precision:: r0
! variables for the numerical integration
    integer:: i,k,itmp
    double precision:: l2,cl2
    complex(kind(0d0)):: ya(4),yb(4),yc(4),yd(4),yn(4)
    double precision:: xs,xe,dr,ra(2)
    complex(kind(0d0)):: work(4,2),aa(160),bb(160),a(64),b(64)
    integer ktmp1,ktmp2,ktmp3,ktmp4

    ! -----------------------<< common variables >>-----------------------
    cl2 = l2
    ! ---------------------<< numerical integration >>---------------------
    ! integration from the lower boundary
    ya=0
    ya(2) = dcmplx( 1.d0 )
    yb=0
    yb(4) = dcmplx( 1.d0 )
    dr = ( r0 - ra(1) ) / 2.d0
    if ( dr>0.d0 ) then
        xs = ra(1)
        do  k=1,2
            ktmp1 = 32 * ( k-1 )
            ktmp2 = ktmp1 - 16
            do i=1,64
                if ( i<=32 ) then
                    a(i) = aa(i+ktmp1)
                    b(i) = bb(i+ktmp1)
                else
                    a(i) = aa(i+ktmp2)
                    b(i) = bb(i+ktmp2)
                endif
            enddo
            itmp = 0
            xe = xs + dr
            call rk3(4,xs,xe,1,ya,yn,4,work)
            itmp = 0
            xs = xe - dr
            call rk3(4,xs,xe,1,yb,yn,4,work)
        enddo
    endif
! integration from the upper boundary
    yc = 0
    yc(2) = dcmplx( 1.d0 )
    yd = 0
    yd(4) = dcmplx( 1.d0 )
    dr = ( ra(2) - r0 ) / 2.d0
    if ( dr>0.d0 ) then
        xs = ra(2)
        do  k=1,2
            ktmp1 = 144 - 32 * ( k-1 )
            ktmp2 = ktmp1 - 32
            ktmp3 = ktmp1 - 48
            ktmp4 = ktmp1 - 80
            do  i=1,64
                if ( i<=32 ) then
                    if ( i<=16 ) then
                        a(i) = aa(i+ktmp1)
                        b(i) = bb(i+ktmp1)
                    else
                        a(i) = aa(i+ktmp2)
                        b(i) = bb(i+ktmp2)
                    endif
                else
                    if ( i<=48 ) then
                        a(i) = aa(i+ktmp3)
                        b(i) = bb(i+ktmp3)
                    else
                        a(i) = aa(i+ktmp4)
                        b(i) = bb(i+ktmp4)
                    endif
                endif
            enddo
            itmp = 0
            xe = xs - dr
            call rk3(4,xs,xe,1,yc,yn,4,work)
            itmp = 0
            xs = xe + dr
            call rk3(4,xs,xe,1,yd,yn,4,work)
        enddo
    endif
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calg( l,m,coef1,coef2,lsq,ecC0,ecF0,ecL0,&
          ya,yb,yc,yd,ra,r0,mt,g )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computing the excitation vector using Geller and Hatori(1995).
!                                1995.7     N.Takeuchi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
! ---------------------------<< constants >>---------------------------
    double precision:: pi
    parameter ( pi=3.1415926535897932d0 )
! ---------------------------<< variables >>---------------------------
! variables for the structure
    complex(kind(0d0)):: coef1,coef2
! variables for the source
    double precision:: r0,mt(3,3),ecC0,ecF0,ecL0
    complex(kind(0d0)):: dd,ee,s1,s2,s(4)
! variables for the numerical integration
    integer l,m
    complex(kind(0d0)):: ya(4),yb(4),yc(4),yd(4)
    complex(kind(0d0)):: g(4)
    double precision:: ra(2)
    ! other variables
    integer ip(4),ier,i
    complex(kind(0d0)):: a(4,4),b(4),wk(4),xtmp
    double precision:: eps,sgn,b1,b2,lsq,r03,dtmp(4)
    data eps / -1.d0 /

! ---------------------<< parameter computation >>---------------------
! computation of the discontinuity
    if ( m>=0 ) then
        sgn = 1.d0
    else
        sgn = - 1.d0
    endif
    b1 = dsqrt( dble(2*l+1)/(16.d0*pi) )
    if ( l/=0 ) then
        b2 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2)/(64.d0*pi) )
    endif
    if ( iabs(m)==2 ) then
        dd = dcmplx( 0.d0 )
        ee = dcmplx( 0.d0 )
    endif
    if ( iabs(m)==1 ) then
        dd = dcmplx( 0.d0 )
        ee = dcmplx( - b1 * sgn * mt(1,2), b1 * mt(1,3) )    / ( dcmplx( r0 * r0 * ecL0 ) * coef1 )
    endif
    if ( iabs(m)==0 ) then
        dd = dcmplx( 2.d0 * b1 * mt(1,1) / ( r0 * r0 ) ) &
            / (   dcmplx( ecC0 - 4.d0 / 3.d0 * ecL0 ) * coef2&
            + dcmplx( 4.d0/3.d0*ecL0 ) * coef1 )
        ee = dcmplx( 0.d0 )
    endif
!
    if ( iabs(m)==0 ) then
        r03 = r0 * r0 * r0
        xtmp =   ( ( ecF0 + 2.d0 / 3.d0 * ecL0) * coef2 &
            - 2.d0 / 3.d0 * ecL0 * coef1 )&
            / ( (ecC0 - 4.d0 / 3.d0 * ecL0)* coef2 &
            + 4.d0 / 3.d0 * ecL0 * coef1 )
        s1 = dcmplx( - 2.d0 * b1 * ( mt(2,2) + mt(3,3) ) / r03 ) + dcmplx( 4.d0 * b1 * mt(1,1) / r03 ) * xtmp
        s2 = dcmplx( b1 * lsq * ( mt(2,2) + mt(3,3) ) / r03 ) - dcmplx( 2.d0 * b1 * lsq * mt(1,1) / r03 ) * xtmp
    endif
    if ( iabs(m)==1 ) then
        s1 = dcmplx( 0.d0 )
        s2 = dcmplx( 0.d0 )
    endif
    if ( iabs(m)==2 ) then
        r03 = r0 * r0 * r0
        s1 = dcmplx( 0.d0 )
        s2 = dcmplx( - b2 * ( mt(2,2) - mt(3,3) ) / r03 ,&
            sgn * 2.d0 * b2 * mt(2,3) / r03 )
    endif
    s(1) = dd
    s(2) = s1
    if ( l/=0 ) then
        s(3) = dcmplx( dble(ee)/lsq, dimag(ee)/lsq )
        s(4) = dcmplx( dble(s2)/lsq, dimag(s2)/lsq )
    else
        s(3) = dcmplx( 0.d0 )
        s(4) = dcmplx( 0.d0 )
    endif
! consideration of the boundary conditions
! determination of the analytical solution
    if ( l/=0 ) then
        call sab1( ya,yb,yc,yd,s,a,b )
        call glu(a,4,4,b,eps,wk,ip,ier)
    else
        call sab2( ya,yc,s,a,b )
        call glu(a,2,4,b,eps,wk,ip,ier)
        b(3) = b(2)
        b(2) = dcmplx( 0.d0 )
        b(4) = dcmplx( 0.d0 )
    endif
! computation of the excitation vector
    dtmp(1) = - ra(1) * ra(1)
    dtmp(2) = dtmp(1) * lsq
    dtmp(3) = ra(2) * ra(2)
    dtmp(4) = dtmp(3) * lsq
    do i=1,4
        xtmp = b(i)
        g(i) = dcmplx( dble(xtmp)*dtmp(i), dimag(xtmp)*dtmp(i) )
    enddo
!
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine eqmotion1(r,y,f)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Equation of the motion of the solid medium.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    common a,b,l2,itmp
    double precision:: r,l2
    complex(kind(0d0)):: y(4),f(4)
    integer i,j,itmp,mtmp
    complex(kind(0d0)):: a(64),b(64),c(4,4)
!
! computation of itmp
    itmp = itmp + 1
! computation of the differential coefficients
    mtmp = 16 * ( itmp - 1 )
    do j=1,4
        do i=1,4
            mtmp = mtmp + 1
            c(i,j) = a(mtmp)   + dcmplx( l2*dble(b(mtmp)),  l2*dimag(b(mtmp)) )
        enddo
    enddo
    do i=1,4
        f(i) = dcmplx( 0.d0 )
        do  j=1,4
            f(i) = f(i) + c(i,j) * y(j)
        enddo
    enddo
!
    return
    end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine sab1( ya,yb,yc,yd,s,a,b )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! determination of the matrix imposing the boundary conditions.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    complex(kind(0d0)):: ya(4),yb(4),yc(4),yd(4),s(4),a(4,4),b(4)
    integer i
!
    do i=1,4
        a(i,1) = - ya(i)
        a(i,2) = - yb(i)
        a(i,3) =   yc(i)
        a(i,4) =   yd(i)
        b(i) = s(i)
    enddo
!
    return
    end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine sab2( ya,yc,s,a,b )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! determination of the matrix imposing the boundary conditions.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    complex(kind(0d0)):: ya(4),yc(4),s(4),a(4,4),b(4)
    integer i
!
    do i=1,2
        a(i,1) = - ya(i)
        a(i,2) =   yc(i)
        b(i) = s(i)
    enddo
    return
    end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine rea2( nn,a,g,c,d,nzone,iphase,kkdr,spn,kkdr0,nn0 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! rearranging the matrix elements for l=0.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nn,nzone
    integer:: iphase(*),kkdr(*),spn,kkdr0,nn0
    complex(kind(0d0)):: a(4,*),c(2,*),g(*),d(*)
    integer:: izone,i1,i2,itmp,mtmp,i
!
    itmp = 0
    do izone=1,nzone
        if ( iphase(izone)==1 ) then
            if ( izone==spn ) kkdr0 = itmp + 1
            i1 = itmp + 1
            if ( izone/=nzone ) then
                i2 = i1 + ( kkdr(izone+1) - kkdr(izone) ) / 2 - 1
            else
                i2 = i1 + ( nn + 1 - kkdr(izone) ) / 2 - 1
            endif
            d(i1) = g(kkdr(izone))
            if ( izone/=1 ) then
                if ( iphase(izone-1)==2 ) then
                    c(1,i1) = a(3,kkdr(izone))
                else
                    c(1,i1) = a(2,kkdr(izone))
                endif
            endif
            c(2,i1) = a(4,kkdr(izone))
            do i=i1+1,i2
                mtmp = kkdr(izone) + 2 * ( i - i1 )
                d(i) = g(mtmp)
                c(1,i) = a(2,mtmp)
                c(2,i) = a(4,mtmp)
            enddo
            itmp = i2
        else
            i1 = itmp + 1
            i2 = i1 + ( kkdr(izone+1) - kkdr(izone) ) - 1
            d(i1) = g(kkdr(izone))
            if ( izone/=1 ) then
                if ( iphase(izone-1)==1 ) then
                    c(1,i1) = a(2,kkdr(izone))
                else
                    c(1,i1) = a(3,kkdr(izone))
                endif
            endif
            c(2,i1) = a(4,kkdr(izone))
            do i=i1+1,i2
                mtmp = kkdr(izone) + ( i - i1 )
                d(i) = g(mtmp)
                c(1,i) = a(3,mtmp)
                c(2,i) = a(4,mtmp)
            enddo
            itmp = i2
        endif
    enddo
    nn0 = itmp

    return
    end
