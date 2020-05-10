!***********************************************************************
!                           TRIAL FUNCTION                             *
!***********************************************************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calbvec( l,theta,phi,plm,bvec )
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! Evaluating the value of toroidal harmonics (fully normalized)
    ! at each station whose latitude and longitude are theta and phi.
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    double precision:: pi
    parameter ( pi=3.1415926535897932d0 )
    !
    integer:: l,m,i
    double precision:: theta,phi,x,plm(3,0:3),fact,coef
    complex(kind(0d0)):: bvec(3,-2:2),expimp
    !
    x = dcos( theta )
    do  m=0,min0(l,3)
        call calplm( l,m,x,plm(1,m) )
    enddo
    do  m=0,min0(l,2)
        fact = 1.d0
        if ( m/=0 ) then
            do  i=l-m+1,l+m
                fact = fact * dble(i)
            enddo
        endif
        coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
        expimp = cdexp( dcmplx( 0.d0, dble(m)*phi ) )
        bvec(1,m)  = coef * plm(1,m) * expimp
        bvec(1,-m) = dconjg( bvec(1,m) )
        bvec(2,m) = coef * ( dble(m) * x / dsin( theta ) * plm(1,m)&
            + plm(1,m+1) )  * expimp
        bvec(2,-m) = dconjg( bvec(2,m) )
        bvec(3,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta )&
            * coef * plm(1,m) * expimp
        bvec(3,-m) = dconjg( bvec(3,m) )
        if ( mod(m,2)==1 ) then
            bvec(1,-m) = - bvec(1,-m)
            bvec(2,-m) = - bvec(2,-m)
            bvec(3,-m) = - bvec(3,-m)
        endif
    enddo
    !
    return
end
!

subroutine caldveczero( l,bvec )
    implicit none
    double precision, parameter ::  pi=3.1415926535897932d0
    integer:: l,m,i
    double precision :: fact,coef
    complex(kind(0d0)) :: bvec(1:3,-2:2)
    double precision :: xl2

    bvec = cmplx(0.d0)
    xl2 = dble(l) * dble(l+1)

    do m=0,min0(l,1)
        fact = 1.d0
        if ( m/=0 ) then
            do i=l-m+1,l+m
                fact = fact * dble(i)
            enddo
        endif
        coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )

        if(m==0) then
            bvec(1,m) = cmplx(coef)
        endif
        if(m==1) then
            bvec(2,m) = dcmplx(dble(m),0.d0) * xl2 * coef / 2.d0
            bvec(2,-m) = conjg(bvec(2,m))
            bvec(3,m) = dcmplx( 0.d0, dble(m)) *  xl2 *coef / 2.d0
            bvec(3,-m) = conjg(bvec(3,m))
        endif

        if(mod(m,2)==1) then
            bvec(2,-m) = -bvec(2,-m)
            bvec(3,-m) = -bvec(3,-m)
        endif
    enddo
    return
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calplm( l,m,x,plm )
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: l,m,i
    double precision:: x,plm(3),pmm,somx2,fact
    !
    if ( m<0.or.m>l.or.dabs(x)>1.d0 ) stop 'bad arguments'
    if ( l==m ) then
        pmm = 1.d0
        if ( m>0 ) then
            somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
            fact = 1.d0
            do i=1,m
                pmm = -pmm * fact * somx2
                fact = fact + 2.d0
            enddo
        endif
        plm(3) = 0.d0
        plm(2) = 0.d0
        plm(1) = pmm
    else
        plm(3) = plm(2)
        plm(2) = plm(1)
        if ( l==m+1 ) then
            plm(1) = x * dble(2*m+1) * plm(2)
        else
            plm(1) = ( x * dble(2*l-1) * plm(2)  - dble(l+m-1) * plm(3)  ) / dble(l-m)
        endif
    endif
    !
    return
end


