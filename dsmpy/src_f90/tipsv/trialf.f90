!----------------------------------------------------------
!                           TRIAL FUNCTION
!----------------------------------------------------------
subroutine calbvec( l,theta,phi,plm,bvec )
!----------------------------------------------------------
! Evaluating the value of toroidal harmonics (fully normalized)
! at each station whose latitude and longitude are theta and phi.
!----------------------------------------------------------
    use parameters
    implicit none
    integer,intent(in):: l
    integer:: m,i
    double precision:: theta,phi,x,plm(3,0:3),fact,coef
    complex(dp):: bvec(3,-2:2),expimp

    x = dcos( theta )
    do m=0,min0(l,3)
        call calplm( l,m,x,plm(1,m) )
    enddo
    do m=0,min0(l,2)
        fact = 1
        if ( m/=0 ) then
            do i=l-m+1,l+m
                fact = fact * dble(i)
            enddo
        endif
        coef = dsqrt( (2*l+1)/(4.d0*pi) / fact )
        expimp = cdexp( dcmplx( 0, m*phi ) )
        bvec(1,m)  = coef * plm(1,m) * expimp
        bvec(1,-m) = dconjg( bvec(1,m) )
        bvec(2,m) = coef *( m*x / dsin( theta ) * plm(1,m) + plm(1,m+1) ) * expimp
        bvec(2,-m) = dconjg( bvec(2,m) )
        bvec(3,m)  = dcmplx( 0, m ) / dsin( theta ) * coef * plm(1,m) * expimp
        bvec(3,-m) = dconjg( bvec(3,m) )
        if ( mod(m,2)==1 ) bvec(1:3,-m) = - bvec(1:3,-m)
    enddo
    return
end

!----------------------------------------------------------
subroutine calplm( l,m,x,plm )
!----------------------------------------------------------
    implicit none
    integer:: l,m,i
    double precision:: x,plm(3),pmm,somx2,fact

    if ( m<0 .or. m>l .or. dabs(x)>1.d0 ) stop 'bad arguments'
    if ( l==m ) then
        pmm = 1
        if ( m>0 ) then
            somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
            fact = 1
            do i=1,m
                pmm = -pmm * fact * somx2
                fact = fact + 2
            enddo
        endif
        plm(2:3) = 0
        plm(1) = pmm
    else
        plm(3) = plm(2)
        plm(2) = plm(1)
        if ( l==m+1 ) then
            plm(1) = x * (2*m+1) * plm(2)
        else
            plm(1) = ( x * (2*l-1) * plm(2) - (l+m-1) * plm(3) ) / (l-m)
        endif
    endif
    return
end



