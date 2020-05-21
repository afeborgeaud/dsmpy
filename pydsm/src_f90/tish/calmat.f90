! calmat.f for wcalprem.f
!------------------------------------------------------------------------
	subroutine calmatc( nlayer,vnp,vra, &
     	                    con,rpow,w1dn,w2dn,ra,m,work )
!------------------------------------------------------------------------
! Computing \int r^rpow con W_p^(w1dn) W_q^(w2dn) dr.
!------------------------------------------------------------------------
	use parameters

	integer, parameter :: maxrpow = 2
	integer, intent(in) :: nlayer,vnp,rpow,w1dn,w2dn
	real*8, intent(in) :: vra(vnp),con(vnp)
	real(dp), dimension(maxnlay+maxnzone+1), intent(in) :: ra
	real*8, dimension(4*nlayer), intent(out) :: m
	real*8, dimension(4*nlayer) :: work
	integer i,j,k,kk,l,nn,snp
	real*8 a(2,2),b(2,2),c(5),rh
! parameter check
	if ( rpow .gt. maxrpow ) then
		write(*,*) 'Invalid arguments.(calmatc)'
		return
	endif
! computing of the matrix elements
	snp = 1
	do 140 i=1,nlayer
	  rh = ra(i+1) - ra(i)
	  if ( w1dn .eq. 0 ) then
	    a(2,1) = -1.d0 / rh
	    a(1,1) = ra(i+1) / rh
	    a(2,2) = 1.d0 / rh
	    a(1,2) = -ra(i) / rh
	  else
	    if ( w1dn .eq. 1 ) then
	      a(2,1) = 0.d0
	      a(1,1) = -1.d0 / rh
	      a(2,2) = 0.d0
	      a(1,2) = 1.d0 / rh
		else
		  write(*,*) 'Invalid arguments.(calmatc)'
		  return
	    endif
	  endif
	  if ( w2dn .eq. 0 ) then
	    b(2,1) = -1.d0 / rh
	    b(1,1) = ra(i+1) / rh
	    b(2,2) = 1.d0 / rh
	    b(1,2) = -ra(i) / rh
	  else
	    if ( w2dn .eq. 1 ) then
	      b(2,1) = 0.d0
	      b(1,1) = -1.d0 / rh
	      b(2,2) = 0.d0
	      b(1,2) = 1.d0 / rh
		else
			write(*,*) 'Invalid arguments.(calmatc)'
			return
	    endif
	  endif
	  do 130 j=1,2
	    do 120 k=1,2
	      do 100 kk=1,5
	        c(kk) = 0.d0
  100	      continue
	      call pmulti( 2,a(1,j),2,b(1,k),3,c )
	      do 110 l=3,1,-1
	        c(l+rpow) = c(l)
	        if ( rpow.gt.0 ) c(l)=0.d0
  110	      continue
	      nn = 4 * (i-1) + 2 * (j-1) + k
	      call pinteg( snp,5,c,ra(i),ra(i+1), &
     	                   vnp,vra,con,work(nn) )
  120	    continue
  130	  continue
  140	continue
	if ( w1dn .ne. w2dn ) then
	  do 150 i=1,4*nlayer
	    if ( (mod(i,4).eq.0) .or. (mod(i,4).eq.1) ) then
	      m(i) = 2.d0 * work(i)
	    else
	      if ( mod(i,4).eq.2 ) then
	        m(i) = work(i) + work(i+1)
	      else
	        m(i) = work(i-1) + work(i)
	      endif
	    endif
  150	  continue
	else
	  do 160 i=1,4*nlayer
	    m(i) = work(i)
  160	  continue
	endif
!
	return
	end
!
!------------------------------------------------------------------------
	subroutine caltl( nlayer,vnp,vra,rho,ra,tl)
!------------------------------------------------------------------------
! Computing of lumped mass matrix.
!------------------------------------------------------------------------
	integer nlayer,vnp
	real*8 vra(vnp),rho(vnp),ra(*),tl(*)
	integer i,nn,snp
	real*8 c(3),from,to
!
	snp = 1
	do 100 i=1,nlayer
	  c(1) = 0.d0
	  c(2) = 0.d0
	  c(3) = 1.d0
!
	  from = ra(i)
	  to = ( ra(i) + ra(i+1) ) / 2.d0
	  nn = 4 * (i-1)
	  call pinteg( snp,3,c,from,to,vnp,vra,rho,tl(nn+1) )
!
	  tl(nn+2) = 0.d0
	  tl(nn+3) = 0.d0
!
	  from = to
	  to = ra(i+1)
	  call pinteg( snp,3,c,from,to,vnp,vra,rho,tl(nn+4) )
!
  100	continue
!
	return
	end
!
!------------------------------------------------------------------------
	subroutine calhl( nlayer,vnp,vra,mu,ra,hl)
!------------------------------------------------------------------------
! Computing of lumped rigidity matrix.
!------------------------------------------------------------------------
	integer nlayer,vnp
	real*8 vra(vnp),mu(vnp),ra(*),hl(*)
	integer i,nn,snp
	real*8 c(1),from,to
!
	snp = 1
	do 100 i=1,nlayer
	  c(1) = 1.d0
!
	  from = ra(i)
	  to = ( ra(i) + ra(i+1) ) / 2.d0
	  nn = 4 * (i-1)
	  call pinteg( snp,1,c,from,to,vnp,vra,mu,hl(nn+1) )
!
	  hl(nn+2) = 0.d0
	  hl(nn+3) = 0.d0
!
	  from = to
	  to = ra(i+1)
	  call pinteg( snp,1,c,from,to,vnp,vra,mu,hl(nn+4) )
!
  100	continue
!
	return
	end
!
!------------------------------------------------------------------------
	subroutine calt( nlayer, tl, tc, t )
!------------------------------------------------------------------------
	integer nlayer, i
	real*8 t(*), tl(*), tc(*)
!
	do 100 i=1,4*nlayer
	  t(i) = ( tl(i) + tc(i) ) / 2.d0
  100	continue
!
	return
	end
!
!------------------------------------------------------------------------
	subroutine pmulti(n,a,m,b,l,c)
!------------------------------------------------------------------------
! Computing the (l-1) degrees polynomial c(n) which is the product of
! (n-1) degrees polynomial a(n) and (m-1) degrees polynomial b(n).
!------------------------------------------------------------------------
	integer n,m,l,i,j
	real*8 a(n),b(m),c(l)
!
	if (n+m-1.ne.l) then
		write(*,*) 'Invalid arguments.(pmulti)'
		return
	endif
	do 100 i=1,l
	  c(i)=0.d0
  100	continue
	do 210 i=1,n
	  do 200 j=1,m
	    c(i+j-1) = c(i+j-1) + a(i)*b(j)
  200	  continue
  210	continue
!
      return
      end
!
!------------------------------------------------------------------------
      subroutine pinteg(snp,n,p,from,to,vnp,vra,con,pint)
!------------------------------------------------------------------------
! Evaluating the integrated value pint from 'from' to 'to' of p(n)*con
! which is the product of (n-1) degrees polynomial 'p(n)' and 'con'.
!------------------------------------------------------------------------
	integer maxn
	parameter ( maxn = 5 )
!
	integer snp,n,vnp
	real*8 from,to,p(n),vra(vnp),con(vnp),pint
	real*8 x1,x2,q(2),pq(maxn+1),psint
!
	if ( n.gt.maxn ) then
		write(*,*) 'Degrees of a polynomial is too large.(pinteg)'
		return
	endif
	if ( snp.ge.vnp ) snp = 1
!
	pint = 0.d0
	x1 = from
  100	continue
	  if ( (vra(snp).le.x1).and.(vra(snp+1).gt.x1) ) then
	    x2 = dmin1( to, vra(snp+1) )
	  else
	    snp = snp + 1
	    if ( snp.eq.vnp ) snp = 1
	    goto 100
	  endif
! evaluating the integrated value
	q(2) = ( con(snp+1)-con(snp) ) / ( vra(snp+1)-vra(snp) )
	q(1) = - q(2) * vra(snp) + con(snp)
	call pmulti( n,p,2,q,n+1,pq )
	call polint( n+1,pq,x1,x2,psint )
	pint = pint + psint
!
	if ( x2.eq.to ) then
	  continue
	else
	  x1 = x2
	  goto 100
	endif
!
	return
	end
!
!------------------------------------------------------------------------
	subroutine polint( n,p,x1,x2,pint )
!------------------------------------------------------------------------
! Evaluating the integrated value 'pint' from 'x1' to 'x2'
! of (n-1) degrees polynomial 'p(n)'.
!------------------------------------------------------------------------
	integer maxn
	parameter ( maxn = 6 )
!
	integer n
	real*8 p(*),x1,x2,pint
	integer i,j
	real*8 a(maxn),b(maxn),dx,xx
!
	if ( n.gt.maxn ) then
		write(*,*) 'Degrees of a polynomial is too large.(polint)'
		return
	endif
!
	a(1) = 1.d0
	b(1) = 1.d0
	if ( n.ge.2 ) then
	  do 100 i=2,n
	    a(i) = a(i-1) * x1
	    b(i) = b(i-1) * x2
  100	  continue
	endif
	dx = x2 - x1
!
	pint = 0.d0
	do 120 i=1,n
	  xx = 0.d0
	  do 110 j=1,i
	    xx = xx + a(j) * b(i-j+1) / dble(i)
  110	  continue
	  pint = pint + p(i) * dx * xx
  120	continue
!
	return
	end
!
!------------------------------------------------------------------------
	subroutine cala0( nlayer,omega,omegai,t,h1,h2,h3,h4,coef,a0 )
!------------------------------------------------------------------------
! Computing the coefficient matrix 'a' in the solid part.
!------------------------------------------------------------------------
	integer nlayer
	real*8 omega,omegai
	real*8 t(*)
	real*8 h1(*),h2(*),h3(*),h4(*)
	complex*16 comega2,coef,a0(*)
	integer i
	real*8 h
!
	comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
	do 110 i=1,4*nlayer
	  h = h1(i) - h2(i) + h3(i) - 2.d0 * h4(i)
	  a0(i) = comega2 * dcmplx( t(i) ) - coef * dcmplx( h )
  110	continue
!
	return
	end
!
!------------------------------------------------------------------------
	subroutine cala2( nlayer,h4,coef,a2 )
!------------------------------------------------------------------------
! Computing the coefficient matrix 'a' in the solid part.
!------------------------------------------------------------------------
	integer nlayer
	real*8 h4(*)
	complex*16 coef,a2(*)
	integer i
!
	do 110 i=1,4*nlayer
	  a2(i) = - coef * dcmplx( h4(i) )
  110	continue
!
	return
	end
!
!------------------------------------------------------------------------
	subroutine cala( nn,l,lda,a0,a2,a )
!------------------------------------------------------------------------
! Computing the coefficient matrix 'a' in the solid part.
!------------------------------------------------------------------------
	integer nn,l,lda
	complex*16 a0(lda,*),a2(lda,*),a(lda,*)
	integer i,j
	real*8 xl2
!
	xl2 = dble(l) * dble(l+1)
	do 110 j=1,nn
	  do 100 i=1,2
	    a(i,j) = a0(i,j) + dcmplx(xl2) * a2(i,j)
  100	  continue
  110	continue
	return
	end
!
!------------------------------------------------------------------------
	subroutine calga(nlayer,omega,omegai,l,t,h1,h2,h3,h4,coef,a)
!------------------------------------------------------------------------
! Computing the coefficient matrix 'a' in the solid part.
!------------------------------------------------------------------------
	implicit none
	integer nlayer,l
	real*8 omega,omegai
	real*8 t(*)
	real*8 h1(*),h2(*),h3(*),h4(*)
	complex*16 comega2,coef,a(*)
	integer i
	real*8 h,xl2m2
!
	h = 0.d0
	do 100 i=1,4*nlayer
	   a(i) = 0.d0
 100	continue
	comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
	xl2m2 = dble(l) * dble(l+1) -2.d0
	do 110 i=1,4*nlayer
	  h = h1(i) - h2(i) + h3(i) + xl2m2 * h4(i)
	  a(i) = comega2 * dcmplx( t(i) ) - coef * dcmplx( h )
  110	continue
	return
	end
!------------------------------------------------------------------------
	subroutine calga2(nlayer,omega,omegai,l,t,h1,h2,h3,coef,a)
!------------------------------------------------------------------------
! Computing the coefficient matrix 'a' in the solid part.
!------------------------------------------------------------------------
	implicit none
	integer nlayer,l
	real*8 omega,omegai
	real*8 t(*)
	real*8 h1(*),h2(*),h3(*)
	complex*16 comega2,coef,a(*)
	integer i
	real*8 h,xl2m1
!
	h = 0.d0
	do 100 i=1,4*nlayer
	   a(i) = 0.d0
 100	continue
	comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
	xl2m1 = dble(l) * dble(l+1) -1.d0
	do 110 i=1,4*nlayer
	  h = h1(i) - h2(i) + xl2m1 * h3(i)
	  a(i) = comega2 * dcmplx( t(i) ) - coef * dcmplx( h )
  110	continue
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
