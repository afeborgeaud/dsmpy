program tipsv
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ************** mpi-tipsv.f ****************
! Computation of PSV synthetic seismograms
! in transversely isotropic for anisotropic PREM
! using modified DSM operators & modified source representation.
! Synthetics for shallow events can be computed.
!                                                 2002.12  K.Kawai
!       2009. ?      Kensuke Konishi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
! ---------------------------<< constants >>---------------------------
double precision pi,lmaxdivf,shallowdepth
integer maxnlay,maxnslay,maxnllay
integer maxnzone,maxnr,maxlmax,ilog
parameter ( pi=3.1415926535897932d0 )
parameter ( maxnlay = 80880 )
parameter ( maxnslay = 48840 )
parameter ( maxnllay = 32040 )
parameter ( maxnzone = 20 )
parameter ( maxnr = 500 )
parameter ( maxlmax = 80000 )
parameter ( lmaxdivf = 2.d4 )
parameter ( shallowdepth = 100.d0 )
! ---------------------------<< variables >>---------------------------
! variable for the trial function
integer nnlayer,nlayer(maxnzone)
integer nslay,nllay
integer inlayer,jnlayer,jnslay,jnllay
integer l,m
double precision ra(maxnlay+maxnzone+1),plm(3,0:3,maxnr)
complex(kind(0d0)) bvec(3,-2:2,maxnr)
! variable for the structure
integer nzone,isl,ill,nsl,nll
integer iphase(maxnzone),ndc,vnp
double precision rmin,rmax
double precision vrmin(maxnzone),vrmax(maxnzone)
double precision rrho(4,maxnzone)
double precision vpv(4,maxnzone),vph(4,maxnzone)
double precision vsv(4,maxnzone),vsh(4,maxnzone),eta(4,maxnzone)
double precision qmu(maxnzone),qkappa(maxnzone)
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
complex(kind(0d0)) coef1(maxnzone),coef2(maxnzone),coef(maxnzone)
! variable for the periodic range
integer np,imin,imax
double precision tlen,omega,omegai
complex(kind(0d0)) u(3,maxnr)
! variable for the source
integer spn,ns
double precision r0,mt(3,3),spo,eqlat,eqlon
double precision ecC0,ecF0,ecL0
complex(kind(0d0)) ya(4),yb(4),yc(4),yd(4)
! variable for the station
integer nr,ir
double precision theta(maxnr),phi(maxnr)
double precision lat(maxnr),lon(maxnr)
! variable for the matrix elements
complex(kind(0d0)) a0(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
complex(kind(0d0)) a1(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
complex(kind(0d0)) a2(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
complex(kind(0d0)) a( 4, 2*(maxnslay+1) + (maxnllay+1) )
complex(kind(0d0)) c( 2, (maxnslay+1) + (maxnllay+1) )
double precision  t( 8*maxnslay )
double precision h1x( 8*maxnslay ),h1y( 8*maxnslay ),h1z( 8*maxnslay )
double precision h2L( 8*maxnslay ),h2N( 8*maxnslay )
double precision h3ax( 8*maxnslay )
double precision h3ay( 8*maxnslay ),h3az( 8*maxnslay )
double precision h4aL( 8*maxnslay ),h4aN( 8*maxnslay )
double precision h5ax( 8*maxnslay ),h5ay( 8*maxnslay ),h5az( 8*maxnslay )
double precision h6aL( 8*maxnslay ),h6aN( 8*maxnslay )
double precision h3x( 8*maxnslay ),h3y( 8*maxnslay ),h3z( 8*maxnslay )
double precision h4L( 8*maxnslay ),h4N( 8*maxnslay )
double precision h5x( 8*maxnslay ),h5y( 8*maxnslay ),h5z( 8*maxnslay )
double precision h6L( 8*maxnslay ),h6N( 8*maxnslay )
double precision h7x( 8*maxnslay ),h7y( 8*maxnslay ),h7z( 8*maxnslay )
double precision h8L( 8*maxnslay ),h8N( 8*maxnslay )
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
complex(kind(0d0)) g(2*(maxnslay+1) + (maxnllay+1) )
complex(kind(0d0)) d((maxnslay+1) + (maxnllay+1) )
! variable for the file
character*80 output(maxnr)
! variable for the stack point
integer isp(maxnzone)
integer issp(maxnzone)
integer ilsp(maxnzone),jssp(maxnzone)
integer jsp(maxnzone)
integer ksp(maxnzone),lsp(maxnzone)
integer isdr,jsdr,ildr
! variables for the gridding
integer jjdr(maxnzone),kkdr(maxnzone)
integer jdr,kdr
double precision vmin(maxnzone),gridpar(maxnzone),dzpar(maxnzone)
! variables for l cut off
complex(kind(0d0)) tmpc(2*(maxnslay+1) + (maxnllay+1))
integer sufzone,ismall,kc,lsuf,llog
double precision maxamp,ratc,ratl,re
! variables for the numerical integration
complex(kind(0d0)) anum(4,4,10),bnum(4,4,10)
! variables for tapering
double precision dlmax0,ctaper
! other variables
integer i,j,ii,jj,nn,lda,ier,itmp,jtmp,mtmp,kkdr0,nn0
integer ll(12),lli(12),llj(12)
double precision eps,work(8*maxnslay),l2,lsq
complex(kind(0d0)) z( 2*(maxnslay+1) + (maxnllay+1) )
complex(kind(0d0)) w( 2*(maxnslay+1) + (maxnllay+1) )
complex(kind(0d0)) cwork( 2*(16*maxnslay + 4*maxnllay) )
integer ltmp(2),iimax

data lda/ 4 /
data eps/ -1.d0 /
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     MPI
!ccccccccccccccccccccccccc
include 'mpif.h'
integer :: petot, my_rank, ista, ierr
complex(kind(0d0)), allocatable, dimension(:,:,:) :: outputu
!     when the values to be output use memory over outputmemory MB,&
!     they are written in output files. The interval is outputinterval
!     memoryperomega is the quantity of memory used for one omega step
integer, allocatable, dimension (:) :: mpimin, mpimax
call mpi_init (ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
!     write(*,*) "myrank",my_rank
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
allocate(mpimin(PETOT), mpimax(PETOT))
!ccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! *************** Inputting and computing the parameters ***************
! --- inputting parameter ---
if (my_rank ==0) then
    call pinput2(maxnzone,maxnr,re,ratc,ratl,tlen,np,omegai,imin,imax,&
        nzone,vrmin,vrmax,rrho,vpv,vph,vsv,vsh,eta,qmu,qkappa,&
        r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)
endif

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!      write(*,*) my_rank
!      call mpi_finalize(ierr)
!      stop
!cc  nnlayer vpv vph eta qkappa
call MPI_BCAST (nnlayer, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
call MPI_BCAST (qkappa, maxnzone, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)

call MPI_BCAST (vpv, 4*maxnzone, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (vph, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (eta, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)

!cc
call MPI_BCAST (re, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
call MPI_BCAST (ratc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
call MPI_BCAST (ratl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
call MPI_BCAST (tlen, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (np, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (omegai, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (imin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (imax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (nzone, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (vrmin, maxnzone, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (vrmax, maxnzone, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (rrho, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (vsv, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (vsh, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (qmu, maxnzone, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (r0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
call MPI_BCAST (eqlat, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (eqlon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (mt, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
call MPI_BCAST (nr, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (theta, maxnr, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (phi, maxnr, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (lat(1), maxnr, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (lon(1), maxnr, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
call MPI_BCAST (output, 80*maxnr, MPI_CHARACTER,0, MPI_COMM_WORLD, ierr)
!cccccccccccccccccccccccccccccccccccccccccccc
allocate (outputu(3,nr,0:np))

! --- computing the required parameters ---
! counting of the nsl and nll
call calnl( nzone,vsv,iphase,nsl,nll )
! computing and checking the parameters
rmin = vrmin(1)
rmax = vrmax(nzone)
ndc = nzone - 1

    theta= theta / 180 * pi
    phi= phi   / 180 * pi

if ( (r0<rmin) .or. (r0>rmax) ) stop 'Location of the source is improper.'

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
iimax = imax
!ccccccccccccccccccccccccccccccccccccccccccccccccc

if( (rmax-r0)<shallowdepth) then ! option for shallow events
    ! computing of the number and the location of grid points
    iimax = int(tlen * 2.d0)
    call calgrid( nzone,vrmin,vrmax,vpv,vsv,rmin,rmax,&
               iimax,1,tlen,vmin,gridpar,dzpar )
    call calra( maxnlay,maxnslay,maxnllay,maxnzone,&
         nnlayer,inlayer,jnlayer,jnslay,jnllay,&
       gridpar,dzpar,nzone,vrmin,vrmax,iphase,&
       rmin,rmax,r0,nslay,nllay,nlayer,ra,re )
    ! --- checking the parameter
    if ( inlayer>maxnlay ) stop 'The number of grid points is too large.'
    if ( nslay>maxnslay ) stop 'The number of the grid points in the solid is too large.'
    if ( nllay>maxnllay ) stop 'The number of the grid points in the liquid is too large.'
    if ( jnlayer>2*maxnlay ) stop 'The number of the total grid points is too large.'
    if ( ( jnslay>2*maxnslay ).or.( jnllay>2*maxnllay ) ) stop 'The number of the total grid points is too large.'
    ! computing the stack points
    call calsp( maxnzone,ndc,nsl,nll,iphase,nlayer,nslay,nllay,&
        isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr,ildr,jdr,kdr )
    ! computing the source location
    call calspo( maxnlay,maxnzone,ndc,vrmax,iphase,inlayer,r0,rmin,rmax,ra,isp,spo,spn )
    ! ******************* Computing the matrix elements *******************
    ! data initialization
    a=0
    t=0
    h1x=0
    h1y=0
    h1z=0
    h2L=0
    h2N=0
    h3ax=0
    h3ay=0
    h3az=0
    h4aL=0
    h4aN=0
    h5ax=0
    h5ay=0
    h5az=0
    h6aL=0
    h6aN=0
    h3x=0
    h3y=0
    h3z=0
    h4L=0
    h4N=0
    h5x=0
    h5y=0
    h5z=0
    h6L=0
    h6N=0
    h7x=0
    h7y=0
    h7z=0
    h8L=0
    h8N=0
    h3mx=0
    h3my=0
    h3mz=0
    h4m2L=0
    h4m2N=0
    h6m2L=0
    h6m2N=0
    h5mx=0
    h5my=0
    h5mz=0
    h4m1L=0
    h4m1N=0
    h6m1L=0
    h6m1N=0
    p1=0
    p2=0
    p3=0

    ! computing the structure grid points
    call calstg( maxnlay,maxnzone,nzone,iphase,rrho,vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,&
        vnp,vra,rho,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN,r0,spn,ecC0,ecF0,ecL0 )
    rhoinv(1:vnp)   = 1.d0 / rho(1:vnp)
    kappainv(1:vnp) = 1.d0 / kappa(1:vnp)
    isl = 0
    ill = 0
    do i=1,ndc+1
        if ( iphase(i)==1 ) then
            isl = isl + 1
            itmp = isdr+issp(isl)
            call calmatc( nlayer(i),vnp,vra,rho,2,0,0,ra(isp(i)), t(itmp) )
            call caltl( nlayer(i),vnp,vra,rho,ra(isp(i)),work(itmp) )
            call calt( nlayer(i), t(itmp),work(itmp), t(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKx,0,0,0,ra(isp(i)),h1x(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKx,ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h1x(itmp),work(itmp),h1x(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKy,0,0,0,ra(isp(i)),h1y(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKy,ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h1y(itmp),work(itmp),h1y(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecKz,0,0,0,ra(isp(i)),h1z(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKz,ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h1z(itmp),work(itmp),h1z(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecL,0,0,0,ra(isp(i)),h2L(itmp) )
            call calhl( nlayer(i),vnp,vra,ecL,ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h2L(itmp),work(itmp),h2L(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h2N(itmp) )
            call calhl( nlayer(i),vnp,vra,ecN,ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h2N(itmp),work(itmp),h2N(itmp) )
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
            call calmatc( nlayer(i),vnp,vra,ecL ,2,1,1,ra(isp(i)), h8L(itmp) )
            call calmatc( nlayer(i),vnp,vra,ecN ,2,1,1,ra(isp(i)), h8N(itmp) )
        else
            ill = ill + 1
            itmp = ildr+ilsp(ill)
            call calmatc( nlayer(i),vnp,vra,rhoinv,2,1,1,ra(isp(i)),p1(itmp) )
            call calmatc( nlayer(i),vnp,vra,rhoinv,0,0,0,ra(isp(i)),p2(itmp) )
            call calhl( nlayer(i),vnp,vra,rhoinv,ra(isp(i)),work(itmp) )
            call calt( nlayer(i),p2(itmp),work(itmp),p2(itmp) )
            call calmatc( nlayer(i),vnp,vra,kappainv,2,0,0,ra(isp(i)),p3(itmp) )
            call caltl( nlayer(i),vnp,vra,kappainv,ra(isp(i)),work(itmp) )
            call calt( nlayer(i),p3(itmp),work(itmp),p3(itmp) )
        endif
    enddo
    ! Computing the modified operator of the 1st derivative
    call caltstg( maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,&
        nlayer,ra,rmax,vra,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN )
    isl = 0
    do i=1,ndc+1
        if ( iphase(i)==1 ) then
            isl = isl + 1
            itmp = isdr+issp(isl)
            jtmp = isp(i)+i-1
            call calh5( nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h5ax(itmp),work(itmp),h5x(itmp) )
            call calh5( nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h5ay(itmp),work(itmp),h5y(itmp) )
            call calh5( nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h5az(itmp),work(itmp),h5z(itmp) )
            call calh5( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h6aL(itmp),work(itmp),h6L(itmp) )
            call calh5( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h6aN(itmp),work(itmp),h6N(itmp) )
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
    !
    do ii=1,2              ! omega-loop
        !
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
            call calabnum( omega,omegai,rmax,&
                rrho(1,spn),vpv(1,spn),vph(1,spn),&
                vsv(1,spn),vsh(1,spn),eta(1,spn),&
                ra(mtmp),r0,coef1(spn),coef2(spn),&
                anum(1,1,1),bnum(1,1,1) )
            ! computing the matrix elements independent of l
            isl = 0
            ill = 0
            do j=1,ndc+1
                if ( iphase(j)==1 ) then
                    isl = isl + 1
                    itmp = isdr+issp(isl)
                    jtmp = jdr+jsp(j)
                    mtmp = kdr+ksp(j)
                    call cala0( nlayer(j),omega,omegai,&
                        t(itmp), h1x(itmp), h1y(itmp),h1z(itmp),&
                                            h2L(itmp), h2N(itmp),&
                                            h3ax(itmp),h3ay(itmp),h3az(itmp),&
                                            h4aL(itmp),h4aN(itmp),&
                                            h5ax(itmp),h5ay(itmp),h5az(itmp),&
                                            h6aL(itmp),h6aN(itmp),&
                                            h7x(itmp),h7y(itmp),h7z(itmp),&
                                            h8L(itmp), h8N(itmp),&
                                            coef1(j),coef2(j),cwork(jtmp) )
                    call overlapa( nlayer(j),cwork(jtmp),a0(1,mtmp) )
                    call cala1( nlayer(j),&
                                            h1x(itmp),h1y(itmp),h1z(itmp),&
                                            h2L(itmp),h2N(itmp),&
                        h3x(itmp), h3y(itmp), h3z(itmp), h4L(itmp), h4N(itmp), &
                        h5x(itmp), h5y(itmp), h5z(itmp),&
                                            h6L(itmp), h6N(itmp),&
                                            coef1(j),coef2(j),cwork(jtmp) )
                    call overlapa( nlayer(j),cwork(jtmp),a1(1,mtmp) )
                    call cala2( nlayer(j), h1x(itmp), h1y(itmp),h1z(itmp),&
                                            h2L(itmp),h2N(itmp),&
                        coef1(j),coef2(j),cwork(jtmp) )
                    call overlapa( nlayer(j),cwork(jtmp),a2(1,mtmp) )
                    jtmp = jsdr+jssp(isl)
                    call calhml( nlayer(j),coef1(j),coef2(j),&
                             h3mx(-2,jtmp),h3my(-2,jtmp),h3mz(-2,jtmp),&
                             h5mx(-1,jtmp),h5my(-1,jtmp),h5mz(-1,jtmp),&
                             h4m1L(-1,jtmp),h4m1N(-1,jtmp),&
                             h4m2L(-2,jtmp),h4m2N(-2,jtmp),&
                             h6m1L(-1,jtmp),h6m1N(-1,jtmp),&
                             h6m2L(-2,jtmp),h6m2N(-2,jtmp),&
                             a1(1,mtmp) )
                else
                    ill = ill + 1
                    itmp = ildr+ilsp(ill)
                    jtmp = jdr+jsp(j)
                    mtmp = kdr+ksp(j)
                    call calb0( nlayer(j),omega,omegai,p1(itmp),p3(itmp),coef(j),cwork(jtmp) )
                    call overlapb( nlayer(j),cwork(jtmp),a0(1,mtmp) )
                    call calb2( nlayer(j),omega,omegai,p2(itmp),coef(j),cwork(jtmp) )
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
                ! computing the coefficient matrix elements
                ! --- renewing  mdr
                if (  mod(l,50)==0  ) then
                    call calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,rmax,sufzone )
                    call calspdr( maxnzone,nzone,iphase,nlayer,jjdr,kkdr )
                    nn = kkdr(nzone) + 2 * nlayer(nzone) + 1
                endif
                ! computing the matrix elements
                call cala( maxnzone,ndc,iphase,nlayer,kkdr,kdr,ksp,l2,lsq,nn,a0,a1,a2,a )
                ! computing the boundary condition elements
                call calbc( maxnzone,ndc,vrmax,iphase,kkdr,a )

                jtmp = kkdr(spn) + 2 * int(spo)
                mtmp = isp(spn) + int(spo)
                if ( spo==int(spo) ) then
                    jtmp = jtmp - 2
                    mtmp = mtmp - 1
                endif
                call calya( anum(1,1,1),bnum(1,1,1),l2,ra(mtmp),r0,ya,yb,yc,yd )
                do m=-2,2     ! m-loop
                    if ( iabs(m)<=iabs(l) ) then
                        ! computing the excitation vector
                        g( 1:nn )=0
                        call calg( l,m,coef1(spn),coef2(spn),lsq,&
                                 ecC0,ecF0,ecL0,ya,yb,yc,yd,&
                                 ra(mtmp),r0,mt,g(jtmp) )
                        if ( l==0 ) then
                            ! rearranging the matrix elements for l=0
                            call rea2( nn,a,g,c,d,nzone,iphase,kkdr,spn,kkdr0,nn0 )
                            ! computing the expansion coefficient
                            itmp=1
                            if ( rmin==0.d0 ) itmp=2
                            ns = kkdr0 + ( nint(spo) - 1 )
                            call dcsymbdl( c(1,itmp),1,nn0-itmp+1,1,eps,z(itmp),w(itmp),ll,lli,llj,ier )
                            call dcsbdlv( c(1,itmp),d(itmp),1,nn0-itmp+1,ns-itmp+1,eps,z(itmp),ier )
                        else
                            ! computing the expansion coefficient
                            itmp=1
                            if ( rmin==0.d0 ) itmp=3
                            ns = kkdr(spn) + 2 * ( nint(spo) - 1 )
                            if( mod(l,100)==0) then
                                if ( ( m==-2 ).or.( m==-l ) )&
                                 call dcsymbdl0( a(1,itmp),3,nn-itmp+1,6,eps,z(itmp),w(itmp),ll,lli,llj,ier )

                                call dcsbdlv0( a(1,itmp),g(itmp),3,nn-itmp+1,eps,z(itmp),ier )
                                ! sum up c of the same l
                                tmpc(1:nn) = tmpc(1:nn) + g(1:nn)

                            else
                                if ( kc<3 ) kc = 3
                                itmp = kc
                                if ( ( m==-2 ).or.( m==-l ) ) then
                                    call dcsymbdl( a(1,itmp),3,nn-itmp+1,6,eps, z(itmp),w(itmp),ll,lli,llj,ier )
                                endif
                                call dcsbdlv( a(1,itmp),g(itmp),3,nn-itmp+1,ns-itmp+1,eps,z(itmp),ier )
                            endif
                            ! computing ampratio
                            call calamp( g(nn-1),l,lsuf,maxamp,ismall,ratl )

                            if( mod(l,100)==0) &
                             call calcutd(nzone,nlayer,tmpc,ratc,nn,iphase,spo,spn,ra,kkdr,kc)

                        ! computing the displacement
                        endif
                    endif
                enddo         ! m-loop
            enddo            ! l-loop
        endif
    enddo                  ! omega-loop
    iimax = dble(max(ltmp(1),ltmp(2))) * tlen / lmaxdivf
endif                     ! option for shallow events
call calgrid( nzone,vrmin,vrmax,vpv,vsv,rmin,rmax,&
         iimax,1,tlen,vmin,gridpar,dzpar )
call calra( maxnlay,maxnslay,maxnllay,maxnzone,nnlayer,inlayer,jnlayer,jnslay,jnllay,&
         gridpar,dzpar,nzone,vrmin,vrmax,iphase,rmin,rmax,r0,nslay,nllay,nlayer,ra,re )
! --- checking the parameter
if ( inlayer>maxnlay ) stop 'The number of grid points is too large.'
if ( nslay>maxnslay ) stop 'The number of the grid points in the solid is too large.'
if ( nllay>maxnllay ) stop 'The number of the grid points in the liquid is too large.'
if ( jnlayer>2*maxnlay ) stop 'The number of the total grid points is too large.'
if ( ( jnslay>2*maxnslay ).or.( jnllay>2*maxnllay ) ) stop 'The number of the total grid points is too large.'
! computing the stack points
call calsp( maxnzone,ndc,nsl,nll,iphase,nlayer,nslay,nllay,&
    isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr,ildr,jdr,kdr )
! computing the source location
call calspo( maxnlay,maxnzone,ndc,vrmax,iphase,inlayer,r0,rmin,rmax,ra,isp,spo,spn )
! ******************* Computing the matrix elements *******************
! data initialization
    a=0
    t=0
    h1x=0
    h1y=0
    h1z=0
    h2L=0
    h2N=0
    h3ax=0
    h3ay=0
    h3az=0
    h4aL=0
    h4aN=0
    h5ax=0
    h5ay=0
    h5az=0
    h6aL=0
    h6aN=0
    h3x=0
    h3y=0
    h3z=0
    h4L=0
    h4N=0
    h5x=0
    h5y=0
    h5z=0
    h6L=0
    h6N=0
    h7x=0
    h7y=0
    h7z=0
    h8L=0
    h8N=0
    h3mx=0
    h3my=0
    h3mz=0
    h4m2L=0
    h4m2N=0
    h6m2L=0
    h6m2N=0
    h5mx=0
    h5my=0
    h5mz=0
    h4m1L=0
    h4m1N=0
    h6m1L=0
    h6m1N=0
    p1=0
    p2=0
    p3=0

! computing the structure grid points
call calstg( maxnlay,maxnzone,nzone,iphase,rrho,&
         vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,&
         vnp,vra,rho,kappa,ecKx,ecKy,ecKz,&
         mu,ecL,ecN,r0,spn,ecC0,ecF0,ecL0 )
rhoinv(1:vnp)   = 1.d0 / rho(1:vnp)
kappainv(1:vnp) = 1.d0 / kappa(1:vnp)
isl = 0
ill = 0
do i=1,ndc+1
    if ( iphase(i)==1 ) then
        isl = isl + 1
        itmp = isdr+issp(isl)
        call calmatc( nlayer(i),vnp,vra,rho ,2,0,0,ra(isp(i)), t(itmp) )
        call caltl( nlayer(i),vnp,vra,rho,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), t(itmp),work(itmp), t(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKx,0,0,0,ra(isp(i)),h1x(itmp) )
        call calhl( nlayer(i),vnp,vra,ecKx,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h1x(itmp),work(itmp),h1x(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKy,0,0,0,ra(isp(i)),h1y(itmp) )
        call calhl( nlayer(i),vnp,vra,ecKy,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h1y(itmp),work(itmp),h1y(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKz,0,0,0,ra(isp(i)),h1z(itmp) )
        call calhl( nlayer(i),vnp,vra,ecKz,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h1z(itmp),work(itmp),h1z(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecL ,0,0,0,ra(isp(i)),h2L(itmp) )
        call calhl( nlayer(i),vnp,vra,ecL,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h2L(itmp),work(itmp),h2L(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h2N(itmp) )
        call calhl( nlayer(i),vnp,vra,ecN,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h2N(itmp),work(itmp),h2N(itmp) )
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
        call calt( nlayer(i),p2(itmp),work(itmp),p2(itmp) )
        call calmatc( nlayer(i),vnp,vra,kappainv,2,0,0,ra(isp(i)),p3(itmp) )
        call caltl( nlayer(i),vnp,vra,kappainv,ra(isp(i)),work(itmp) )
        call calt( nlayer(i),p3(itmp),work(itmp),p3(itmp) )
    endif
enddo
! Computing the modified operator of the 1st derivative
call caltstg( maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,&
         nlayer,ra,rmax,vra,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN )
isl = 0
do i=1,ndc+1
    if ( iphase(i)==1 ) then
        isl = isl + 1
        itmp = isdr+issp(isl)
        jtmp = isp(i)+i-1
        call calh5( nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h5ax(itmp),work(itmp),h5x(itmp) )
        call calh5( nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h5ay(itmp),work(itmp),h5y(itmp) )
        call calh5( nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h5az(itmp),work(itmp),h5z(itmp) )
        call calh5( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h6aL(itmp),work(itmp),h6L(itmp) )
        call calh5( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h6aN(itmp),work(itmp),h6N(itmp) )
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
! ******************** Computing the displacement *********************
!
llog = 0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      call simplesplit (imin, imax, PETOT, mpimin, mpimax)
call trianglesplit (imin, imax, PETOT, mpimin, mpimax)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
do i=mpimin(my_rank+1), mpimax(my_rank+1)            ! omega-loop
    u=0
    if ( i/=0 ) then
        omega = 2.d0 * pi * dble(i) / tlen
        call callsuf(omega,nzone,vrmax,vsv,lsuf)
        plm(1:3,0:3,1:nr)=0
        call calcoef( nzone,omega,qmu,qkappa,coef1,coef2,coef )
        mtmp = isp(spn) + int(spo)
        if ( spo==int(spo) ) mtmp = mtmp - 1
        call calabnum( omega,omegai,rmax,&
                       rrho(1,spn),vpv(1,spn),vph(1,spn),&
                       vsv(1,spn),vsh(1,spn),eta(1,spn),&
                       ra(mtmp),r0,coef1(spn),coef2(spn),&
                       anum(1,1,1),bnum(1,1,1) )
        ! computing the matrix elements independent of l
        isl = 0
        ill = 0
        do j=1,ndc+1
            if ( iphase(j)==1 ) then
                isl = isl + 1
                itmp = isdr+issp(isl)
                jtmp = jdr+jsp(j)
                mtmp = kdr+ksp(j)
                call cala0( nlayer(j),omega,omegai,&
                    t(itmp),h1x(itmp), h1y(itmp),h1z(itmp),&
                    h2L(itmp), h2N(itmp),&
                    h3ax(itmp),h3ay(itmp),h3az(itmp),&
                    h4aL(itmp),h4aN(itmp),&
                    h5ax(itmp),h5ay(itmp),h5az(itmp),&
                    h6aL(itmp),h6aN(itmp),&
                    h7x(itmp),h7y(itmp),h7z(itmp),&
                    h8L(itmp), h8N(itmp),&
                    coef1(j),coef2(j),cwork(jtmp) )
                call overlapa( nlayer(j),cwork(jtmp),a0(1,mtmp) )
                call cala1( nlayer(j),&
                    h1x(itmp),h1y(itmp),h1z(itmp),&
                    h2L(itmp),h2N(itmp),&
                    h3x(itmp), h3y(itmp),h3z(itmp),h4L(itmp),h4N(itmp), &
                    h5x(itmp), h5y(itmp), h5z(itmp), h6L(itmp), h6N(itmp),&
                    coef1(j),coef2(j),cwork(jtmp) )
                call overlapa( nlayer(j),cwork(jtmp),a1(1,mtmp) )
                call cala2( nlayer(j),h1x(itmp), h1y(itmp),h1z(itmp),&
                    h2L(itmp),h2N(itmp),&
                    coef1(j),coef2(j),cwork(jtmp) )
                call overlapa( nlayer(j),cwork(jtmp),a2(1,mtmp) )
                jtmp = jsdr+jssp(isl)
                call calhml( nlayer(j),coef1(j),coef2(j),&
                    h3mx(-2,jtmp),h3my(-2,jtmp),h3mz(-2,jtmp),&
                    h5mx(-1,jtmp),h5my(-1,jtmp),h5mz(-1,jtmp),&
                    h4m1L(-1,jtmp),h4m1N(-1,jtmp),&
                    h4m2L(-2,jtmp),h4m2N(-2,jtmp),&
                    h6m1L(-1,jtmp),h6m1N(-1,jtmp),&
                    h6m2L(-2,jtmp),h6m2N(-2,jtmp), a1(1,mtmp) )
            else
                ill = ill + 1
                itmp = ildr+ilsp(ill)
                jtmp = jdr+jsp(j)
                mtmp = kdr+ksp(j)
                call calb0( nlayer(j),omega,omegai,p1(itmp),p3(itmp),coef(j),cwork(jtmp) )
                call overlapb( nlayer(j),cwork(jtmp),a0(1,mtmp) )
                call calb2( nlayer(j),omega,omegai,p2(itmp),coef(j),cwork(jtmp) )
                call overlapb( nlayer(j),cwork(jtmp),a2(1,mtmp) )
            endif
        enddo
        !
        maxamp = -1.0
        ismall = 0
        llog = maxlmax
        kc = 1
        do l=0,maxlmax      ! l-loop
   ! initialize
                tmpc(1:maxnlay+1) = dcmplx(0.d0)
            if ( ismall>20 ) then
                if(llog>l) llog = l
                exit
            endif
            l2 = dble(l)*dble(l+1)
            lsq = dsqrt( l2 )
            ! computing the coefficient matrix elements
            ! --- renewing  mdr
            if (  mod(l,50)==0  ) then
                call calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,rmax,sufzone )
                call calspdr( maxnzone,nzone,iphase,nlayer,jjdr,kkdr )
                nn = kkdr(nzone) + 2 * nlayer(nzone) + 1
            endif
            ! ***** Computing the trial function *****
            do ir=1,nr
                call calbvec( l,theta(ir),phi(ir),plm(1,0,ir),bvec(1,-2,ir) )
            enddo
            ! computing the matrix elements
            call cala( maxnzone,ndc,iphase,nlayer,kkdr,kdr,ksp,l2,lsq,nn,a0,a1,a2,a )
            ! computing the boundary condition elements
            call calbc( maxnzone,ndc,vrmax,iphase,kkdr,a )
            !
            jtmp = kkdr(spn) + 2 * int(spo)
            mtmp = isp(spn) + int(spo)
            if ( spo==int(spo) ) then
                jtmp = jtmp - 2
                mtmp = mtmp - 1
            endif
            call calya( anum(1,1,1),bnum(1,1,1),l2,ra(mtmp),r0,ya,yb,yc,yd )
            do m=-2,2        ! m-loop
                if ( iabs(m)<=iabs(l) ) then
                    ! computing the excitation vector
                    g( 1:nn )=0
                    call calg( l,m,coef1(spn),coef2(spn),lsq,&
                             ecC0,ecF0,ecL0,ya,yb,yc,yd,ra(mtmp),r0,mt,g(jtmp) )
                    if ( l==0 ) then
                        ! rearranging the matrix elements for l=0
                        call rea2( nn,a,g,c,d,nzone,iphase,kkdr,spn,kkdr0,nn0 )
                        ! computing the expansion coefficient
                        itmp=1
                        if ( rmin==0.d0 ) itmp=2
                        ns = kkdr0 + ( nint(spo) - 1 )
                        call dcsymbdl( c(1,itmp),1,nn0-itmp+1,1,eps,z(itmp),w(itmp),ll,lli,llj,ier )
                        call dcsbdlv( c(1,itmp),d(itmp),1,nn0-itmp+1,ns-itmp+1,eps,z(itmp),ier )
                        ! computing the displacement
                        do ir=1,nr
                            call calu0( d(nn0),bvec(1,m,ir),u(1,ir) )
                        enddo
                    else
                        ! computing the expansion coefficient
                        itmp=1
                        if ( rmin==0.d0 ) itmp=3
                        ns = kkdr(spn) + 2 * ( nint(spo) - 1 )
                        if( mod(l,100)==0) then
                            if ( ( m==-2 ).or.( m==-l ) ) &
                             call dcsymbdl0( a(1,itmp),3,nn-itmp+1,6,eps, z(itmp),w(itmp),ll,lli,llj,ier)

                            call dcsbdlv0( a(1,itmp),g(itmp),3,nn-itmp+1,eps,z(itmp),ier )
                            ! sum up c of the same l
                            tmpc(1:nn) = tmpc(1:nn) + g(1:nn)
                        else
                            if ( kc<3 ) kc = 3
                            itmp = kc
                            if ( ( m==-2 ).or.( m==-l ) )&
                             call dcsymbdl( a(1,itmp),3,nn-itmp+1,6,eps, z(itmp),w(itmp),ll,lli,llj,ier )

                            call dcsbdlv( a(1,itmp),g(itmp),3,nn-itmp+1,ns-itmp+1,eps,z(itmp),ier )
                        endif
                        ! computing ampratio
                        call calamp( g(nn-1),l,lsuf,maxamp,ismall,ratl )
                        !
                        if( mod(l,100)==0)&
                         call calcutd(nzone,nlayer,tmpc,ratc,nn,iphase,spo,spn,ra,kkdr,kc)

                        ! computing the displacement
                        do ir=1,nr
                            call calu( g(nn-1),lsq,bvec(1,m,ir),u(1,ir) )
                        enddo
                    endif
                endif
            enddo            ! m-loop
        enddo               ! l-loop
    endif
    ! ************************** Files Handling **************************

    outputu(1:3,1:nr,i)=u(1:3,1:nr)
enddo                   ! omega-loop

if(my_rank /=0) then
    call mpi_send(outputu(1,1,mpimin(my_rank+1)),&
        3*nr*(mpimax(my_rank+1)-mpimin(my_rank+1)+1),&
        MPI_DOUBLE_complex,0,my_rank, MPI_COMM_WORLD, ierr)
endif

if(my_rank==0)then
    do i = 1,petot-1
        call mpi_recv(outputu(1,1,mpimin(i+1)),&
            3*nr*(mpimax(i+1)-mpimin(i+1)+1),&
            MPI_DOUBLE_complex,i,i,MPI_COMM_WORLD, ista,ierr)
    enddo
    ! ************************** Files Handling **************************
    do ir = 1 ,nr
        open(unit=10,file=trim(output(ir)),&
            status='unknown', form='unformatted',&
           access='stream', convert='big_endian')
        write(10) tlen,np,1,3,omegai,lat(ir),lon(ir)
        write(10) eqlat,eqlon,r0
        do i= imin, imax
            write(10) i,dble(outputu(1,ir,i)),dimag(outputu(1,ir,i))
            write(10) dble(outputu(2,ir,i)),dimag(outputu(2,ir,i))
            write(10) dble(outputu(3,ir,i)),dimag(outputu(3,ir,i))
        enddo
        close(10)
    enddo
endif

!	call mpi_finalize(ierr)
!	stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
write(*,*) my_rank, "Ivalice looks to the horizon"
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
call mpi_finalize (ierr)
stop
end

