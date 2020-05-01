module parameters
    integer, parameter :: dp = kind(1.d0)
    real(dp), parameter :: pi=3.1415926535897932d0

    integer, parameter :: maxnzone = 15
    integer, parameter :: maxnr = 600
    integer, parameter :: maxlmax = 80000
    integer, parameter :: maxnlay = 88300
    integer, parameter :: ilog = 0
    real(dp), parameter :: lmaxdivf = 2.d4
    real(dp), parameter :: shallowdepth = 100.d0
    

    type input_params
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
        character*80 :: output(maxnr)
    end type input_params

contains

    subroutine input_params_cpy(params, &
        re,ratc,ratl,tlen,np,omegai,imin,imax, &
        nzone,vrmin,vrmax,rho,vsv,vsh,qmu, &
        r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)
  
        implicit none
  
        type(input_params), intent(in) :: params
        integer, intent(out) :: np
        integer, intent(out) :: imin,imax
        integer, intent(out) :: nzone,nr
        real*8, intent(out) :: tlen,omegai,re,ratc,ratl
        real*8, dimension(maxnzone), intent(out) :: vrmin,vrmax,qmu
        real*8, dimension(4,maxnzone), intent(out) :: rho,vsv,vsh
        real*8, intent(out) :: r0,mt(3,3)
        real*8, dimension(maxnr), intent(out) :: theta,phi,lat,lon
        real*8, intent(out) :: eqlat,eqlon
        character*80, intent(out) :: output(maxnr)
  
        re = params%re
        np = params%np
        imin = params%imin
        imax = params%imax
        nzone = params%nzone
        nr = params%nr
        tlen = params%tlen
        omegai = params%omegai
        re = params%re
        ratc = params%ratc
        ratl = params%ratl
        vrmin = params%vrmin
        vrmax = params%vrmax
        qmu = params%qmu
        rho = params%rho
        vsv = params%vsv
        vsh = params%vsh
        r0 = params%r0
        mt = params%mt
        theta = params%theta
        phi = params%phi
        lat = params%lat
        lon = params%lon
        eqlat = params%eqlat
        eqlon = params%eqlon
        output = params%output
  
        return
    end

end module parameters