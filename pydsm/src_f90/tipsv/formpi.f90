!cccccccccc       these subroutines are for mpi
!ccc
subroutine simplesplit(istart, iend, n, imin, imax)
    !c     this routine returns the imin(i) and imax(i) (i=1,...,n)
    !c      separates istart--iend into n parts
    !c     each part contains imin(i) -- imax(i) i= 1, ..., n
    !c
    !c     iend-istart+1 = inum = n * deltai + remainder
    !c     remainder is included in i th irange
    !c     iend = n*deltai + remainder + istart -1
    !c     istart = iend-remainder-n*deltai+1
    !c     n * deltai = iend-(start+remainder) +1
    !c
    !c     remain :: istart, istart+1, ..., (istart + remainder-1)
    !c     the others ::
    !c     i = (istart+iamari), ...., imax ------------- (deltai*n)
    !c     i(i) = istart+remainder+(i-1)*deltai,
    !c     ......istart+remainder+n*deltai-1
    !cccccccccccccccccccccccccccccccccccccc
    !c     input parameters
    implicit none
    integer, intent(in) :: istart, iend
    integer, intent(in) :: n  ! the number of processors
    !c
    integer, dimension(*), intent(out) :: imin, imax
    !c
    integer :: remainder
    integer :: inum
    integer :: deltai
    !c
    integer :: i
    !c
    inum = iend-istart+1
    remainder = mod (inum, n)
    deltai = (inum-remainder)/n
    !c
    imin (1) = istart
    imax (1) = istart + remainder-1 + deltai
    do i = 2, n
        imin(i)= istart + remainder + (i-1) * deltai
        imax(i)= istart + remainder + i * deltai -1
    enddo
    return
end
      
subroutine trianglesplit (istart, iend, n, imin, imax)
    !c     this routine returns the imin(i) and imax(i) (i=1,...,n)
    !c      separates istart--iend into n parts
    !c     each part contains imin(i) -- imax(i) i= 1, ..., n
    !c
    !c     iend-istart+1 = inum
    !c
    !c     Assume that the cpu time t for i th omega is a * i
    !c     then we divide a* 0.5 * i **2 into n parts.
    !c
    !c     return imin imax which satisfy above assumption
    !c     input parameters
    implicit none
    integer, intent(in) :: istart, iend
    integer, intent(in) :: n  ! the number of processors
    !c
    integer, dimension(*), intent(out) :: imin, imax
    !c
    integer :: remainder
    integer :: inum
    integer :: deltai
    !c
    integer :: i
    integer:: x(0:n)
    double precision :: s  !  0.5 *  iend **2 / n
    double precision ::p
    !c
    inum = iend-istart+1
    s =  iend *iend / n
    !c
    x(0) = istart
    do i = 1, n
        x(i)= s+x(i-1)**2
        x(i) = x(i)**0.5
    enddo

    do i = 1, n
        imin(i)= x(i-1)
        imax(i)= x(i)-1
    enddo
    imax (n) = iend
    return

end
      
      
