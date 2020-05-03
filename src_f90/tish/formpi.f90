cccccccccc       these subroutines are for mpi
ccc
      subroutine simplesplit(istart, iend, n, imin, imax,
     & omegaindex)
c     this routine returns the imin(i) and imax(i) (i=1,...,n)
c      separates istart--iend into n parts
c     each part contains imin(i) -- imax(i) i= 1, ..., n
c
c     iend-istart+1 = inum = n * deltai + remainder 
c     remainder is included in i th irange
c     iend = n*deltai + remainder + istart -1
c     istart = iend-remainder-n*deltai+1
c     n * deltai = iend-(start+remainder) +1
c     
c     remain :: istart, istart+1, ..., (istart + remainder-1)
c     the others ::
c     i = (istart+iamari), ...., imax ------------- (deltai*n)
c     i(i) = istart+remainder+(i-1)*deltai,
c     ......istart+remainder+n*deltai-1
cccccccccccccccccccccccccccccccccccccc
c     input parameters
      implicit none 
      integer, intent(in) :: istart, iend
      integer, intent(in) :: n  ! the number of processors
c
      integer, dimension(*), intent(out) :: imin, imax, omegaindex
c
      integer :: remainder
      integer :: inum
      integer :: deltai
c
      integer :: i
c
      inum = iend-istart+1
      remainder = mod (inum, n)
      deltai = (inum-remainder)/n
c
      imin (1) = istart
      imax (1) = istart + remainder-1 + deltai
      do i = 2, n
         imin(i)= istart + remainder + (i-1) * deltai
         imax(i)= istart + remainder + i * deltai -1
      enddo

c omegaindex has is actually used only
c for the eventcputimesplit subroutine
      do i = 1,inum
         omegaindex(i) = i
      enddo
      return
      end
      
      subroutine trianglesplit (istart, iend, n, imin, imax,
     & omegaindex)

c     this routine returns the imin(i) and imax(i) (i=1,...,n)
c      separates istart--iend into n parts
c     each part contains imin(i) -- imax(i) i= 1, ..., n
c
c     iend-istart+1 = inum 
c
c     Assume that the cpu time t for i th omega is a * i
c     then we divide a* 0.5 * i **2 into n parts.
c     
c     return imin imax which satisfy above assumption
c
c     
c





cccccccccccccccccccccccccccccccccccccc
c     input parameters
      implicit none 
      integer, intent(in) :: istart, iend
      integer, intent(in) :: n  ! the number of processors
c
      integer, dimension(*), intent(out) :: imin, imax, omegaindex
c
      integer :: remainder
      integer :: inum
      integer :: deltai
c
      integer :: i
      integer, dimension(0:n) :: x
      real(8) :: s  !  0.5 *  iend **2 / n
      real(8) ::p
c
      inum = iend-istart+1
      s =  iend *iend / n
c
      x(0) = istart
      do i = 1, n
         x(i)= s+x(i-1)**2
      !   print *,x(i)
         x(i) = x(i)**0.5
c         p= x(i)**2-x(i-1)**2
c         print *,p
c         print *,x(i)
      enddo

c      print *,x(0)
      do i = 1, n
         imin(i)= x(i-1)
         imax(i)= x(i)-1
c         print *, imin(i),imax(i)
      enddo
      imax (n) = iend
c      print *, x(1), s,n,iend

c omegaindex has is actually used only 
c for the eventcputimesplit subroutine
      do i = 1,inum
         omegaindex(i) = i
      enddo
      return
      end

      subroutine evencputimesplit(istart, iend, n, imin, imax,
     & omegaindex)
c
c Split the frequency points between processors in a cylic fashion
c with period equal to the number of processors (hereafter nproc). Assuming that
c l-cutoff(omega) is not varying significantly within nproc omega points
c after omega, this ensure that the computation time is evenly splitted
c between the processors.
c
      implicit none
      integer, intent(in) :: istart, iend
      integer, intent(in) :: n  ! the number of processors
c
      integer, dimension(*), intent(out) :: imin, imax, omegaindex
c
      integer :: remainder
      integer :: inum
      integer :: deltai, inumround
      integer :: i, j

      inum = iend-istart+1
      deltai = inum/n
      inumround = n*deltai

      do i = 1,n
         imin(i) = istart + (i-1)*deltai + 1
         imax(i) = i*deltai
      enddo
      
      do i = 1,n
         do j = imin(i),imax(i)
            omegaindex(j) = istart + (j-imin(i))*n + (i-1)
         enddo
      enddo
      do i = inumround+1,inum
         omegaindex(i) = istart + i-1
      enddo
      imax(n) = inum
      
c      do i=1,n
c         do j = imin(i),imax(i)
c            write(*,*) i-1,omegaindex(j)
c         enddo
c      enddo

      return
      end
