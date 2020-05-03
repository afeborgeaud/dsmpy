program main

    use spctime
    implicit none
    
    integer, parameter :: n=100
    integer :: i
    complex(dpc) :: in(n)
    complex(dpc) :: out(n)

    do i=1,n
        in(i) = (0.d0, 0.d0)
        out(i) = (0.d0, 0.d0)
    enddo

    call fft(in,out,n)

end program main