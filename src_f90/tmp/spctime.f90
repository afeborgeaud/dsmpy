module spctime

    integer, parameter :: dpc = kind((1.d0,1.d0))
    
contains

    subroutine fft(in, out, n)
        integer, intent(in) :: n
        complex(dpc), intent(in) :: in(n)
        complex(dpc), intent(out) :: out(n)
        integer :: plan, i

        call dfftw_plan_dft_1d(plan,n,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,in,out)
        call dfftw_destroy_plan(plan)

        do i = 1,n
            write(*,*) in(i), out(i)
        enddo

        return
    end subroutine

end module spctime