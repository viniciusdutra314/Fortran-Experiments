program measure_performance
    use fourier_procedures
    complex(kind=8), dimension(:), allocatable :: input
    complex(kind=8), dimension(:), allocatable :: output
    write(*,*) "          N   DFT  (s)         FFT  (s)         speed_up  "
    do i=6,11 
        iarray_size=2**i
        allocate(input(iarray_size),output(iarray_size))
        
        call cpu_time(start_time_DFT)
        output=DFT(input)
        call cpu_time(end_time_DFT)
        DFT_time=end_time_DFT-start_time_DFT

        call cpu_time(start_time_FFT)
        output=FFT(input)
        call cpu_time(end_time_FFT)
        FFT_time=end_time_FFT-start_time_FFT

        write(*,*) iarray_size, DFT_time,FFT_time,DFT_time/FFT_time
        deallocate(input,output)
    end do
end program measure_performance