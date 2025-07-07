module fourier_procedures
    real(kind=8) :: pi = 3.1415926535897932384626433
    contains
        function DFT(input_domain) result(output_domain)
            complex(kind=8), dimension(:), intent(in) :: input_domain
            complex(kind=8), dimension(size(input_domain)) :: output_domain
            N = size(input_domain)
            do k = 0, N-1
                output_domain(k+1) = (0.0, 0.0)
                do j = 0, N-1
                    output_domain(k+1) = output_domain(k+1) + & 
                    input_domain(j+1) *exp(cmplx(0.0,(-2.0d0 * pi * k * j) / N,kind=8))
                end do
            end do
        end function 

        function DFT_inverse(input_domain) result(output_domain)
            complex(kind=8), dimension(:), intent(in) :: input_domain
            complex(kind=8), dimension(size(input_domain)) :: output_domain
            N = size(input_domain)
            do k = 0, N-1
                output_domain(k+1) = (0.0, 0.0)
                do j = 0, N-1
                    output_domain(k+1) = output_domain(k+1) + &
                    input_domain(j+1) *exp(cmplx(0.0,(2.0 * pi * k * j) / N,kind=8))
                end do
                output_domain(k+1)=output_domain(k+1)/N
            end do
        end function

        recursive function FFT(x) result(y)
            complex(kind=8), dimension(:), intent(in) :: x
            complex(kind=8), dimension(size(x)) :: y
            complex(kind=8), dimension(size(x)/2) :: y_even,y_odd
            complex(kind=8) :: omega
            n=size(x)
            if (n==1) then
                y=x
                return 
            end if
            omega=exp(cmplx(0,-2*pi/n,kind=8))
            y_even=FFT(x(1::2))
            y_odd=FFT(x(2::2))
            y=(0.0,0.0)
            do j=1,n/2
                y(j)=y_even(j) +(omega**(j-1))*y_odd(j)
                y(j+n/2)=y_even(j) -(omega**(j-1))*y_odd(j)
            end do 
        end function
        
        recursive function FFT_inverse_not_normalized(x) result(y)
            complex(kind=8), dimension(:), intent(in) :: x
            complex(kind=8), dimension(size(x)) :: y
            complex(kind=8), dimension(size(x)/2) :: y_even,y_odd
            complex(kind=8) :: omega
            n=size(x)
            if (n==1) then
                y=x
                return 
            end if
            omega=exp(cmplx(0,2*pi/n,kind=8))
            y_even=FFT_inverse_not_normalized(x(1::2))
            y_odd=FFT_inverse_not_normalized(x(2::2))
            y=(0.0,0.0)
            do j=1,n/2
                y(j)=(y_even(j) +(omega**(j-1))*y_odd(j))
                y(j+n/2)=(y_even(j) -(omega**(j-1))*y_odd(j))
            end do 
        end function

        function FFT_inverse(x) result (y)
            complex(kind=8), dimension(:), intent(in) :: x
            complex(kind=8), dimension(size(x)) :: y
            y=FFT_inverse_not_normalized(x)/size(x)
        end function

end module fourier_procedures