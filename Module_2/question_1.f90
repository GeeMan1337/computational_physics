program trapezoidal_integral
    implicit none 

    integer :: i, n, index, exponent
    real*8 :: x, dx
    real*8 :: function_1, upper_lim, lower_lim, integral
    real*8 :: integral_list(10), error_list(10), n_list(10), dx_list(10)
    real*8, parameter :: pi = 2.0d0 * asin(1.0d0)

    !!! question 1a and 1b

    index = 0
    upper_lim = 1.0d0
    lower_lim = 0.0d0

    do exponent = 1, 7
        
        index = index + 1
        dx = 10.0d0**(-real(exponent))
        integral = 0d0

        ! number of subdivisions is n

        n = int((upper_lim - lower_lim)/dx)

        if (n <= 0) then
            print *, "Error! Number of subdivisions is not positive"
            print *, dx
        end if

        do i = 1, n
            x = lower_lim + i*dx
            function_1 = 4/(1+x**2)
            integral = integral + dx*function_1
        end do

        x = lower_lim
        function_1 = 4/(1+x**2)
        integral = integral + (dx/2) * function_1 

        x = upper_lim
        function_1 = 4/(1+x**2)
        integral = integral + ((upper_lim - lower_lim - n*dx)/2) * function_1 

        integral_list(index) = integral
        error_list(index) = abs(integral - pi)
        n_list(index) = n
    
    end do


    !!! question 1c

    index = 0
    upper_lim = pi
    lower_lim = 0.0d0

    do exponent = 1, 6
        
        index = index + 1
        dx = 10.0d0**(-real(exponent))
        integral = 0d0

        ! number of subdivisions is n

        n = int((upper_lim - lower_lim)/dx)

        if (n <= 0) then
            print *, "Error! Number of subdivisions is not positive"
            print *, dx
        end if

        do i = 1, n
            x = lower_lim + i*dx
            function_1 = sin(x)
            integral = integral + dx*function_1
        end do

        x = lower_lim
        function_1 = sin(x)
        integral = integral + (dx/2) * function_1 

        x = upper_lim
        function_1 = sin(x)
        integral = integral + ((upper_lim - lower_lim - n*dx)/2) * function_1 

        integral_list(index) = integral
        error_list(index) = abs(integral - 2.0d0)
        n_list(index) = n
    
    end do


    !!! question 1d

    index = 0
    upper_lim = 3.0d0
    lower_lim = -3.0d0

    do exponent = 1, 7
        
        index = index + 1
        dx = 10.0d0**(-real(exponent))
        integral = 0d0

        ! number of subdivisions is n

        n = int((upper_lim - lower_lim)/dx)
        print *, n
        if (n <= 0) then
            print *, "Error! Number of subdivisions is not positive"
            print *, dx
        end if

        do i = 1, n
            x = lower_lim + i*dx
            function_1 = 1/sqrt(2*pi) * exp(-(x**2)/2.0d0)
            integral = integral + dx*function_1
        end do

        x = lower_lim
        function_1 = 1/sqrt(2*pi) * exp(-(x**2)/2.0d0)
        integral = integral + (dx/2) * function_1 

        x = upper_lim
        function_1 = 1/sqrt(2*pi) * exp(-(x**2)/2.0d0)
        integral = integral + ((upper_lim - lower_lim - n*dx)/2) * function_1 

        integral_list(index) = integral
        error_list(index) = abs(integral - 0.9973d0)
        n_list(index) = n
        dx_list(index) = dx

    end do
    print *, error_list(1:index)
    print *, 0.9973d0
end program trapezoidal_integral

!subroutine graph_data()