program trapezoidal_integral
    implicit none

    integer :: i, n
    real*8 :: x, dx
    real*8 :: function_1, function_2, upper_lim, lower_lim, integral
    real*8 :: integral_list(4), error_list(4)
    real*8, parameter :: pi = 2*asin(1.0d0)
    
    upper_lim = 1.0d0
    lower_lim = 0.0d0

    !!! Question 1a
    ! dx = 0.01

    dx = 0.01d0
    integral = 0d0

    ! number of subdivisions is n

    n = int((upper_lim - lower_lim)/dx)

    if (n <= 0) then
        print *, "Error! Number of subdivisions is not positive"
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

    integral_list(1) = integral
    error_list(1) = abs(integral - pi)


    ! dx = 0.001

    dx = 0.001d0
    integral = 0d0

    ! number of subdivisions is n

    n = int((upper_lim - lower_lim)/dx)

    if (n <= 0) then
        print *, "Error! Number of subdivisions is not positive"
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

    integral_list(2) = integral
    error_list(2) = abs(integral - pi)


    ! dx = 0.0001

    dx = 0.0001d0
    integral = 0d0

    ! number of subdivisions is n

    n = int((upper_lim - lower_lim)/dx)

    if (n <= 0) then
        print *, "Error! Number of subdivisions is not positive"
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

    integral_list(3) = integral
    error_list(3) = abs(integral - pi)


    ! dx = 0.00001

    dx = 0.00001d0
    integral = 0d0

    ! number of subdivisions is n

    n = int((upper_lim - lower_lim)/dx)

    if (n <= 0) then
        print *, "Error! Number of subdivisions is not positive"
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

    integral_list(4) = integral
    error_list(4) = abs(integral - pi)

    print *, integral_list
    print *, error_list

end program trapezoidal_integral