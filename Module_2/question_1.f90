real function trapezoidal_integral(function, step_size, upper_lim, lower_lim) result(integral)
    implicit none

    integer :: i, j
    real*8 :: x
    real*8 :: function, step_size, upper_lim, lower_lim, integral

    integral = 0

    do i = 1, int((upper_lim - lower_lim)/step_size)
        integral = integral + 



end function trapezoidal_lintegral


real function fun(in_function, in_value) result(out_value)
    implicit none 

    real*8 :: in_value, function, x
    
    x = in_value
    out_value = in_function

end function fun