real function fun(in_function, x) result(out_value)
    implicit none 

    real :: in_function, x
    
    out_value = in_function

end function fun


program output 
    implicit none 
    real :: x, fun
    integer :: value
    character :: expression(50)

    expression = "3+3"

    read(expression, *) value
    print *, value

end program output