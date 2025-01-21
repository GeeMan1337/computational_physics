real*8 function fun(in_function, in_value) result(out_value)
    implicit none 

    real*8 :: in_value, in_function, x
    
    x = in_value
    out_value = in_function

end function fun


program output 
    implicit none 
    real*8 :: fun, x

    call fun(1+x, 3)


end program output