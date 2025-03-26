real*8 function rhs_expression(y_left, y_right, x_val, step_size)
    implicit none

    real*8 :: y_left, y_right, x_val, step_size

    rhs_expression = (-5*y_right*step_size + 2*y_right + 5*y_left*step_size + 2*y_left - 20*(step_size**2)*x_val) &
                    /(4-20*(step_size**2))

end function rhs_expression


program gauss_seidel
    implicit none

    integer :: i, count, flag
    real*8 :: temp, temp_2, temp_3, converge_cond
    real*8 :: step_size, start_x, end_x
    real*8, allocatable :: x_val(:), y_val(:), y_old(:)
    real*8 :: rhs_expression  !these are function variables

    start_x = 0.0d0
    end_x = 1.0d0
    step_size = 0.01d0
    converge_cond = 0.0001d0

    allocate(x_val(nint((end_x - start_x)/step_size) + 1))
    allocate(y_val(nint((end_x - start_x)/step_size) + 1))
    allocate(y_old(nint((end_x - start_x)/step_size) + 1))

    temp = start_x/step_size
    temp_2 = end_x/step_size
    temp_3 = step_size/step_size

    do i = 1, size(x_val)
        x_val(i) = (i-1)*step_size
        y_val(i) = 2*(i-1)*step_size        !using straight line as starting condition
    end do

    y_old = y_val
    count = 0

    do
        flag = 1

        do i = 2, size(x_val) - 1
            y_val(i) = rhs_expression(y_val(i-1), y_val(i+1), x_val(i), step_size)
        end do

        do i = 1, size(y_val)
            if (abs(y_val(i) - y_old(i)) > converge_cond) then
                flag = 0
                exit
            end if
        end do

        if (flag == 1) then
            exit
        end if

        y_old = y_val
        count = count + 1

    end do

    print *, "Number of iterations required =", count
    
    open(file = "E:\computational_physics\Module_4_out\question_9.dat", unit = 10)
    do i = 1, size(x_val)
        write(10,*) x_val(i), " , ", 2*(i-1)*step_size, " , ", y_val(i)
    end do
    close(10)

    do i = 1, size(x_val)
        if (x_val(i) == 0.8d0) then
            print *, "y(0.8) =", y_val(i)
        end if
    end do

end program gauss_seidel