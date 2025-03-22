!!! this function implements y^2 + 1 on the right hand side of the ODE: y' = y^2 + 1
real*8 function ode_rhs(argument_val)
    implicit none

    real*8 :: argument_val

    ode_rhs = argument_val**2 + 1

end function ode_rhs


program question_1_to_5
    implicit none

    call euler_method(0.001d0, 0.0d0, 1.55d0, 0.0d0, 1, "E:\computational_physics\Module_4_out\question_1.dat")
    call modified_euler_method(0.001d0, 0.0d0, 1.55d0, 0.0d0, 2, "E:\computational_physics\Module_4_out\question_2.dat")
    call improved_euler_method(0.001d0, 0.0d0, 1.55d0, 0.0d0, 3, "E:\computational_physics\Module_4_out\question_3.dat")
    call rk_4(0.01d0, 0.0d0, 1.55d0, 0.0d0, 4, "E:\computational_physics\Module_4_out\question_4.dat")

end program question_1_to_5


!!! Euler Method
subroutine euler_method(step_size, start_x, end_x, y_initial, unit_num, path)
    implicit none

    integer :: i, count, unit_num
    character(len=*) :: path
    real*8 :: temp, temp_2, temp_3
    real*8 :: y_val
    real*8 :: step_size, start_x, end_x, y_initial
    real*8, allocatable :: x_values_list(:), y_values_list(:)
    real*8 :: ode_rhs       !these are function variables

    temp = start_x/step_size
    temp_2 = end_x/step_size
    temp_3 = step_size/step_size
    
    allocate(x_values_list(nint(temp_2 - temp + 1)))
    allocate(y_values_list(nint(temp_2 - temp + 1)))

    x_values_list(1) = start_x      ! saving initial values
    y_values_list(1) = y_initial

    y_val = y_initial + step_size*ode_rhs(y_initial)

    x_values_list(2) = start_x + step_size      ! saving one step after. Beyond this is in the loop
    y_values_list(2) = y_val

    count = 2

    do i = nint(temp), nint(temp_2) - 2, nint(temp_3)
        y_val = y_val + step_size*ode_rhs(y_val)

        count = count + 1
        x_values_list(count) = start_x + (count - 1)*step_size
        y_values_list(count) = y_val
    end do

    open(file = path, unit = unit_num)

    do i = 1, count
        write(unit_num,"(f10.5, a3, f10.5)") x_values_list(i), " , ", y_values_list(i)
    end do

    close(unit_num)

    print*, new_line("a"), "For Euler Method (dx = ", step_size, "), Y_A - Y_E =", 48.078d0 - y_values_list(count)

end subroutine euler_method


!!! Modified Euler Method

! The slope (function of x and y in general) is calculated not at the next x value but 
! the mid point (x_i + x_i+1)/2 of the current x value (say x_i) and the next one (say x_i+1). 
! This slope is then used to calculate slope value at x_i+1 
subroutine modified_euler_method(step_size, start_x, end_x, y_initial, unit_num, path)
    implicit none

    integer :: i, count, unit_num
    character(len=*) :: path
    real*8 :: temp, temp_2, temp_3
    real*8 :: y_val
    real*8 :: step_size, start_x, end_x, y_initial
    real*8 :: x_mid, y_mid
    real*8, allocatable :: x_values_list(:), y_values_list(:)
    real*8 :: ode_rhs       !these are function variables

    temp = start_x/step_size
    temp_2 = end_x/step_size
    temp_3 = step_size/step_size
    
    allocate(x_values_list(nint(temp_2 - temp + 1)))
    allocate(y_values_list(nint(temp_2 - temp + 1)))

    x_values_list(1) = start_x
    y_values_list(1) = y_initial

    y_mid = y_initial + (step_size/2)*ode_rhs(y_initial)
    y_val = y_initial + step_size*ode_rhs(y_mid)

    x_values_list(2) = start_x + step_size
    y_values_list(2) = y_val

    count = 2

    do i = nint(temp), nint(temp_2) - 2, nint(temp_3)
        y_mid = y_val + (step_size/2)*ode_rhs(y_val)
        y_val = y_val + step_size*ode_rhs(y_mid)

        count = count + 1
        x_values_list(count) = start_x + (count - 1)*step_size
        y_values_list(count) = y_val
    end do

    open(file = path, unit = unit_num)

    do i = 1, count
        write(unit_num,"(f10.5, a3, f10.5)") x_values_list(i), " , ", y_values_list(i)
    end do

    close(unit_num)

    print*, new_line("a"), "For Modified Euler Method (dx = ", step_size, "), Y_A - Y_ME =", 48.078d0 - y_values_list(count)

end subroutine modified_euler_method


!!! Improved Euler Method

! The slope (function of x and y in general) is calculated not at the next x value but 
! the mid point (x_i + x_i+1)/2 of the current x value (say x_i) and the next one (say x_i+1). 
! This slope is then used to calculate slope value at x_i+1 
subroutine improved_euler_method(step_size, start_x, end_x, y_initial, unit_num, path)
    implicit none

    integer :: i, count, unit_num
    character(len=*) :: path
    real*8 :: temp, temp_2, temp_3
    real*8 :: y_val
    real*8 :: step_size, start_x, end_x, y_initial
    real*8 :: x_step, y_step
    real*8, allocatable :: x_values_list(:), y_values_list(:)
    real*8 :: ode_rhs       !these are function variables

    temp = start_x/step_size
    temp_2 = end_x/step_size
    temp_3 = step_size/step_size
    
    allocate(x_values_list(nint(temp_2 - temp + 1)))
    allocate(y_values_list(nint(temp_2 - temp + 1)))

    x_values_list(1) = start_x
    y_values_list(1) = y_initial

    y_step = y_initial + step_size*ode_rhs(y_initial)
    y_val = y_initial + step_size*(ode_rhs(y_step) + ode_rhs(y_initial))/2

    x_values_list(2) = start_x + step_size
    y_values_list(2) = y_val

    count = 2

    do i = nint(temp), nint(temp_2) - 2, nint(temp_3)
        y_step = y_val + step_size*ode_rhs(y_val)
        y_val = y_val + step_size*(ode_rhs(y_step) + ode_rhs(y_val))/2

        count = count + 1
        x_values_list(count) = start_x + (count - 1)*step_size
        y_values_list(count) = y_val
    end do

    open(file = path, unit = unit_num)

    do i = 1, count
        write(unit_num,"(f10.5, a3, f10.5)") x_values_list(i), " , ", y_values_list(i)
    end do

    close(unit_num)

    print*, new_line("a"), "For Improved Euler Method (dx = ", step_size, "), Y_A - Y_IE =", 48.078d0 - y_values_list(count)

end subroutine improved_euler_method


!!! Runge Kutta 4 Method
subroutine rk_4(step_size, start_x, end_x, y_initial, unit_num, path)
    implicit none

    integer :: i, count, unit_num
    character(len=*) :: path
    real*8 :: temp, temp_2, temp_3
    real*8 :: y_val
    real*8 :: step_size, start_x, end_x, y_initial
    real*8 :: y_mid, y_mid_better, y_step
    real*8, allocatable :: x_values_list(:), y_values_list(:)
    real*8 :: ode_rhs       !these are function variables

    temp = start_x/step_size
    temp_2 = end_x/step_size
    temp_3 = step_size/step_size
    
    allocate(x_values_list(nint(temp_2 - temp + 1)))
    allocate(y_values_list(nint(temp_2 - temp + 1)))

    x_values_list(1) = start_x
    y_values_list(1) = y_initial

    y_mid = y_initial + (step_size/2)*ode_rhs(y_initial)
    y_mid_better = y_initial + (step_size/2)*ode_rhs(y_mid)
    y_step = y_initial + step_size*ode_rhs(y_mid_better)

    y_val = y_initial + step_size*(2*ode_rhs(y_mid) + 2*ode_rhs(y_mid_better) &
                                    + ode_rhs(y_initial) + ode_rhs(y_step))/6

    x_values_list(2) = start_x + step_size
    y_values_list(2) = y_val

    count = 2

    do i = nint(temp), nint(temp_2) - 2, nint(temp_3)
        y_mid = y_val + (step_size/2)*ode_rhs(y_val)
        y_mid_better = y_val + (step_size/2)*ode_rhs(y_mid)
        y_step = y_val + step_size*ode_rhs(y_mid_better)

        y_val = y_val + step_size*(2*ode_rhs(y_mid) + 2*ode_rhs(y_mid_better) &
                                    + ode_rhs(y_val) + ode_rhs(y_step))/6

        count = count + 1
        x_values_list(count) = start_x + (count - 1)*step_size
        y_values_list(count) = y_val
    end do

    open(file = path, unit = unit_num)

    do i = 1, count
        write(unit_num,"(f10.5, a3, f10.5)") x_values_list(i), " , ", y_values_list(i)
    end do

    close(unit_num)

    print*, new_line("a"), "For Runge Kutta 4 (dx = ", step_size, "), Y_A - Y_RK4 =", 48.078d0 - y_values_list(count)

end subroutine rk_4
