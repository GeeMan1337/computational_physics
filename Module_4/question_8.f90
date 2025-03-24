real*8 function ode_rhs_1(argument_val)
    implicit none

    real*8 :: argument_val

    ode_rhs_1 = argument_val 

end function ode_rhs_1


real*8 function ode_rhs_2(argument_val_1, argument_val_2, argument_val_3)
    implicit none

    real*8 :: argument_val_1, argument_val_2, argument_val_3

    ode_rhs_2 = argument_val_1 + argument_val_3 - 2*argument_val_2

end function ode_rhs_2


!!! Runge Kutta 4 Method for 2nd order ODE. We split that ODE up into coupled 1st order ODEs.
program rk_4_2nd_order_50_equations
    implicit none

    integer :: i, j, count
    real*8 :: temp, temp_2, temp_3
    real*8 :: y_val(50)
    real*8 :: step_size, start_x, end_x
    integer :: left_index, right_index
    real*8 :: y_mid(50), y_mid_better(50), y_step(50), z_mid(50), z_mid_better(50), z_step(50)
    real*8, allocatable :: x_values_list(:), y_values_list(:,:), z_values_list(:,:) 
    real*8 :: ode_rhs_1, ode_rhs_2       !these are function variables

    step_size = 0.02d0
    start_x = 0.0d0
    end_x = 40.0d0

    temp = start_x/step_size
    temp_2 = end_x/step_size
    temp_3 = step_size/step_size

    allocate(x_values_list(nint(temp_2 - temp + 1)))
    allocate(y_values_list(nint(temp_2 - temp + 1), 50))
    allocate(z_values_list(nint(temp_2 - temp + 1), 50))

    x_values_list(1) = start_x
    y_values_list = 0.0d0               ! 1st row of y_values_list and z_values_list are the initial conditions
    y_values_list(1,1) = 0.8d0
    y_values_list(1,26) = 0.8d0
    z_values_list = 0.0d0
   
    do i = 1, 50

        if (i == 1) then
            left_index = 50 
            right_index = 2
        elseif (i == 50) then
            left_index = 49 
            right_index = 1
        else
            left_index = i-1
            right_index = i+1
        end if

        y_mid(i) = y_values_list(1,i) + (step_size/2)*ode_rhs_1(z_values_list(1,i))
        z_mid(i) = z_values_list(1,i) + (step_size/2)* &
                    ode_rhs_2(y_initial(1,left_index),y_initial(1,i),y_initial(1,right_index))

        y_mid_better = y_initial + (step_size/2)*ode_rhs_1(z_mid)
        z_mid_better = z_initial + (step_size/2)*ode_rhs_2(y_mid)

        y_step = y_initial + step_size*ode_rhs_1(z_mid_better)
        z_step = z_initial + step_size*ode_rhs_2(y_mid_better)

        y_val = y_initial + step_size*(2*ode_rhs_1(z_mid) + 2*ode_rhs_1(z_mid_better) &
                                        + ode_rhs_1(z_initial) + ode_rhs_1(z_step))/6

        z_val = z_initial + step_size*(2*ode_rhs_2(y_mid) + 2*ode_rhs_2(y_mid_better) &
                                        + ode_rhs_2(y_initial) + ode_rhs_2(y_step))/6



        x_values_list(2) = start_x + step_size
        y_values_list(2) = y_val
        z_values_list(2) = z_val

    end do  

    count = 2

    do i = nint(temp), nint(temp_2) - 2, nint(temp_3)
        y_mid = y_val + (step_size/2)*ode_rhs_1(z_val)
        z_mid = z_val + (step_size/2)*ode_rhs_2(y_val)

        y_mid_better = y_val + (step_size/2)*ode_rhs_1(z_mid)
        z_mid_better = z_val + (step_size/2)*ode_rhs_2(y_mid)

        y_step = y_val + step_size*ode_rhs_1(z_mid_better)
        z_step = z_val + step_size*ode_rhs_2(y_mid_better)

        y_val = y_val + step_size*(2*ode_rhs_1(z_mid) + 2*ode_rhs_1(z_mid_better) &
                                    + ode_rhs_1(z_val) + ode_rhs_1(z_step))/6
                                    
        z_val = z_val + step_size*(2*ode_rhs_2(y_mid) + 2*ode_rhs_2(y_mid_better) &
                                    + ode_rhs_2(y_val) + ode_rhs_2(y_step))/6

        count = count + 1
        x_values_list(count) = start_x + (count - 1)*step_size
        y_values_list(count) = y_val
        z_values_list(count) = z_val

    end do

    open(file = path, unit = unit_num)

    do i = 1, count
        write(unit_num,"(f10.5, a3, f10.5, a3, f10.5, a3, f10.5, a3, f10.5, a3, f10.5)") &
                    x_values_list(i), " , ", y_values_list(i), &
                    " , ", z_values_list(i), " , ", kinetic_energy(i), " , ", potential_energy(i), &
                    " , ", total_energy(i)
    end do

    close(unit_num)

end program rk_4_2nd_order_50_equations
