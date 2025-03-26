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
    real*8 :: y_val(50), z_val(50)
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

    y_val = 0.0d0
    y_val(1) = 0.8d0
    y_val(26) = 0.8d0
    z_val = 0.0d0

    x_values_list(1) = start_x  
    y_values_list(1,1:50) = y_val   ! 1st row of y_values_list and z_values_list are the initial conditions
    z_values_list(1,1:50) = z_val
   
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

        y_mid(i) = y_val(i) + (step_size/2)*ode_rhs_1(z_val(i))
        z_mid(i) = z_val(i) + (step_size/2)* &
                    ode_rhs_2(y_val(left_index),y_val(i),y_val(right_index))

        y_mid_better(i) = y_val(i) + (step_size/2)*ode_rhs_1(z_mid(i))
        z_mid_better(i) = z_val(i) + (step_size/2)*ode_rhs_2(y_mid(left_index),y_mid(i),y_mid(right_index))

        y_step(i) = y_val(i) + step_size*ode_rhs_1(z_mid_better(i))
        z_step(i) = z_val(i) + step_size*ode_rhs_2(y_mid_better(left_index),y_mid_better(i),y_mid_better(right_index))

        y_val(i) = y_val(i) + step_size*(2*ode_rhs_1(z_mid(i)) + 2*ode_rhs_1(z_mid_better(i)) &
                                        + ode_rhs_1(z_val(i)) + ode_rhs_1(z_step(i)))/6

        z_val(i) = z_val(i) + step_size*(2*ode_rhs_2(y_mid(left_index),y_mid(i),y_mid(right_index)) &
                                        + 2*ode_rhs_2(y_mid_better(left_index),y_mid_better(i),y_mid_better(right_index)) &
                                        + ode_rhs_2(y_val(left_index),y_val(i),y_val(right_index)) &
                                        + ode_rhs_2(y_step(left_index),y_step(i),y_step(right_index)))/6
    end do  

    x_values_list(2) = start_x + step_size
    y_values_list(2,1:50) = y_val
    z_values_list(2,1:50) = z_val

    count = 2

    do j = nint(temp), nint(temp_2) - 2, nint(temp_3)
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
            
            y_mid(i) = y_val(i) + (step_size/2)*ode_rhs_1(z_val(i))
            z_mid(i) = z_val(i) + (step_size/2)* &
                        ode_rhs_2(y_val(left_index),y_val(i),y_val(right_index))
    
            y_mid_better(i) = y_val(i) + (step_size/2)*ode_rhs_1(z_mid(i))
            z_mid_better(i) = z_val(i) + (step_size/2)*ode_rhs_2(y_mid(left_index),y_mid(i),y_mid(right_index))
    
            y_step(i) = y_val(i) + step_size*ode_rhs_1(z_mid_better(i))
            z_step(i) = z_val(i) + step_size*ode_rhs_2(y_mid_better(left_index),y_mid_better(i),y_mid_better(right_index))
    
            y_val(i) = y_val(i) + step_size*(2*ode_rhs_1(z_mid(i)) + 2*ode_rhs_1(z_mid_better(i)) &
                                            + ode_rhs_1(z_val(i)) + ode_rhs_1(z_step(i)))/6
    
            z_val(i) = z_val(i) + step_size*(2*ode_rhs_2(y_mid(left_index),y_mid(i),y_mid(right_index)) &
                                            + 2*ode_rhs_2(y_mid_better(left_index),y_mid_better(i),y_mid_better(right_index)) &
                                            + ode_rhs_2(y_val(left_index),y_val(i),y_val(right_index)) &
                                            + ode_rhs_2(y_step(left_index),y_step(i),y_step(right_index)))/6
        end do
        count = count + 1
        x_values_list(count) = start_x + step_size
        y_values_list(count,1:50) = y_val
        z_values_list(count,1:50) = z_val
    end do

    open(file = "E:\computational_physics\Module_4_out\question_8.xyz", unit = 10)

    do i = 1, count
        write(10,*) 50, new_line("a"), "Timestep = ", i-1
        
        do j = 1, 50
            write(10,*) j, 5*sin(2*3.14159265359*j/50), y_values_list(i,j), 5*cos(2*3.14159265359*j/50)
        end do 
    end do

    close(10)

end program rk_4_2nd_order_50_equations
