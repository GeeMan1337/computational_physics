program question_5
    implicit none
    
    call rk_4_2nd_order(0.01d0, 0.0d0, 50.0d0, 0.1d0, 1.9d0, 1, "E:\computational_physics\Module_4_out\question_5.dat")
    call rk_4_2nd_order(0.01d0, 0.0d0, 50.0d0, 0.0d0, 1.999d0, 2, "E:\computational_physics\Module_4_out\question_6.dat")
    call rk_4_2nd_order(0.001d0, 0.0d0, 50.0d0, 0.0d0, 1.999d0, 3, "E:\computational_physics\Module_4_out\question_6b.dat")
    call rk_4_2nd_order(0.01d0, 0.0d0, 50.0d0, 0.0d0, 2.3d0, 4, "E:\computational_physics\Module_4_out\question_6c.dat")

end program question_5
    

!!! this function implements z on the right hand side of the ODE: y' = z
real*8 function ode_rhs_1(argument_val)
    implicit none

    real*8 :: argument_val

    ode_rhs_1 = argument_val    

end function ode_rhs_1


!!! this function implements -sin(y) on the right hand side of the ODE: z' = -sin(y)
real*8 function ode_rhs_2(argument_val)
    implicit none

    real*8 :: argument_val

    ode_rhs_2 = -sin(argument_val)

end function ode_rhs_2



!!! Runge Kutta 4 Method for 2nd order ODE. We split that ODE up into 2 coupled 1st order ODEs.
! here we are solving for y'' = -sin(y)
! break it up into y' = z and z' = -sin(y)
subroutine rk_4_2nd_order(step_size, start_x, end_x, y_initial, z_initial, unit_num, path)
    implicit none

    integer :: i, count, unit_num
    character(len=*) :: path
    real*8 :: temp, temp_2, temp_3
    real*8 :: y_val, z_val
    real*8 :: step_size, start_x, end_x, y_initial, z_initial
    real*8 :: y_mid, y_mid_better, y_step, z_mid, z_mid_better, z_step
    real*8, allocatable :: x_values_list(:), y_values_list(:), z_values_list(:) 
    real*8, allocatable :: kinetic_energy(:), potential_energy(:), total_energy(:)
    real*8 :: ode_rhs_1, ode_rhs_2       !these are function variables

    temp = start_x/step_size
    temp_2 = end_x/step_size
    temp_3 = step_size/step_size
    
    allocate(x_values_list(nint(temp_2 - temp + 1)))
    allocate(y_values_list(nint(temp_2 - temp + 1)))
    allocate(z_values_list(nint(temp_2 - temp + 1)))
    allocate(kinetic_energy(nint(temp_2 - temp + 1)))
    allocate(potential_energy(nint(temp_2 - temp + 1)))
    allocate(total_energy(nint(temp_2 - temp + 1)))


    x_values_list(1) = start_x
    y_values_list(1) = y_initial
    z_values_list(1) = z_initial
    kinetic_energy(1) = 0.5d0*(z_values_list(1)**2)
    potential_energy(1) = -cos(y_values_list(1))
    total_energy(1) = 0.5d0*(z_values_list(1)**2) - cos(y_values_list(1))

    y_mid = y_initial + (step_size/2)*ode_rhs_1(z_initial)
    z_mid = z_initial + (step_size/2)*ode_rhs_2(y_initial)

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
    kinetic_energy(2) = 0.5d0*(z_values_list(2)**2)
    potential_energy(2) = -cos(y_values_list(2))
    total_energy(2) = 0.5d0*(z_values_list(2)**2) - cos(y_values_list(2))


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
        kinetic_energy(count) = 0.5d0*(z_values_list(count)**2)
        potential_energy(count) = -cos(y_values_list(count))
        total_energy(count) = 0.5d0*(z_values_list(count)**2) - cos(y_values_list(count))

    end do

    open(file = path, unit = unit_num)

    do i = 1, count
        write(unit_num,"(f10.5, a3, f10.5, a3, f10.5, a3, f10.5, a3, f10.5, a3, f10.5)") &
                       x_values_list(i), " , ", y_values_list(i), &
                       " , ", z_values_list(i), " , ", kinetic_energy(i), " , ", potential_energy(i), &
                       " , ", total_energy(i)
    end do

    close(unit_num)

    print *, "The value of x after 5000 iterations is =", y_values_list(count)
    print *, "x_0, v_0 and dt are", y_initial, z_initial, step_size, new_line("a")

end subroutine rk_4_2nd_order
