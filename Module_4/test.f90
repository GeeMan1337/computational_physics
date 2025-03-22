!!! this function implements y^2 + 1 on the right hand side of the ODE: y' = y^2 + 1
real*8 function ode_rhs(argument_val)
implicit none

real*8 :: argument_val

ode_rhs = 2!argument_val**2 + 1

end function ode_rhs


program question_1_to_5
implicit none

real*8 :: step_size, start_x, end_x, y_initial

step_size = 0.001d0
start_x = 0.0d0
end_x = 1.55d0
y_initial = 0.0d0

call euler_method(step_size, start_x, end_x, y_initial, 1, "E:\computational_physics\Module_4_out\question_1.dat")

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
    write(1,"(f10.5, a3, f10.5)") x_values_list(i), " , ", y_values_list(i)
end do

close(unit_num)

end subroutine euler_method
