real*8 function rhs_expression(x_left, x_right, y_up, y_down, x_step, y_step)
   implicit none 

   real*8 :: x_left, x_right, y_up, y_down, x_step, y_step

   rhs_expression = (x_step**2)*(y_step**2)*((x_left + x_right)/x_step**2 + (y_down + y_up)/y_step**2) &
                  /(2*(x_step**2 + y_step**2))

end function rhs_expression


program laplace_solver_dirichlet_boundary_conditions
   implicit none 

   integer :: i, j, flag, count
   real*8 :: temp_old(34,34), temp_val(34,34)
   real*8 :: x_step, y_step, converge_cond
   real*8 :: rhs_expression      !function variable

   x_step = 1
   y_step = 1
   converge_cond = 0.0001d0

   temp_val = 0
   temp_val(1,1:34) = 3.7
   temp_val(34,1:34) = 0.4

   do i = 1, 34
      temp_val(i,1) = i*(temp_val(34,1) - temp_val(1,1))/(34 - 1) &
                     + (temp_val(1,1)*34 - temp_val(34,1)*1)/(34 - 1)

      temp_val(i,34) = i*(temp_val(34,34) - temp_val(1,34))/(34 - 1) &
                        + (temp_val(1,34)*34 - temp_val(34,34)*1)/(34 - 1)
   end do

   temp_old = temp_val

   open(file = "E:\computational_physics\Module_4_out\pde_question_3a.dat", unit = 10)
   do i = 34, 1, -1
      write (10,*) temp_val(1:34,i)
   end do
   close(10)

   count = 0
   
   do

      flag = 1

      do i = 2, 33
         do j = 2, 33
            temp_val(i,j) = rhs_expression(temp_val(i-1,j), temp_val(i+1,j), temp_val(i,j+1), temp_val(i,j-1), 1.0d0, 1.0d0)
         end do
      end do

      do i = 1, 34
         do j = 1, 34
            if (abs(temp_val(i,j) - temp_old(i,j)) > converge_cond) then
               flag = 0
               go to 99
            end if
         end do
      end do

99    if (flag == 1) then
         exit
      end if

      count = count + 1
      temp_old = temp_val

   end do


   open(file = "E:\computational_physics\Module_4_out\pde_question_3b.dat", unit = 11)
   do i = 1, 34
      write (11,*) temp_val(1:34,i)
   end do
   close(11)

   print *, "Number of iterations required (34 X 34 lattice)=", count
   print *, new_line("a"), "Temperature(20,20) =", temp_val(20,20)

end program laplace_solver_dirichlet_boundary_conditions