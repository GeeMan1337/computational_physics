program question_2_and_3
    use gas_sim
    implicit none
    
    integer :: i
    real*8, allocatable :: out(:,:)

!subroutine gas_sim_microcanon(num_particles, mass, r_cutoff, epsilon, sigma, temperature, dt, run_time, &
                        !    sample_iter, len, output_values)


call gas_simulation(2197, 1.0d0, 2.5d0, 1.0d0, 1.0d0, 1.0d0, 0.005d0, 20000, 50, 20.0d0, out, "off")

open(file = "E:\computational_physics\Module_6_out\question_2.dat", unit = 1)
do i = 1, size(out(1,:))
    write(1,*) out(1,i), out(2,i), out(3,i), out(4,i), out(5,i), out(6,i), out(7,i)
end do
close(1)

print *, "Question 2 done"

deallocate(out)

call gas_simulation(2197, 1.0d0, 2.5d0, 1.0d0, 1.0d0, 1.0d0, 0.005d0, 20000, 50, 20.0d0, out, "on")

open(file = "E:\computational_physics\Module_6_out\question_3.dat", unit = 2)
do i = 1, size(out(1,:))
    write(2,*) out(1,i), out(2,i), out(3,i), out(4,i), out(5,i), out(6,i), out(7,i)
end do
close(2)

print *, "Question 3 done"

end program question_2_and_3


