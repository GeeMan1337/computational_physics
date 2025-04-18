program question_4
    use gas_sim_better
    implicit none
    
    integer :: i
    real*8, allocatable :: out(:,:)

!subroutine gas_sim_microcanon(num_particles, mass, r_cutoff, epsilon, sigma, temperature, dt, run_time, &
                        !    sample_iter, len, output_values, corr_path, neigh_path)


call gas_simulation_better(1200, 1.0d0, 2.5d0, 1.0d0, 1.0d0, 1.0d0, 0.0025d0, 250000, &
     50, 20.0d0, out, "on", "E:\computational_physics\Module_6_out\question_4_corr.dat", &
     "E:\computational_physics\Module_6_out\question_4_neigh.dat")

open(file = "E:\computational_physics\Module_6_out\question_4.dat", unit = 1)
do i = 1, size(out(1,:))
    write(1,*) out(1,i), out(2,i), out(3,i), out(4,i), out(5,i), out(6,i), out(7,i)
end do
close(1)

print *, "Question 4 done"

deallocate(out)

end program question_4


