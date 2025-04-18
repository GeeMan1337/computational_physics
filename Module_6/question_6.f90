program question_6
    use gas_sim_better
    implicit none
    
    integer :: i
    real*8, allocatable :: out(:,:)

!subroutine gas_sim_microcanon(num_particles, mass, r_cutoff, epsilon, sigma, temperature, dt, run_time, &
                        !    sample_iter, len, output_values, corr_path, neigh_path)


call gas_simulation_better(2400, 1.0d0, 2.5d0, 1.0d0, 1.0d0, 1.0d0, 0.0025d0, 50000, &
     100, 20.0d0, out, "on", "E:\computational_physics\Module_6_out\question_6_corr.dat", &
     "E:\computational_physics\Module_6_out\question_6_neigh.dat")

open(file = "E:\computational_physics\Module_6_out\question_6.dat", unit = 1)
do i = 1, size(out(1,:))
    write(1,*) out(1,i), out(2,i), out(3,i), out(4,i), out(5,i), out(6,i), out(7,i)
end do
close(1)

print *, "Question 6 done"

deallocate(out)

end program question_6


