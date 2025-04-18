module gas_sim
    implicit none 
    
    contains


    subroutine gas_simulation(num_particles, mass, r_cutoff, epsilon, sigma, temperature, dt, run_time, &
                                  sample_iter, len, output_values, thermostat)
        
        integer :: i, j, k
        integer :: num_particles, run_time, count
        real*8 :: temp, temp_2, temp_3
        real*8 :: len
        real*8 :: epsilon, sigma, mass, r_cutoff, dt                ! dt is time step, run_time is number of iterations
        real*8 :: temperature, vel_coeff                            ! temperature = KT
        real*8, allocatable :: r_list(:), v_list(:), f_new_list(:), f_list(:) ! x=pos, v=vel, f=force
        real*8, allocatable :: output_values(:,:)
        real*8 :: sigma_6, sigma_12, dt2by2m                     ! this are vars to simplify calcs
        real*8 :: force, r_cutoff_12, r_cutoff_6, force_cutoff   ! this are vars to simplify calcs
        real*8 :: x_dist, y_dist, z_dist, r_dist
        integer :: sample_iter, num_cells
        character(len=*) :: thermostat

        sigma_6 = sigma**6; sigma_12 = sigma**12
        r_cutoff_6 = r_cutoff**6; r_cutoff_12 = r_cutoff**12
        force_cutoff = 4*epsilon*(12*sigma_12/(r_cutoff**13) - 6*sigma_6/(r_cutoff**7))
        
        dt2by2m = 0.5d0*dt*dt/mass
        vel_coeff = dsqrt(12.0d0*temperature/mass)

        allocate(output_values(7, run_time/sample_iter))    
        allocate(r_list(3*num_particles))
        allocate(v_list(3*num_particles))
        allocate(f_list(3*num_particles))
        allocate(f_new_list(3*num_particles))
        f_new_list = 0.0d0; f_list = 0.0d0

        !!! initializing x values
        count = 0
        num_cells = int(num_particles**(1.0/3)) + 1     ! this is the num of cells in one direction 

        do i = 1, num_cells
            do j = 1, num_cells
                do k = 1, num_cells

                    count = count + 1
                    r_list(3*count) = modulo(k*len/num_cells, len)
                    r_list(3*count-1) = modulo(j*len/num_cells, len)
                    r_list(3*count-2) = modulo(i*len/num_cells, len)

                    if (count == num_particles) then
                        goto 63
                    end if

                end do
            end do
        end do

        !!! initializing v values
63      do i = 1, 3*num_particles
            call random_number(temp)
            v_list(i) = temp - 0.5d0
        end do

        temp = 0.0d0; temp_2 = 0.0d0; temp_3 = 0.0d0

        do i = 1, num_particles
            temp = temp + v_list(3*i-2)/num_particles       !these are average vel in x,y,z directions
            temp_2 = temp_2 + v_list(3*i-1)/num_particles
            temp_3 = temp_3 + v_list(3*i)/num_particles
        end do
    
        do i = 1, num_particles
            v_list(3*i) = vel_coeff*(v_list(3*i) - temp_3)      ! this is done to make the net velocit zero 
            v_list(3*i-1) = vel_coeff*(v_list(3*i-1) - temp_2)
            v_list(3*i-2) = vel_coeff*(v_list(3*i-2) - temp)
        end do                       

        !!! this loop calculates initial force 
        do j = 1, num_particles - 1
            do k = j + 1, num_particles

                x_dist = r_list(3*j-2) - r_list(3*k-2)
                y_dist = r_list(3*j-1) - r_list(3*k-1)
                z_dist = r_list(3*j  ) - r_list(3*k  )

                x_dist = x_dist - len * nint(x_dist / len) 
                y_dist = y_dist - len * nint(y_dist / len)
                z_dist = z_dist - len * nint(z_dist / len)
                
                r_dist = dsqrt(x_dist**2 + y_dist**2 + z_dist**2)

                if (r_dist <= r_cutoff) then
                    force = 4*epsilon*(12*sigma_12/(r_dist**13) - 6*sigma_6/(r_dist**7)) - force_cutoff
                    f_list(3*j-2) = f_list(3*j-2) + (x_dist/r_dist)*force
                    f_list(3*j-1) = f_list(3*j-1) + (y_dist/r_dist)*force
                    f_list(3*j)   = f_list(3*j)   + (z_dist/r_dist)*force
                    f_list(3*k-2) = f_list(3*k-2) - (x_dist/r_dist)*force
                    f_list(3*k-1) = f_list(3*k-1) - (y_dist/r_dist)*force
                    f_list(3*k)   = f_list(3*k)   - (z_dist/r_dist)*force
                end if

            end do
        end do

        count = 0
        !!! this is where all the iterations happen
        do i = 1, run_time

            do j = 1, 3*num_particles
                r_list(j) = r_list(j) + v_list(j)*dt + f_list(j)*dt2by2m
                r_list(j) = modulo(r_list(j), len)
            end do

            !!! this loop calculates force for the next time step
            do j = 1, num_particles - 1
                do k = j + 1, num_particles

                    x_dist = r_list(3*j-2) - r_list(3*k-2)
                    y_dist = r_list(3*j-1) - r_list(3*k-1)
                    z_dist = r_list(3*j) - r_list(3*k)
                    
                    x_dist = x_dist - len * nint(x_dist / len)
                    y_dist = y_dist - len * nint(y_dist / len)
                    z_dist = z_dist - len * nint(z_dist / len)
                    
                    r_dist = dsqrt(x_dist**2 + y_dist**2 + z_dist**2)

                    if (r_dist <= r_cutoff) then
                        force = 4*epsilon*(12*sigma_12/(r_dist**13) - 6*sigma_6/(r_dist**7)) - force_cutoff
                        f_new_list(3*j-2) = f_new_list(3*j-2) + (x_dist/r_dist)*force
                        f_new_list(3*j-1) = f_new_list(3*j-1) + (y_dist/r_dist)*force
                        f_new_list(3*j)   = f_new_list(3*j)   + (z_dist/r_dist)*force
                        f_new_list(3*k-2) = f_new_list(3*k-2) - (x_dist/r_dist)*force
                        f_new_list(3*k-1) = f_new_list(3*k-1) - (y_dist/r_dist)*force
                        f_new_list(3*k)   = f_new_list(3*k)   - (z_dist/r_dist)*force
                    end if

                end do
            end do

            do j = 1, 3*num_particles
                v_list(j) = v_list(j) + 0.5d0*dt*(f_new_list(j) + f_list(j))/mass
            end do

            if (mod(i, 100) == 0) then
                call thermal_equil(thermostat, num_particles, v_list, mass, temperature)
            end if

            if (mod(i,sample_iter) == 0) then
                count = count + 1
                output_values(1,count) = i
                call calc_kinetic_energy(num_particles, v_list, mass, output_values(2,count))
                call calc_potential_energy(num_particles, r_list, len, r_cutoff, epsilon, sigma, output_values(3,count))
                output_values(4,count) = output_values(2,count) + output_values(3,count)
                call calc_momentum(num_particles, v_list, mass, &
                output_values(5,count), output_values(6,count), output_values(7,count))
            end if 

            f_list = f_new_list
            f_new_list = 0.0d0


        end do

    end subroutine gas_simulation


    subroutine calc_kinetic_energy(num_particles, v_list, mass, output)
        
        integer :: i, num_particles
        real*8 :: v_list(3*num_particles), mass
        real*8 :: output

        output = 0.0d0

        do i = 1, 3*num_particles
            output = output + 0.5d0*mass*(v_list(i)**2)
        end do

        output = output/num_particles
        
    end subroutine calc_kinetic_energy


    subroutine calc_potential_energy(num_particles, r_list, len, r_cutoff, epsilon, sigma, output)
        integer, intent(in) :: num_particles
        real*8, dimension(:), intent(in) :: r_list
        real*8, intent(in) :: len, r_cutoff, epsilon, sigma
        real*8, intent(out) :: output

        integer :: j, k
        real*8 :: sigma_6, sigma_12, x_dist, y_dist, z_dist, r_dist
        real*8 :: U, U_cut, F_cut, V_shifted

        sigma_6 = sigma**6
        sigma_12 = sigma**12

        U_cut = 4.d0*epsilon*(sigma_12/r_cutoff**12 - sigma_6/r_cutoff**6)
        F_cut = 24.d0*epsilon*(2.d0*sigma_12/r_cutoff**13 - sigma_6/r_cutoff**7)

        output = 0.d0

        do j = 1, num_particles - 1
            do k = j + 1, num_particles

                x_dist = r_list(3*j-2) - r_list(3*k-2)
                y_dist = r_list(3*j-1) - r_list(3*k-1)
                z_dist = r_list(3*j) - r_list(3*k)

                x_dist = x_dist - len * nint(x_dist / len)
                y_dist = y_dist - len * nint(y_dist / len)
                z_dist = z_dist - len * nint(z_dist / len)

                r_dist = dsqrt(x_dist**2 + y_dist**2 + z_dist**2)

                if (r_dist < r_cutoff) then
                    U = 4.d0*epsilon*(sigma_12/r_dist**12 - sigma_6/r_dist**6)
                    V_shifted = U - U_cut - F_cut*(r_dist - r_cutoff)
                    output = output + V_shifted
                end if
            end do
        end do

        output = output / num_particles

    end subroutine calc_potential_energy


    subroutine calc_momentum(num_particles, v_list, mass, p_x, p_y, p_z)

        integer :: i, num_particles
        real*8 :: v_list(3*num_particles), mass
        real*8 :: p_x, p_y, p_z

        p_x = 0.0d0; p_y = 0.0d0; p_z = 0.0d0

        do i = 1, num_particles
            p_x = p_x + mass*v_list(3*i-2)
            p_y = p_y + mass*v_list(3*i-1)
            p_z = p_z + mass*v_list(3*i)
        end do

    end subroutine calc_momentum


    subroutine thermal_equil(thermostat, num_particles, v_list, mass, temperature)

        integer :: num_particles
        character(len=*) :: thermostat
        real*8 :: theroy_KE_per_particle, KE_per_particle, v_list(3*num_particles), temperature
        real*8 :: mass, scale_factor

        theroy_KE_per_particle = 1.5d0*temperature
        call calc_kinetic_energy(num_particles, v_list, mass, KE_per_particle)
        scale_factor = dsqrt(theroy_KE_per_particle/KE_per_particle)

        if (thermostat == "on") then
            v_list = scale_factor*v_list
        end if

    end subroutine thermal_equil


end module gas_sim