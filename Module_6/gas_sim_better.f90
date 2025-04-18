module gas_sim_better
    implicit none 
    
    contains


    subroutine gas_simulation_better(num_particles, mass, r_cutoff, epsilon, sigma, temperature, dt, run_time, &
                    sample_iter, len, output_values, thermostat, corr_path, neigh_path)
        
        integer :: i, j, k
        integer :: num_particles, run_time, count, count_2, count_3
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
        character(len=*) :: thermostat, corr_path, neigh_path
        real*8 :: r_nearby 
        integer, allocatable :: num_neighbour(:), neighbour_list(:,:)
        real*8 :: dr, r_max
        real*8, allocatable :: correlation_list(:,:)
        integer :: num_bins 
        integer, allocatable :: max_neighbours(:,:)

        r_nearby = 4.5d0
        sigma_6 = sigma**6; sigma_12 = sigma**12
        r_cutoff_6 = r_cutoff**6; r_cutoff_12 = r_cutoff**12
        force_cutoff = 4*epsilon*(12*sigma_12/(r_cutoff**13) - 6*sigma_6/(r_cutoff**7))

        dr = 0.1d0
        r_max = 0.5d0*len
        num_bins = int(r_max/dr)
        
        dt2by2m = 0.5d0*dt*dt/mass
        vel_coeff = dsqrt(12.0d0*temperature/mass)

        allocate(output_values(7, run_time/sample_iter))    
        allocate(r_list(3*num_particles))
        allocate(v_list(3*num_particles))
        allocate(f_list(3*num_particles))
        allocate(f_new_list(3*num_particles))
        allocate(num_neighbour(num_particles))
        allocate(neighbour_list(num_particles, num_particles))
        allocate(correlation_list(2, num_bins))
        allocate(max_neighbours(2, 251))


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

        call neighbours(num_particles, r_list, len, r_nearby, num_neighbour, neighbour_list)

        do j = 1, num_particles 
            do k = 1, num_neighbour(j)

                x_dist = r_list(3*j-2) - r_list(3*neighbour_list(j,k)-2)
                y_dist = r_list(3*j-1) - r_list(3*neighbour_list(j,k)-1)
                z_dist = r_list(3*j) - r_list(3*neighbour_list(j,k))

                x_dist = x_dist - len * nint(x_dist / len) 
                y_dist = y_dist - len * nint(y_dist / len)
                z_dist = z_dist - len * nint(z_dist / len)
                
                r_dist = dsqrt(x_dist**2 + y_dist**2 + z_dist**2)

                if (r_dist <= r_cutoff) then
                    force = 4*epsilon*(12*sigma_12/(r_dist**13) - 6*sigma_6/(r_dist**7)) - force_cutoff
                    f_list(3*j-2) = f_list(3*j-2) + (x_dist/r_dist)*force
                    f_list(3*j-1) = f_list(3*j-1) + (y_dist/r_dist)*force
                    f_list(3*j)   = f_list(3*j)   + (z_dist/r_dist)*force
                end if

            end do
        end do

!!!================================= loops strart ===================================!!!

        count = 0; count_2 = 0; count_3 = 0
        correlation_list = 0.0d0
        !!! this is where all the iterations happen
        do i = 1, run_time

            if (mod(i,1000) == 0) print *, i

            do j = 1, 3*num_particles
                r_list(j) = r_list(j) + v_list(j)*dt + f_list(j)*dt2by2m
                r_list(j) = modulo(r_list(j), len)
            end do

            !!! this loop calculates force for the next time step
            if (mod(i,40) == 0) then
                call neighbours(num_particles, r_list, len, r_nearby, num_neighbour, neighbour_list)
                
                if (i >= 20000 .and. i <= 30000) then
                    count_3 = count_3 + 1
                    max_neighbours(1,count_3) = i
                    max_neighbours(2,count_3) = maxval(num_neighbour)
                end if
            end if

            do j = 1, num_particles 
                do k = 1, num_neighbour(j)

                    x_dist = r_list(3*j-2) - r_list(3*neighbour_list(j,k)-2)
                    y_dist = r_list(3*j-1) - r_list(3*neighbour_list(j,k)-1)
                    z_dist = r_list(3*j) - r_list(3*neighbour_list(j,k))
                    
                    x_dist = x_dist - len * nint(x_dist / len)
                    y_dist = y_dist - len * nint(y_dist / len)
                    z_dist = z_dist - len * nint(z_dist / len)
                    
                    r_dist = dsqrt(x_dist**2 + y_dist**2 + z_dist**2)

                    if (r_dist <= r_cutoff) then
                        force = 4*epsilon*(12*sigma_12/(r_dist**13) - 6*sigma_6/(r_dist**7)) - force_cutoff
                        f_new_list(3*j-2) = f_new_list(3*j-2) + (x_dist/r_dist)*force
                        f_new_list(3*j-1) = f_new_list(3*j-1) + (y_dist/r_dist)*force
                        f_new_list(3*j)   = f_new_list(3*j)   + (z_dist/r_dist)*force
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
            
            if (i > 20000 .and. mod(i,100) == 0) then
                count_2 = count_2 + 1
                call pair_correlation(num_particles, r_list, len, num_bins, correlation_list)
            end if 

            f_list = f_new_list
            f_new_list = 0.0d0

        end do

        do i = 1, num_bins
            correlation_list(2,i) = correlation_list(2,i)/count_2
        end do
        
        open(file = corr_path, unit = 1)
        do i = 1, num_bins
            write(1,*) correlation_list(1,i), correlation_list(2,i)
        end do
        close(1)

        open(file = neigh_path, unit = 2)
        do i = 1, size(max_neighbours(1,:))
            write(2,*) max_neighbours(1,i), max_neighbours(2,i)
        end do
        close(2)

    end subroutine gas_simulation_better


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


    subroutine neighbours(num_particles, r_list, len, r_nearby, num_neighbour, neighbour_list)

        integer :: num_particles
        integer :: j, k
        real*8 :: len, r_nearby
        real*8 :: r_list(3*num_particles)
        integer :: num_neighbour(num_particles)
        integer :: neighbour_list(num_particles,num_particles)
        real*8 :: x_dist, y_dist, z_dist, r_dist

        num_neighbour = 0; neighbour_list = 0

        do j = 1, num_particles - 1
            do k = j + 1, num_particles

                x_dist = r_list(3*j-2) - r_list(3*k-2)
                y_dist = r_list(3*j-1) - r_list(3*k-1)
                z_dist = r_list(3*j) - r_list(3*k)

                x_dist = x_dist - len * nint(x_dist / len)
                y_dist = y_dist - len * nint(y_dist / len)
                z_dist = z_dist - len * nint(z_dist / len)

                r_dist = dsqrt(x_dist**2 + y_dist**2 + z_dist**2)

                if (r_dist < r_nearby) then
                    num_neighbour(j) = num_neighbour(j) + 1
                    num_neighbour(k) = num_neighbour(k) + 1
                    neighbour_list(j,num_neighbour(j)) = k
                    neighbour_list(k,num_neighbour(k)) = j
                end if

            end do
        end do

    end subroutine neighbours


    subroutine pair_correlation(num_particles, r_list, len, num_bins, correlation_list)
    
        integer :: num_particles
        integer :: j, k
        real*8 :: len, x_dist, y_dist, z_dist, r_dist
        real*8 :: pi
        integer :: bin, num_bins
        real*8 ::  dr , r_max, vol, rho, shell_vol
        real*8 :: g_of_r(num_bins), r_list(3*num_particles), correlation_list(2,num_bins)
        integer :: hist(num_bins)
    
        
        dr = 0.1d0
        r_max = 0.5d0*len
        num_bins = int(r_max/dr)
        pi = 3.14159265359d0

        g_of_r = 0.0d0
        hist = 0
    
        do j = 1, num_particles - 1
            do k = j + 1, num_particles

                x_dist = r_list(3*j-2) - r_list(3*k-2)
                y_dist = r_list(3*j-1) - r_list(3*k-1)
                z_dist = r_list(3*j) - r_list(3*k)

                x_dist = x_dist - len * nint(x_dist / len)
                y_dist = y_dist - len * nint(y_dist / len)
                z_dist = z_dist - len * nint(z_dist / len)

                r_dist = dsqrt(x_dist**2 + y_dist**2 + z_dist**2)
                
                if (r_dist < r_max) then
                    bin = int(r_dist/dr) + 1      
                    if (bin <= num_bins) hist(bin) = hist(bin) + 2   
                end if
            end do
        end do
    
        vol = len**3
        rho = dfloat(num_particles)/vol
    
        ! normalizing g(r_dist)
        do j = 1, num_bins
            r_dist = (j - 0.5d0)*dr
            shell_vol = 4.0d0/3.0d0 * pi * ((r_dist + dr)**3 - (r_dist)**3)
            g_of_r(j) = dfloat(hist(j)) / (rho * shell_vol * dfloat(num_particles))

            correlation_list(1,j) = r_dist
            correlation_list(2,j) = correlation_list(2,j) + g_of_r(j)
        end do        

    end subroutine pair_correlation


end module gas_sim_better