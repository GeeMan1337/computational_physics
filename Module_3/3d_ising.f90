program question_1
    implicit none

    !real :: j_ising, kt, start_config
    !integer :: len, niter, unit_num
    !character(len=*) :: path

    call ising_model(20, 1.0, 60.0, 10**5, "minus", 1, "E:\computational_physics\Module_3_out\question_1_data.dat")
    call ising_model(10, 1.0, 1.0, 10**5, "plus", 2, "E:\computational_physics\Module_3_out\question_2_data.dat")
    call ising_model(10, 1.0, 0.001, 10**5, "plus", 3, "E:\computational_physics\Module_3_out\question_3_data.dat")
    call ising_model(10, 1.0, 60.0, 10**5, "random", 4, "E:\computational_physics\Module_3_out\question_4_data.dat")
    call ising_model(10, 1.0, 0.001, 10**5, "random", 5, "E:\computational_physics\Module_3_out\question_5_data.dat")


end program question_1


subroutine ising_model(len, j_ising, kt, niter, start_config, unit_num, path)
    implicit none 

    integer :: i, j, k      !iterating variables
    real :: temp, temp_2, temp_3, temp_4, temp_5
    real :: a, b, c, d, e, f      !a, b, c... are neighbouring spin variables
    real :: avg_energy, avg_magnet, heat_cap, chi       !Macroscopic observables
    real :: energy, magnet, energy_per_spin, magnet_per_spin, config_energy !magnet is the total magnetic moment of the entire lattice
    real :: scaled_energy, scaled_magnet
    integer :: len, niter, num_spin, index, unit_num
    real :: j_ising, kt     !kt is temperature multiplied by boltzmann const. 
    real :: lattice_energy      !these are function variables
    integer :: spin(len, len, len)       !3D lattice with spins
    real :: energy_list(niter*(len**3)+1), magnet_list(niter*(len**3)+1)
    character(len=*) :: start_config, path

    num_spin = len**3
    index = 1
    magnet = 0

    if (start_config == "plus") then        !this is initializing the spins in the lattice
        do i = 1, len 
            do j = 1, len 
                do k = 1, len
                    spin(i,j,k) = 1
                    magnet = magnet + spin(i,j,k)      !this is adding up the magnetic moments
                end do
            end do
        end do

    elseif (start_config == "minus") then
        do i = 1, len 
            do j = 1, len 
                do k = 1, len
                    spin(i,j,k) = -1
                    magnet = magnet + spin(i,j,k)      !this is adding up the magnetic moments
                end do
            end do
        end do

    elseif (start_config == "random") then      !we do this to initialize random spins
        do i = 1, len 
            do j = 1, len 
                do k = 1, len
                    call random_number(temp)
                    if (temp<0.5) then
                        spin(i,j,k) = -1
                    else
                        spin(i,j,k) = 1
                    end if
                    magnet = magnet + spin(i,j,k)      !this is adding up the magnetic moments
                end do
            end do
        end do

    end if

    energy = lattice_energy(spin,len,1.0)
    magnet_list(index) = magnet
    energy_list(index) = energy


    !now we flip spins and keep the spin, or discard it with a probability ~ exp(-x)

    do i = 1, niter !we do niter number of monte carlo steps (a single step is explained in the comment below)
        index = index + 1
        do j = 1, num_spin  !when all the j's are done, it amounts to one monte carlo step
            !index = index + 1

            call random_number(temp_2)
            call random_number(temp_3)
            call random_number(temp_4)
            temp_2 = int(temp_2*len + 1); temp_3 = int(temp_3*len + 1); temp_4 = int(temp_4*len + 1)

            a=temp_2-1; b=temp_2+1; c=temp_3-1; d=temp_3+1; e=temp_4-1; f=temp_4+1 !neighbouring spin indices
            if (temp_2==1) a=len     !these impose periodic boundary conditions
            if (temp_3==1) c=len
            if (temp_4==1) e=len
            if (temp_2==len) b=1
            if (temp_3==len) d=1
            if (temp_4==len) f=1

            temp = 0
            temp = temp - j_ising*spin(int(temp_2),int(temp_3),int(temp_4))* &  !temp is the initial energy before flip
            (spin(int(a),int(temp_3),int(temp_4)) + spin(int(b),int(temp_3),int(temp_4)) &
            + spin(int(temp_2),int(c),int(temp_4)) + spin(int(temp_2),int(d),int(temp_4)) &
            + spin(int(temp_2),int(temp_3),int(e)) + spin(int(temp_2),int(temp_3),int(f)))
            
            !flipping a random spin
            spin(int(temp_2),int(temp_3),int(temp_4)) = -spin(int(temp_2),int(temp_3),int(temp_4))
            
            config_energy = 0
            config_energy = config_energy - j_ising*spin(int(temp_2),int(temp_3),int(temp_4))* &
            (spin(int(a),int(temp_3),int(temp_4)) + spin(int(b),int(temp_3),int(temp_4)) &
            + spin(int(temp_2),int(c),int(temp_4)) + spin(int(temp_2),int(d),int(temp_4)) &
            + spin(int(temp_2),int(temp_3),int(e)) + spin(int(temp_2),int(temp_3),int(f)))
            
            if (config_energy - temp <= 0) then !this condition is energetically favourabale 
                energy = energy + (config_energy - temp)
                magnet = magnet + (2*spin(int(temp_2),int(temp_3),int(temp_4)))
                !magnet_list(index) = magnet
                !energy_list(index) = energy
            else
                call random_number(temp_5)  !here final energy is more than the initial
                if (exp(-(config_energy-temp)/kt) > temp_5) then !we only accept the spin change with some probability
                    energy = energy + (config_energy - temp)
                    magnet = magnet + (2*spin(int(temp_2),int(temp_3),int(temp_4)))
                    !magnet_list(index) = magnet
                    !energy_list(index) = energy
                else    !we revert to the initial spin if we fail this test
                    spin(int(temp_2),int(temp_3),int(temp_4)) = -spin(int(temp_2),int(temp_3),int(temp_4))
                end if
            end if
        end do
        magnet_list(index) = magnet
        energy_list(index) = energy
    end do

    open(unit=unit_num, file=path)
    do i = 1, index
        write(unit_num,*) i, " - ", magnet_list(i), " - ", magnet_list(i)/num_spin, " - ", &
        energy_list(i), " - ", energy_list(i)/num_spin
    end do
    close(unit_num)

end subroutine ising_model


function lattice_energy(lattice_spin, length, j_ising)  !this calculates the energy of the entire lattice 
    implicit none

    integer :: i, j, k, a, b, c, d, e, f
    integer, intent(in) :: length       !this is the lattice length in one direction
    integer, intent(in) :: lattice_spin(length, length, length)
    real, intent(in) ::j_ising
    real :: lattice_energy

    lattice_energy = 0
    
    do i = 1, length
        do j = 1, length 
            do k = 1, length
                a=i-1; b=i+1; c=j-1; d=j+1; e=k-1; f=k+1
                if (i==1) a=length     !these impose periodic boundary conditions
                if (j==1) c=length
                if (k==1) e=length
                if (i==length) b=1
                if (j==length) d=1
                if (k==length) f=1

                lattice_energy = lattice_energy - j_ising*lattice_spin(i,j,k)* &
                (lattice_spin(int(a),j,k) + lattice_spin(int(b),j,k) &
                + lattice_spin(i,int(c),k) + lattice_spin(i,int(d),k) &
                + lattice_spin(i,j,int(e)) + lattice_spin(i,j,int(f)))
            end do
        end do
    end do

    lattice_energy = lattice_energy/2   !dividing by two to compensate for overcounting

end function lattice_energy
