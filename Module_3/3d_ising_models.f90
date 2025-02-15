program questions
    implicit none

    !real :: j_ising, kt, start_config
    !integer :: len, niter, unit_num
    !character(len=*) :: path

    call ising_model(5, 1.0, 0.1, 10**5, "random", 1, "E:\computational_physics\Module_3_out\question_1_data.dat")
    call ising_model(5, 1.0, 0.1, 10**5, "random", 2, "E:\computational_physics\Module_3_out\question_2_data.dat")
    call ising_model(5, 1.0, 0.1, 10**5, "random", 3, "E:\computational_physics\Module_3_out\question_3_data.dat")

    !call ising_model_2(10, 1.0, 9.0, 10**5, "plus", 4, "E:\computational_physics\Module_3_out\question_4_data.dat")
    !call ising_model_2(10, 1.0, 3.0, 10**5, "plus", 5, "E:\computational_physics\Module_3_out\question_5_data.dat")
    !call ising_model_2(10, 1.0, 0.1, 10**5, "plus", 6, "E:\computational_physics\Module_3_out\question_6_data.dat")

    !call ising_model_3(10, 1.0, 9.0, 10**5, "plus", 7, "E:\computational_physics\Module_3_out\question_7_data.dat")
    !call ising_model_3(10, 1.0, 3.0, 10**5, "plus", 8, "E:\computational_physics\Module_3_out\question_8_data.dat")
    !call ising_model_3(10, 1.0, 0.1, 10**5, "plus", 9, "E:\computational_physics\Module_3_out\question_9_data.dat")

end program questions


!this model uses only the nearest neighbours
subroutine ising_model(len, j_ising, kt, niter, start_config, unit_num, path)
    implicit none 

    integer :: i, j, k      !iterating variables
    real :: temp, temp_2, temp_3, temp_4, temp_5
    real :: avg_energy, avg_magnet, heat_cap, chi       !Macroscopic observables
    real :: energy, magnet, energy_per_spin, magnet_per_spin, config_energy !magnet is the total magnetic moment of the entire lattice
    real :: scaled_energy, scaled_magnet
    integer :: len, niter, num_spin, index, unit_num
    real :: j_ising, kt     !kt is temperature multiplied by boltzmann const. 
    integer :: spin(len, len, len)       !3D lattice with spins
    real :: energy_list(niter*(len**3)+1), magnet_list(niter*(len**3)+1)
    character(len=*) :: start_config, path

    real :: lattice_energy      !these are function variables
    integer :: sum_1_neighbour

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

    energy = lattice_energy(spin,len,j_ising)
    magnet_list(index) = magnet
    energy_list(index) = energy
print*, energy/num_spin

    !now we flip spins and keep the spin, or discard it with a probability ~ exp(-x)

    do i = 1, niter !we do niter number of monte carlo steps (a single step is explained in the comment below)
        index = index + 1
        do j = 1, num_spin  !when all the j's are done, it amounts to one monte carlo step

            call random_number(temp_2)
            call random_number(temp_3)
            call random_number(temp_4)
            temp_2 = temp_2*len + 1; temp_3 = temp_3*len + 1; temp_4 = temp_4*len + 1

            temp = 0        !temp is the initial energy before flip
            temp = temp &
            - j_ising*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_1_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4)) 
            
            !flipping a random spin
            spin(int(temp_2),int(temp_3),int(temp_4)) = -spin(int(temp_2),int(temp_3),int(temp_4))
            
            config_energy = 0
            config_energy = config_energy &
            - j_ising*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_1_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4))

            if (config_energy - temp <= 0) then !this condition is energetically favourabale 
                energy = energy + (config_energy - temp)
                magnet = magnet + (2*spin(int(temp_2),int(temp_3),int(temp_4)))
                
            else
                call random_number(temp_5)  !here final energy is more than the initial
                if (exp(-(config_energy-temp)/kt) > temp_5) then !we only accept the spin change with some probability
                    energy = energy + (config_energy - temp)
                    magnet = magnet + (2*spin(int(temp_2),int(temp_3),int(temp_4)))
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


!this model has upto second nearest neighbours
subroutine ising_model_2(len, j_ising, kt, niter, start_config, unit_num, path)
    implicit none 

    integer :: i, j, k      !iterating variables
    real :: temp, temp_2, temp_3, temp_4, temp_5
    real :: avg_energy, avg_magnet, heat_cap, chi       !Macroscopic observables
    real :: energy, magnet, energy_per_spin, magnet_per_spin, config_energy !magnet is the total magnetic moment of the entire lattice
    real :: scaled_energy, scaled_magnet
    integer :: len, niter, num_spin, index, unit_num
    real :: j_ising, kt     !kt is temperature multiplied by boltzmann const. 
    integer :: spin(len, len, len)       !3D lattice with spins
    real :: energy_list(niter*(len**3)+1), magnet_list(niter*(len**3)+1)
    character(len=*) :: start_config, path

    real :: lattice_energy_2      !these are function variables
    integer :: sum_1_neighbour, sum_2_neighbour

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

    energy = lattice_energy_2(spin,len,1.0)
    magnet_list(index) = magnet
    energy_list(index) = energy
print*, energy/num_spin

    !now we flip spins and keep the spin, or discard it with a probability ~ exp(-x)

    do i = 1, niter !we do niter number of monte carlo steps (a single step is explained in the comment below)
        index = index + 1
        call random_seed()
        do j = 1, num_spin  !when all the j's are done, it amounts to one monte carlo step
            !index = index + 1

            call random_number(temp_2)
            call random_number(temp_3)
            call random_number(temp_4)
            temp_2 = int(temp_2*len + 1); temp_3 = int(temp_3*len + 1); temp_4 = int(temp_4*len + 1)

            temp = 0        !temp is the energy before flip
            temp = temp &
            - j_ising*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_1_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4)) &
            - j_ising*(1.0/(2.0**1.5))*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_2_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4))

            
            !flipping a random spin
            spin(int(temp_2),int(temp_3),int(temp_4)) = -spin(int(temp_2),int(temp_3),int(temp_4))
            
            config_energy = 0       !this is energy after flip
            config_energy = config_energy &
            - j_ising*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_1_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4)) &
            - j_ising*(1.0/(2.0**1.5))*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_2_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4))
            
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

end subroutine ising_model_2


!this model has upto third nearest neighbours
subroutine ising_model_3(len, j_ising, kt, niter, start_config, unit_num, path)
    implicit none 

    integer :: i, j, k      !iterating variables
    real :: temp, temp_2, temp_3, temp_4, temp_5
    real :: avg_energy, avg_magnet, heat_cap, chi       !Macroscopic observables
    real :: energy, magnet, energy_per_spin, magnet_per_spin, config_energy !magnet is the total magnetic moment of the entire lattice
    real :: scaled_energy, scaled_magnet
    integer :: len, niter, num_spin, index, unit_num
    real :: j_ising, kt     !kt is temperature multiplied by boltzmann const. 
    integer :: spin(len, len, len)       !3D lattice with spins
    real :: energy_list(niter*(len**3)+1), magnet_list(niter*(len**3)+1)
    character(len=*) :: start_config, path

    real :: lattice_energy_3      !these are function variables
    integer :: sum_1_neighbour, sum_2_neighbour, sum_3_neighbour

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

    energy = lattice_energy_3(spin,len,1.0)
    magnet_list(index) = magnet
    energy_list(index) = energy
print*, energy/num_spin

    !now we flip spins and keep the spin, or discard it with a probability ~ exp(-x)

    do i = 1, niter !we do niter number of monte carlo steps (a single step is explained in the comment below)
        index = index + 1
        do j = 1, num_spin  !when all the j's are done, it amounts to one monte carlo step
            !index = index + 1

            call random_number(temp_2)
            call random_number(temp_3)
            call random_number(temp_4)
            temp_2 = int(temp_2*len + 1); temp_3 = int(temp_3*len + 1); temp_4 = int(temp_4*len + 1)

            temp = 0        !temp is the initial energy before flip
            temp = temp &
            - j_ising*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_1_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4)) &
            - j_ising*(1.0/(2.0**1.5))*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_2_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4)) &
            - j_ising*(1.0/(3.0**1.5))*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_3_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4))

            
            !flipping a random spin
            spin(int(temp_2),int(temp_3),int(temp_4)) = -spin(int(temp_2),int(temp_3),int(temp_4))
            
            config_energy = 0
            config_energy = config_energy &
            - j_ising*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_1_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4)) &
            - j_ising*(1.0/(2.0**1.5))*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_2_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4)) &
            - j_ising*(1.0/(3.0**1.5))*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_3_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4))
            
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

end subroutine ising_model_3


function lattice_energy(lattice_spin, length, j_ising)  !this calculates the energy of the entire lattice 
    implicit none

    integer :: i, j, k
    integer, intent(in) :: length       !this is the lattice length in one direction
    integer, intent(in) :: lattice_spin(length, length, length)
    real, intent(in) ::j_ising
    real :: lattice_energy
    integer :: sum_1_neighbour      !these are function variables

    lattice_energy = 0
    
    do i = 1, length
        do j = 1, length 
            do k = 1, length

                lattice_energy = lattice_energy &
                - j_ising*lattice_spin(i,j,k)*sum_1_neighbour(length,lattice_spin,i,j,k)

            end do
        end do
    end do

    lattice_energy = lattice_energy/2   !dividing by two to compensate for overcounting

end function lattice_energy


function lattice_energy_2(lattice_spin, length, j_ising)  !this calculates the energy of the entire lattice 
    implicit none

    integer :: i, j, k
    integer, intent(in) :: length       !this is the lattice length in one direction
    integer, intent(in) :: lattice_spin(length, length, length)
    real, intent(in) ::j_ising
    real :: lattice_energy_2

    integer :: sum_1_neighbour, sum_2_neighbour      !these are function variables

    lattice_energy_2 = 0
    
    do i = 1, length
        do j = 1, length 
            do k = 1, length
                
                lattice_energy_2 = lattice_energy_2 & 
                - j_ising*lattice_spin(i,j,k)*sum_1_neighbour(length,lattice_spin,i,j,k) &
                - j_ising*(1.0/(2.0**1.5))*(lattice_spin(i,j,k)*sum_2_neighbour(length,lattice_spin,i,j,k))

            end do
        end do
    end do

    lattice_energy_2 = lattice_energy_2/2   !dividing by two to compensate for overcounting

end function lattice_energy_2


function lattice_energy_3(lattice_spin, length, j_ising)  !this calculates the energy of the entire lattice 
    implicit none

    integer :: i, j, k
    integer, intent(in) :: length       !this is the lattice length in one direction
    integer, intent(in) :: lattice_spin(length, length, length)
    real, intent(in) ::j_ising
    real :: lattice_energy_3

    integer :: sum_1_neighbour, sum_2_neighbour      !these are function variables

    lattice_energy_3 = 0
    
    do i = 1, length
        do j = 1, length 
            do k = 1, length
                
                lattice_energy_3 = lattice_energy_3 & 
                - j_ising*lattice_spin(i,j,k)*sum_1_neighbour(length,lattice_spin,i,j,k) &
                - j_ising*(1.0/(2.0**1.5))*(lattice_spin(i,j,k)*sum_2_neighbour(length,lattice_spin,i,j,k)) &
                - j_ising*(1.0/(3.0**1.5))*(lattice_spin(i,j,k)*sum_2_neighbour(length,lattice_spin,i,j,k))

            end do
        end do
    end do

    lattice_energy_3 = lattice_energy_3/2   !dividing by two to compensate for overcounting

end function lattice_energy_3


function sum_1_neighbour(length, lattice_spin, index_1, index_2, index_3)
    implicit none

    integer :: a, b, c, d, e, f      !a, b, c... are neighbouring spin variables
    integer :: sum_1_neighbour
    integer :: index_1, index_2, index_3, length, lattice_spin(length, length, length)

    a=index_1-1; b=index_1+1; c=index_2-1; d=index_2+1; e=index_3-1; f=index_3+1 !neighbouring spin indices

    if (index_1==1) a=length     !these impose periodic boundary conditions
    if (index_2==1) c=length
    if (index_3==1) e=length
    if (index_1==length) b=1
    if (index_2==length) d=1
    if (index_3==length) f=1

    sum_1_neighbour = 0
    sum_1_neighbour = lattice_spin(a,index_2,index_3) + lattice_spin(b,index_2,index_3) &
                      + lattice_spin(index_1,c,index_3) + lattice_spin(index_1,d,index_3) &
                      + lattice_spin(index_1,index_2,e) + lattice_spin(index_1,index_2,f)

end function sum_1_neighbour


function sum_2_neighbour(length, lattice_spin, index_1, index_2, index_3)
    implicit none

    integer :: a, b, c, d, e, f      !a, b, c... are neighbouring spin variables
    integer :: sum_2_neighbour
    integer :: index_1, index_2, index_3, length, lattice_spin(length, length, length)

    a=index_1-1; b=index_1+1; c=index_2-1; d=index_2+1; e=index_3-1; f=index_3+1 !neighbouring spin indices

    if (index_1==1) a=length     !these impose periodic boundary conditions
    if (index_2==1) c=length
    if (index_3==1) e=length
    if (index_1==length) b=1
    if (index_2==length) d=1
    if (index_3==length) f=1

    sum_2_neighbour = 0
    sum_2_neighbour = lattice_spin(a,c,index_3) + lattice_spin(a,d,index_3) &
        + lattice_spin(a,index_2,e) + lattice_spin(a,index_2,f) &
        + lattice_spin(b,c,index_3) + lattice_spin(b,d,index_3) &
        + lattice_spin(b,index_2,e) + lattice_spin(b,index_2,f) &
        + lattice_spin(index_1,c,e) + lattice_spin(index_1,c,f) &
        + lattice_spin(index_1,d,e) + lattice_spin(index_1,d,f)

end function sum_2_neighbour


function sum_3_neighbour(length, lattice_spin, index_1, index_2, index_3)
    implicit none

    integer :: a, b, c, d, e, f      !a, b, c... are neighbouring spin variables
    integer :: sum_3_neighbour
    integer :: index_1, index_2, index_3, length, lattice_spin(length, length, length)

    a=index_1-1; b=index_1+1; c=index_2-1; d=index_2+1; e=index_3-1; f=index_3+1 !neighbouring spin indices

    if (index_1==1) a=length     !these impose periodic boundary conditions
    if (index_2==1) c=length
    if (index_3==1) e=length
    if (index_1==length) b=1
    if (index_2==length) d=1
    if (index_3==length) f=1

    sum_3_neighbour = 0
    sum_3_neighbour = lattice_spin(a,c,e) + lattice_spin(a,c,f) &
        + lattice_spin(a,d,e) + lattice_spin(a,d,f) &
        + lattice_spin(b,c,e) + lattice_spin(b,c,f) &
        + lattice_spin(b,d,e) + lattice_spin(b,d,f)

end function sum_3_neighbour