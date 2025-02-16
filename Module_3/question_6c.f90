program questions
    implicit none

    !!! avg_magnet, avg_magnet_abs, chi, chi_abs, avg_energy, heat_cap, binders_cum
    integer :: i
    real :: answers(7)
    real :: temp, temp_2, temp_3
    integer :: iterations

    !L = 10 -------------------------------------------
    temp = 0; temp_2 = 0; temp_3 = 0
    iterations = 1
    
    do i = 1, iterations
        call ising_model(10, 1.0, 3.9, 50000, "random", answers, 1, "dont write")
        temp = temp + answers(1)/1000
        temp_2 = temp_2 + answers(2)/1000
        temp_3 = temp_3 + answers(5)/1000
    end do

    temp = temp/iterations
    temp_2 = temp_2/iterations
    temp_3 = temp_3/iterations

    print *, "For L = 10"
    print *, "The instantaneous magnetisation per spin fluctuates around the value:", temp
    print *, "The instantaneous magnetisation (abs value) per spin fluctuates around the value:", temp_2
    print *, "The instantaneous energy per spin fluctuates around the value:", temp_3
    print *, " "
    call ising_model(10, 1.0, 3.9, 50000, "random", answers, 1, "E:\computational_physics\Module_3_out\question_6c_data.dat")

end program questions


!this model uses only the nearest neighbours
subroutine ising_model(len, j_ising, kt, niter, start_config, return_values, unit_num, path)
    implicit none 

    integer :: i, j, k      !iterating variables
    real :: temp, temp_2, temp_3, temp_4, temp_5
    real :: avg_magnet, avg_magnet_abs, chi, chi_abs, avg_energy, heat_cap, binders_cum !Macroscopic observables
    real :: energy, magnet, config_energy !magnet is the total magnetic moment of the entire lattice
    integer :: len, niter, num_spin, index, unit_num
    real :: j_ising, kt     !kt is temperature multiplied by boltzmann const. 
    integer :: spin(len, len, len)       !3D lattice with spins
    real :: energy_list(niter*(len**3)+1), magnet_list(niter*(len**3)+1)
    character(len=*) :: start_config, path
    integer :: relax_time, nstat
    real :: return_values(7)

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

    !now we flip spins and keep the spin, or discard it with a probability ~ exp(-x)

    do i = 1, niter !we do niter number of monte carlo steps (a single step is explained in the comment below)
        index = index + 1
        do j = 1, num_spin  !when all the j's are done, it amounts to one monte carlo step

            call random_number(temp_2)
            call random_number(temp_3)
            call random_number(temp_4)
            temp_2 = temp_2*len + 1; temp_3 = temp_3*len + 1; temp_4 = temp_4*len + 1

            !temp is the initial energy before flip
            temp = - j_ising*spin(int(temp_2),int(temp_3),int(temp_4)) &
            *sum_1_neighbour(len,spin,int(temp_2),int(temp_3),int(temp_4)) 
            
            !flipping a random spin
            spin(int(temp_2),int(temp_3),int(temp_4)) = -spin(int(temp_2),int(temp_3),int(temp_4))
 
            config_energy = - j_ising*spin(int(temp_2),int(temp_3),int(temp_4)) &
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

    if (path == "dont write" .or. path == "don't write") then
        goto 7

    else
        open(unit=unit_num, file=path)
        do i = 1, index
            write(unit_num,*) i, " , ", magnet_list(i), " , ", magnet_list(i)/num_spin, " , ", &
            energy_list(i), " , ", energy_list(i)/num_spin
        end do
        close(unit_num)

    end if

    !this section calculates the averages for the entire lattice
7    temp = 0; temp_2 = 0; temp_3 = 0; temp_4 = 0
    avg_energy = 0; avg_magnet = 0; avg_magnet_abs = 0
    nstat = 10000
    relax_time = 1

    do i = nstat, index, relax_time
        temp = temp + 1
        avg_energy = avg_energy + energy_list(i)        !this will be <E>
        avg_magnet = avg_magnet + magnet_list(i)        !this will be <M>
        avg_magnet_abs = avg_magnet_abs + abs(magnet_list(i))   !this will be <abs(M)> 
        temp_2 = temp_2 + energy_list(i)**2     !this will calculate <E^2>
        temp_3 = temp_3 + magnet_list(i)**2     !this will calculate <M^2>
        temp_4 = temp_4 + magnet_list(i)**4     !this will calculate <M^4>
    end do

    avg_energy = avg_energy/temp
    avg_magnet = avg_magnet/temp
    avg_magnet_abs = avg_magnet_abs/temp
    chi = (temp_3/temp - avg_magnet**2)/kt
    chi_abs = (temp_3/temp - avg_magnet_abs**2)/kt
    heat_cap = (temp_2/temp - avg_energy**2)/(kt**2)
    binders_cum = (temp_4/temp)/(3*(temp_3/temp)**2)

    return_values = [avg_magnet, avg_magnet_abs, chi, chi_abs, avg_energy, heat_cap, binders_cum]

end subroutine ising_model


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