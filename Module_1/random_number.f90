program random_number_generation
    implicit none

    integer :: i, j     !these are variables for loops
    integer :: seed_size
    integer, allocatable :: get_seed(:), give_seed(:)
    integer :: date_time(8)
    real *8 :: rand_num(10), rand_num_10(1:10,1:10)
    call random_seed(size = seed_size)
    allocate(get_seed(seed_size))

    !!! can use this to give custom seed

    !allocate(give_seed(seed_size))
    !call date_and_time(values=date_time)
    !give_seed=date_time     !using date and time as the custom seed
    !call random_seed(put=give_seed)

    !!! using the same seed to generate 10 random numbers
    call random_seed()
    do i=1,10
        call random_number(rand_num(i))    
    end do

    !!! writing in the file
    open(access="sequential", unit=1, file="E:\computational_physics\Module_1_out\test_ran.dat")        !access="sequential" is default and erases everything before writing
    write (1,"(a,/,a)") "Anirban Nath's assignment. Reg. no. 20242019", " "     !used a line break

    do i=1,10
        write (1,*) rand_num(i)
    end do

    write (1,"(a,/,a,/,a)") "", "Changing seed and generating 10 new random numbers:", ""

    !!! changing seed and writing 10 new random numbers
    call random_seed()

    do i=1,10
        call random_number(rand_num(i))
    end do

    do i=1,10
        write (1,*) rand_num(i)
    end do
    
    close(1)

    !!! changing seed 10 times and writing 10 random numbers for each seed   
    do j=1,10
        call random_seed()
        do i=1,10
            call random_number(rand_num_10(i,j))
        end do
    end do

    open(unit=2, file="E:\computational_physics\Module_1_out\test_ran_10_seeds.dat")

    do i = 1, 10
        write(2, "(10f15.10)") rand_num_10(i,:)
    end do
    
    close(2)

end program random_number_generation