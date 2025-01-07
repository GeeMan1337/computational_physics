program random_number_generation
    implicit none

    integer :: i, j, k        !these are variables for loops
    integer :: seed_size
    integer, allocatable :: get_seed(:), give_seed(:)
    integer :: date_time(8)
    real *8 :: rand_num(10)

    call random_seed(size = seed_size)
    allocate(get_seed(seed_size))

    !!! can use this to give custom seed

    !allocate(give_seed(seed_size))
    !call date_and_time(values=date_time)
    !give_seed=date_time     !using date and time as the custom seed
    !call random_seed(put=give_seed)

    !!! using the same seed to generate 10 random numbers
    do i=1,10
        call random_number(rand_num(i))    
    end do

    !!! writing in the file
    open(access="sequential", unit=1, file="E:\computational_physics\Module_1_out\test_ran.dat")        !access="sequential" is default and erases everything before writing
    write (1,"(a,/,a)") "Anirban Nath's assignment. Reg. no. 20242019", ""     !used a line break

    do i=1,10
        write (1,*) rand_num(i)
    end do

    write (1,"(a,/,a,/,a)") "", "Changing seed and generating 10 new random numbers:", ""

    !!! changing seed and writing 10 new random numbers
    do i=1,10
        call random_seed(get=get_seed)
        call random_number(rand_num(i))
    end do

    do i=1,10
        write (1,*) rand_num(i)
    end do

    !!! changing seed 10 times and writing 10 random numbers for each seed
        
  !  integer, allocatable :: get_seed_10(:), give_seed_10(:), rand_num_10(10,10)

 !   call random_seed(size=seed_size)
  ! allocate(get_seed_10(seed_size))
 !   allocate(give_seed_10(seed_size))

  !  do i=1,10
  !      call random_seed(get=get_seed_10)
  !      give_seed_10=get_seed_10
    !    do j=1,10
   !         call random_seed(put=give_seed_10)
    !        call random_number(rand_num(j))
   !     end do
  !  end do



end program random_number_generation


!!! subroutine to make random numbers. requires 4 inputs

subroutine random_number_seed(input_seed, seed_len, out_num, out_len)

    implicit none
    
    integer :: seed_len, out_len, i, min_size
    integer, intent(in) :: input_seed(seed_len)
    real, intent(out) :: out_num(out_len)

    if (seed_len /= size(input_seed) .and. out_len /= size(out_num)) then
        print *, "Error: Subroutine random_number_seed error. Lengths don't match"
        go to 1
    end if

    if (seed_len<8) then
        print *, "Error: Subroutine random_number_seed error. Seed must be an array with atleast 8 elements."
        go to 1
    end if

    min_size=min(seed_len,out_len)

    call random_seed(size=min_size)
    call random_seed(put=input_seed)

    do i=1,out_len
        call random_number(out_num(i))
    end do

1   end subroutine random_number_seed