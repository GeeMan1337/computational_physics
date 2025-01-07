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


program test_1
    implicit none
    integer :: len
    integer:: seed(8)
    real :: output(10)

    call random_seed(get=seed)
    len=size(seed)
    
    !call random_seed(put=seed)
    seed=[22342323,32313,34234234,342423424,454353534,78686786,65464564,53452452]

    call random_number_seed(seed,len-2, output,10)
    
    print *, output

    end program test_1