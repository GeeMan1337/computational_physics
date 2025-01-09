program large_number_trials
    implicit none

    integer :: i, j        !iterating variable
    real *8 :: sum_list(10000), rand_num(10000)
    real *8 :: sum_bigger_list(100000), rand_bigger_num(100000)

    integer :: num_bins, min_val, max_val
    real *8, allocatable :: bins(:), frequency(:)
    real *8 :: bin_size 

    !!! question 1h
    !!! finding the sum of 10^4 trials 10^4 times
    do i=1,10000
        call random_seed()
        do j=1,10000
            call random_number(rand_num(j))
        end do
        sum_list(i)=sum(rand_num)
    end do

    !!! code to construct bins and frequency for histogram
    min_val=int(minval(sum_list)-bin_size)
    max_val=int(maxval(sum_list)+1)

    if (mod(real(max_val),bin_size) /= 0.0) then
        max_val=int(max_val-mod(real(max_val),bin_size)+2*bin_size)
    end if

    num_bins=int(max_val-min_val)

    do i=1,num_bins
        bins(i)=min_val+(i-1)*bin_size
    end do

    print *, bins

end program large_number_trials