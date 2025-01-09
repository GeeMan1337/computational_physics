program large_number_trials
    implicit none

    integer :: k, l        !iterating variable
    real *8 :: sum_list(10000), rand_num(10000)
    real *8 :: sum_bigger_list(100000), rand_bigger_num(100000)

    !integer :: num_bins, min_val, max_val
    !real, allocatable :: bins(:), frequency(:)
    !real, parameter :: bin_size=2

    !!! question 1h
    !!! finding the sum of 10^4 trials 10^4 times
    do k=1,10000
        call random_seed()
        do l=1,10000
            call random_number(rand_num(l))
        end do
        sum_list(k)=sum(rand_num)
    end do

    call graph_data(sum_list,10000,0.5,"file path here",5)
    
end program large_number_trials

!!! code to construct bins and frequency for histogram
subroutine graph_data(data_array, data_size, bin_size, path, unit_num)
    implicit none

    character(len=*), intent(in) :: path
    integer :: i, j
    integer :: num_bins, min_val, max_val
    integer, intent(in) :: data_size, unit_num
    real, allocatable :: bins(:), frequency(:)
    real, intent(in) :: bin_size
    real *8, intent(in) :: data_array(data_size)

    min_val=int(minval(data_array)-bin_size)
    max_val=int(maxval(data_array)+1.0)

    if (mod(real(max_val-min_val),bin_size) /= 0.0) then
        max_val=int(real(max_val)-mod(real(max_val-min_val),bin_size)+bin_size)
    end if

    num_bins=int(real(max_val-min_val)/bin_size)

    allocate(bins(num_bins))
    allocate(frequency(num_bins))

    do i=1,num_bins
        bins(i)=real(min_val)+(i-1)*bin_size
        frequency(i)=0
    end do

    do i=1,num_bins
        do j=1,10000
            if (data_array(j)>=bins(i) .and. data_array(j)<bins(i)+bin_size) then
                frequency(i)=frequency(i)+1
            end if
        end do
    end do
 
    do i=1,num_bins
        bins(i)=bins(i)+bin_size/2
    end do

    open(unit=unit_num, file=path)
    
    do i=1,num_bins
        write (unit_num,*) bins(i), " - ", frequency(i), " - ", frequency(i)/data_size
    end do

    close(unit_num)

end subroutine graph_data