!!! question 2
program random_numbers 
    implicit none

    integer :: i, k, num 
    real, allocatable :: rand_num(:), rand_num_square(:), corr_list(:)
    real :: sigma, avg, sum_val

    num = 10**5

    allocate(rand_num(num))
    allocate(rand_num_square(num))
    allocate(corr_list(num))

    do i = 1, num
        call random_seed()
        call random_number(rand_num(i))
    end do

    !!! list with square of the random numbers
    do i = 1, num
        rand_num_square(i) = rand_num(i)**2
    end do

    !!! this is standard deviation
    avg = (sum(rand_num)/num)
    sigma = sqrt(sum(rand_num_square)/num - avg**2)

    print *, "The average of", num, "random numbers is ", avg 
    print *, "The standard deviation of", num, "random numbers is ", sigma 

    call graph_data(rand_num,num,0.01,"even","E:\computational_physics\Module_2_out\question_2_a_data.dat",1)

    open(unit = 2, file = "E:\computational_physics\Module_2_out\question_2_b_data.dat")
    do i = 1, num - 1
        write(2, "(f7.5,a3,f7.5)") rand_num(i), " - ", rand_num(i+1)
    end do
    close(2)

    !!! calculating correlation functions
    do k = 0, num -1
        sum_val = 0.0d0
        do i = 1, num - k
            sum_val = sum_val + rand_num(i)*rand_num(i+k)
        end do
        corr_list(k+1) = (sum_val/(num-k) - avg**2)/(sigma**2)
    end do

    open(unit = 3, file = "E:\computational_physics\Module_2_out\question_2_c_data.dat")
    do i = 1, num
        write(3, *) i-1, " - ", corr_list(i)
    end do
    close(3)

end program random_numbers


!!! code to write bins and frequency for histogram
subroutine graph_data(data_array, data_size, bin_size, bin_start, path, unit_num)
    implicit none

    character(len=*), intent(in) :: bin_start, path
    integer :: i, j
    integer :: num_bins, min_val, max_val
    integer, intent(in) :: data_size, unit_num
    real, allocatable :: bins(:), frequency(:)
    real, intent(in) :: bin_size
    real, intent(in) :: data_array(data_size)

    min_val=int(minval(data_array)-bin_size)

    if (bin_size==2.0 .and. (bin_start=="Even" .or. bin_start=="even") .and. mod(min_val,2) /= 0) then
        min_val=min_val-1
    else if (bin_size==2.0 .and. (bin_start=="Odd" .or. bin_start=="odd") .and. mod(min_val,2) == 0) then
        min_val=min_val-1
    end if

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
        do j=1,data_size
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
        write (unit_num, "(f7.5,a3,i10,a3,f7.5)") bins(i), " - ", nint(frequency(i)), " - ", frequency(i)/(data_size*bin_size)
    end do

    close(unit_num)

end subroutine graph_data