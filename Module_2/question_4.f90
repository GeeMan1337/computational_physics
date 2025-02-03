program prob_dist_func
    implicit none

    integer :: i, k, num 
    real, allocatable :: rand_num(:), rand_num_square(:), corr_list(:)
    real :: sigma, avg, temp, a_1, a_2
    real, parameter :: pi = 2.0 * asin(1.0)

    num = 10**5

    allocate(rand_num(num))
    allocate(rand_num_square(num))
    allocate(corr_list(num))

    !!! this part generates a the PDF 2 * exp(-2*x)
    do i = 1, num
        call random_seed()
        call random_number(temp)
        !!! this formula generates a PDF of 2 * exp(-2*x)
        rand_num(i) = - 0.5*log(1 - temp)
    end do

    !!! list with square of the random numbers with PDF of 2 * exp(-2*x)
    do i = 1, num
        rand_num_square(i) = rand_num(i)**2
    end do

    !!! this is standard deviation and average of PDF 2 * exp(-2*x)
    avg = (sum(rand_num)/num)
    sigma = sqrt(sum(rand_num_square)/num - avg**2)

    print *, "The average of", num, "random numbers is (distribution is 2*exp(-2x))", avg 
    print *, "The standard deviation of", num, "random numbers is (distribution is 2*exp(-2x))", sigma

    call graph_data(rand_num,num,0.01,"even","E:\computational_physics\Module_2_out\question_4a_data.dat",1)


    !!! this part generates a Gaussian PDF with standard deviation 2
    do i = 1, num
        call random_seed()
        call random_number(a_1)
        call random_number(a_2)
        !!! this formula generates a PDF of 1/sqrt(8*pi) exp(-0.5*(x/2)^2)
        rand_num(i) = 2.0*sqrt(-2*log(a_1))*cos(2*pi*a_2)
    end do

    !!! list with square of the random numbers with PDF of 1/sqrt(8*pi) exp(-0.5*(x/2)^2)
    do i = 1, num
        rand_num_square(i) = rand_num(i)**2
    end do

    !!! this is standard deviation and average of PDF 1/sqrt(8*pi) exp(-0.5*(x/2)^2)
    avg = (sum(rand_num)/num)
    sigma = sqrt(sum(rand_num_square)/num - avg**2)

    print *, "The average of", num, "random numbers is (distribution is Gaussian with SD=2)", avg 
    print *, "The standard deviation of", num, "random numbers is (distribution is Gaussian with SD=2)", sigma

    call graph_data(rand_num,num,0.1,"even","E:\computational_physics\Module_2_out\question_4b_data.dat",2)

end program prob_dist_func


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
        write (unit_num, *) bins(i), " - ", nint(frequency(i)), " - ", frequency(i)/(data_size*bin_size)
    end do

    close(unit_num)

end subroutine graph_data