!!! graph_data subroutine in the end 

program large_number_trials
    implicit none

    integer :: i, j        !iterating variable
    real *8 :: sum_list(10000), rand_num(10000)
    real *8 :: sum_bigger_list(100000), rand_bigger_num(100000)

    !!! question 1h
    !!! finding the sum of 10^4 trials 10^4 times
    do i=1,10000
        call random_seed()
        do j=1,10000
            call random_number(rand_num(j))
        end do
        sum_list(i)=sum(rand_num)
    end do

    !open(unit=1, file="E:\computational_physics\Module_1_out\question_1h.csv")
    
    !do i=1,10000
    !    write (1,*) sum_list(i)
    !end do
    !close(1)

    call graph_data(sum_list,10000,0.5,"even","E:\computational_physics\Module_1_out\graph_data_1h_1.dat",11)
    call graph_data(sum_list,10000,1.0,"even","E:\computational_physics\Module_1_out\graph_data_1h_2.dat",12)
    call graph_data(sum_list,10000,2.0,"even","E:\computational_physics\Module_1_out\graph_data_1h_3.dat",13)

    !!! 10^4 trials with 10^4 random numbers between -1 and 1 
    do i=1,10000
        call random_seed()
        do j=1,10000
            call random_number(rand_num(j))
            rand_num(j)=2.0*rand_num(j)-1.0     !doubling the random number in [0,1] and subtracting 1 to get it in range [-1,1]
        end do
        sum_list(i)=sum(rand_num)
    end do

    !open(unit=2, file="E:\computational_physics\Module_1_out\question_1h_part_b.csv")
    
    !do i=1,10000
    !    write (2,*) sum_list(i)
    !end do
    !close(2)

    call graph_data(sum_list,10000,2.0,"even","E:\computational_physics\Module_1_out\graph_data_1h_4.dat",14)

    !!! 10^5 trials with 10^4 random numbers between -1 and 1 (10 times the previous one)
    do i=1,100000
        call random_seed()
        do j=1,10000
            call random_number(rand_num(j))
            rand_num(j)=2.0*rand_num(j)     !doubling the random number between 0 and 1
            if (rand_num(j)>1) then
                rand_num(j)=rand_num(j)-2       !if the doubled random number is greater than 1 then subrtacting 2 to put it in -1 to 0 range
            end if
        end do
        sum_bigger_list(i)=sum(rand_num)
    end do

    !open(unit=3, file="E:\computational_physics\Module_1_out\question_1h_part_c.csv")
    
    !do i=1,100000
    !    write (3,*) sum_bigger_list(i)
    !end do
    !close(3)

    call graph_data(sum_bigger_list,100000,2.0,"even","E:\computational_physics\Module_1_out\graph_data_1h_5.dat",15)

    !!! question 1i
    !!! 10^4 trials with 10^4 random numbers either -1 and 1
    do i=1,10000
        call random_seed()
        do j=1,10000
            call random_number(rand_num(j))
            if (rand_num(j)>0.5) then
                rand_num(j)=-1
            else
                rand_num(j)=1
            end if
        end do
        sum_list(i)=sum(rand_num)
    end do

    !open(unit=4, file="E:\computational_physics\Module_1_out\question_1i.csv")
    
    !do i=1,10000
    !    write (4,*) sum_list(i)
    !end do
    !close(4)

    call graph_data(sum_list,10000,1.0,"even","E:\computational_physics\Module_1_out\graph_data_1i_1.dat",16)
    call graph_data(sum_list,10000,2.0,"even","E:\computational_physics\Module_1_out\graph_data_1i_2.dat",17)
    call graph_data(sum_list,10000,2.0,"odd","E:\computational_physics\Module_1_out\graph_data_1i_3.dat",18)
    call graph_data(sum_list,10000,5.0,"even","E:\computational_physics\Module_1_out\graph_data_1i_4.dat",19)
    call graph_data(sum_list,10000,10.0,"even","E:\computational_physics\Module_1_out\graph_data_1i_5.dat",20)

    !!! question 1k
    !!! 10^5 trials with 10^5 random numbers either -1 and 1
    do i=1,100000
        call random_seed()
        do j=1,100000
            call random_number(rand_bigger_num(j))
            if (rand_bigger_num(j)>0.5) then
                rand_bigger_num(j)=-1
            else
                rand_bigger_num(j)=1
            end if
        end do
        sum_bigger_list(i)=sum(rand_bigger_num)
    end do

    !open(unit=5, file="E:\computational_physics\Module_1_out\question_1k.csv")
    
    !do i=1,100000
    !    write (5,*) sum_bigger_list(i)
    !end do
    !close(5)

    call graph_data(sum_bigger_list,100000,2.0,"even","E:\computational_physics\Module_1_out\graph_data_1k_1.dat",21)

end program large_number_trials

!!! code to write bins and frequency for histogram
subroutine graph_data(data_array, data_size, bin_size, bin_start, path, unit_num)
    implicit none

    character(len=*), intent(in) :: bin_start, path
    integer :: i, j
    integer :: num_bins, min_val, max_val
    integer, intent(in) :: data_size, unit_num
    real, allocatable :: bins(:), frequency(:)
    real, intent(in) :: bin_size
    real *8, intent(in) :: data_array(data_size)

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