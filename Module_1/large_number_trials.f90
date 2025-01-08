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

    open(unit=1, file="E:\computational_physics\Module_1_out\question_1h.csv")
    
    do i=1,10000
        write (1,*) sum_list(i)
    end do
    close(1)

    !!! 10^4 trials with 10^4 random numbers between -1 and 1 
    do i=1,10000
        call random_seed()
        do j=1,10000
            call random_number(rand_num(j))
            rand_num(j)=2.0*rand_num(j)     !doubling the random number between 0 and 1
            if (rand_num(j)>1) then
                rand_num(j)=rand_num(j)-2       !if the doubled random number is greater than 1 then subrtacting 2 to put it in -1 to 0 range
            end if
        end do
        sum_list(i)=sum(rand_num)
    end do

    open(unit=2, file="E:\computational_physics\Module_1_out\question_1h_part_b.csv")
    
    do i=1,10000
        write (2,*) sum_list(i)
    end do
    close(2)

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

    open(unit=3, file="E:\computational_physics\Module_1_out\question_1h_part_c.csv")
    
    do i=1,100000
        write (3,*) sum_bigger_list(i)
    end do
    close(3)

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

    open(unit=4, file="E:\computational_physics\Module_1_out\question_1i.csv")
    
    do i=1,10000
        write (4,*) sum_list(i)
    end do
    close(4)

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

    open(unit=5, file="E:\computational_physics\Module_1_out\question_1k.csv")
    
    do i=1,100000
        write (5,*) sum_bigger_list(i)
    end do
    close(5)

end program large_number_trials