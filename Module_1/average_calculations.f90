program average_calculations
    implicit none

    integer :: i        !iterating variable
    real *8 :: avg_10, avg_100, avg_10000, avg_1000000
    real *8 :: rand_10(10), rand_100(100), rand_10000(10000), rand_1000000(1000000)
    real *8 :: half=0.5, diff_10, diff_100, diff_10000, diff_1000000

    !!! average of 10 random numbers
    call random_seed()
    do i=1,10
        call random_number(rand_10(i))
    end do

    avg_10=real(sum(rand_10))/size(rand_10)

    !!! average of 100 random numbers
    call random_seed()
    do i=1,100
        call random_number(rand_100(i))
    end do

    avg_100=real(sum(rand_100))/size(rand_100)

    !!! average of 10000 random numbers
    call random_seed()
    do i=1,10000 
        call random_number(rand_10000(i))
    end do

    avg_10000=real(sum(rand_10000))/size(rand_10000)

    !!! average of 1000000 random numbers
    call random_seed()
    do i=1,1000000
        call random_number(rand_1000000(i))
    end do

    avg_1000000=real(sum(rand_1000000))/size(rand_1000000)

    open(unit=1, file="E:\computational_physics\Module_1_out\test_ran.dat", access="append")
    write (1,"(/,a,/)") "Now calculating average of 10, 100, 10000, 1000000 random numbers:"
    write (1,*) "Average of 10 random numbers = ", avg_10
    write (1,*) "Average of 100 random numbers = ", avg_100
    write (1,*) "Average of 10000 random numbers = ", avg_10000
    write (1,*) "Average of 1000000 random numbers = ", avg_1000000

    write (1,"(/,a,/)") "Absolute value of difference between 0.5 and mean of:"

    write (1,*) "10 random numbers = ", abs(half-avg_100)
    write (1,*) "100 random numbers = ", abs(half-avg_100)
    write (1,*) "10000 random numbers = ", abs(half-avg_10000)
    write (1,*) "1000000 random numbers = ", abs(half-avg_1000000)

    write (1,"(/,a)") "It is clear that as we run more trials &
    on the uniform probability distribution, the mean value of the outcomes tend towards the mid point of the distribution."

end program average_calculations