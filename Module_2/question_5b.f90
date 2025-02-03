program multi_dim_integral
    implicit none
    
    integer :: i, j, k, num(7)
    real*8 :: integral, x(6), temp_1, temp_2, sum_val, F, avg
    real*8 :: sigma, sum_temp
    real*8, parameter :: pi = 2.0d0 * asin(1.0d0)

    num = [10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7]

    print *, "Importance sampling method"
    print *, " "

    do k = 1, 7
        sum_val = 0.0d0
        sum_temp = 0.0d0
        avg = 0.0d0
        do i = 1, num(k)
            do j = 1, 6
                call random_number(temp_1)
                call random_number(temp_2)
                x(j) = (1/sqrt(2.0d0))*sqrt(-2*log(temp_1))*cos(2*pi*temp_2)
            end do
            F = exp(-0.5d0*((x(1)-x(4))**2 + (x(2)-x(5))**2 + (x(3)-x(6))**2))
            sum_val = sum_val + F
            sum_temp = sum_temp + F**2
        end do
    avg = sum_val/num(k)
    sigma = sum_temp/num(k) - avg**2
    integral = pi**3 * avg
    sigma= pi**3 * sigma/sqrt(1.0d0*num(k))
    
    print *, "Integral with", num(k), "samples is", integral
    print *, "Error is", sigma
    print *, " "
    
    end do

end program multi_dim_integral