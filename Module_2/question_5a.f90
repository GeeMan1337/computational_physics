program multi_dim_integral
    implicit none
    
    integer :: i, j, k, num(7)
    real*8 :: integral, length, x(6), volume, temp, sum_val, G, avg
    real*8 :: xx, yy, xy
    real*8 :: sigma, sum_temp

    num = [10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7]
    length = 5.0d0
    volume = (2*length)**6

    print *, "Brute force method with integral limits from -5 to 5"
    print *, " "

    do k = 1, 7
        sum_val = 0.0d0
        sum_temp = 0.0d0
        avg = 0.0d0
        do i = 1, num(k)
            do j = 1, 6
                call random_number(temp)
                x(j) = 2*length*temp - length
            end do
            xx= x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
            yy= x(4)*x(4) + x(5)*x(5) + x(6)*x(6)
            xy= (x(1)-x(4))**2 + (x(2)-x(5))**2 + (x(3)-x(6))**2
            G = exp(- xx - yy - 0.5d0*xy)
            sum_temp = sum_temp + G**2
            sum_val = sum_val + G
        end do
    avg = sum_val/num(k)
    sigma = sum_temp/num(k) - avg**2
    integral = volume * avg
    sigma = volume * sigma/sqrt(1.0d0*num(k))
    print *, "Integral with", num(k), "samples is", integral
    print *, "Error is", sigma
    print *, " "
    
    end do

end program multi_dim_integral