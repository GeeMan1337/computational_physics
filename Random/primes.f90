program primes 
    use iso_fortran_env
    implicit none

    integer :: num, i
    integer :: upper_lim
    integer :: count
    integer, allocatable :: prime_list(:)
    integer, external :: is_prime

    count = 0
    upper_lim = 10**5
    allocate(prime_list(upper_lim/2))

    prime_list(1) = 2
    count = count + 1

    do num = 3, upper_lim, 2
        if (is_prime(num) == 1) then
            prime_list(count + 1) = num
            count = count + 1
        end if
    end do

    open(unit = 1, file = "E:\computational_physics\Random\prime_list.dat")

    do i = 1, count 
        write (1,*) prime_list(i)
    end do

    close(1)
    
end program primes

integer function is_prime(number) result(flag)
    implicit none 

    integer :: i
    integer, intent(in) :: number 

    do i = 2, int(sqrt(real(number)))
        if (mod(number, i) == 0) then
            flag = 0
            go to 1
        end if
    end do

    flag = 1

1 end function is_prime 
