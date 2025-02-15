program sieve_primes 
    use iso_fortran_env
    implicit none

    integer :: num, i
    integer :: upper_lim
    integer, allocatable :: prime_list(:)

    upper_lim = 10**3
    
    allocate(prime_list(upper_lim))
    prime_list(1) = 0

    do i = 2, upper_lim
        prime_list(i) = i
    end do

    do num = 2, upper_lim
        if (prime_list(num) /= 0) then
            do i = 2, upper_lim/num
                prime_list(num*i) = 0
            end do
        end if
    end do

    open(unit = 1, file = "E:\computational_physics\Random\prime_list.dat")

    do i = 1, upper_lim
        if (prime_list(i) /= 0) then
            write(1,*) prime_list(i)
        end if
    end do

    close(1)

end program sieve_primes