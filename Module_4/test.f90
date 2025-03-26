program test
    implicit none

    integer :: i, list(3,5), s(5), left
    real*8 :: a, b

    a = 5.342424113242424
    b = 1.344267565352541
    s = [1,2,3,4,5]

    print*, abs(b-a)

end program test