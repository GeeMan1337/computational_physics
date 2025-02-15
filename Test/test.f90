program test
    implicit none

    real :: num =7.0, temp_2, temp_3, temp_4 
    integer :: i = 5, len=5

    call random_seed()
    call random_number(temp_2)
    call random_number(temp_3)
    call random_number(temp_4)
    temp_2 = temp_2*len + 1; temp_3 = temp_3*len + 1; temp_4 = temp_4*len + 1


    print*, temp_2, int(temp_2)


end program test