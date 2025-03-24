program test
    implicit none

    integer :: i, list(3,5), s(5), left

    s = [1,2,3,4,5]

    list = 0
    list(1,1:5) = s
   ! list(2,1:5) = s*10
   ! list(3,1:5) = s*100

    do  i = 1,3
        write(1,*) list(i,1), list(i,2), list(i,3), list(i,4), list(i,5)
    end do

    do i = 1, 7
        if (i==3) then 
            left = 69

        else 
            left = 30
        end if
        print *, i, left
    end do
    


end program test