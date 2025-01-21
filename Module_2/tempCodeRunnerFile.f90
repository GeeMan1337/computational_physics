program cvtExample2
    implicit none
    character(30) :: str1
    integer :: iNum=246873
    real :: pi = 3.14, tau

    !  Convert integer value to a string.
    
        iNum = iNum / 100
        write (str1, '(i1)') iNum
        print*, iNum
    end program cvtExample2