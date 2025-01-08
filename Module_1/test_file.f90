!!!!!!!!!!!!! Question 1.h  !!!!!!!!!!!!!!!!!
program name
    implicit none
    integer :: i, j, no_bins=50
    real*8 :: random(100000), sum_(100000), s, m, range, binwidth 
    
    real*8, allocatable :: bincenter(:),binedges(:), hist(:), norm_hist(:)
    allocate(binedges(no_bins+1))
    allocate(bincenter(no_bins))
    allocate(hist(no_bins))
    allocate (norm_hist(no_bins))
    sum_=0
    do i=1,100000
        call random_number(random)
        end do
        sum_=sum(random)
   
    open(unit=5, file= 'distv.dat')

    s=maxval(sum_)
    print*, s 
    m= minval(sum_)
    print*, m  
    
    range=s-m 
    binwidth=range/no_bins

    do j=1,no_bins
        binedges(j)=m+binwidth*(j-1)
    end do
    binedges(no_bins+1)=s  ! assign max value to last bin point

    
    hist=0

    do i=1,100000
        do j=1,no_bins
        if ( sum_(i)>=binedges(j) .and. sum_(i)<binedges(j+1) ) then
            hist(j)=hist(j)+1
        end if
    end do 
end do
s=sum(hist)
do j=1,no_bins
    norm_hist(j)=hist(j)/s
end do

write(5,*) 'bin center', 'normalized freq'
do i=1,no_bins
    bincenter(i)=binedges(i)+binwidth/2
    write(5,*) bincenter(i), norm_hist(i)
end do




end program name