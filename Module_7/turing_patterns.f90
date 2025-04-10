program turing_patterns
    implicit none 

    integer :: i, j, k
    integer :: left, right, top, down, grid_size
    real*8 :: temp, temp_2, temp_3
    real*8, allocatable :: a_grid(:,:), b_grid(:,:)
    real*8, allocatable :: a_grid_old(:,:), b_grid_old(:,:)
    real*8 :: dt, h, diff_a, diff_b, alpha, beta   !h= dx =dy  !diff_a and diff_b are diffusion constants
    integer :: num_iter, num_snapshot

    num_iter = 10000; num_snapshot = 100
    grid_size = 100

    dt = 0.002d0; h = 1.0d0
    diff_a = 0.35d0; diff_b = 93.0d0
    alpha = 0.025d0; beta = 20.0d0

    allocate(a_grid(grid_size,grid_size))
    allocate(a_grid_old(grid_size,grid_size))
    allocate(b_grid(grid_size,grid_size))
    allocate(b_grid_old(grid_size,grid_size))

    do i = 1, grid_size
        do j = 1, grid_size
            call random_number(temp)
            call random_number(temp_2)
            a_grid(i,j) = temp
            b_grid(i,j) = temp_2
        end do 
    end do

    a_grid_old = a_grid; b_grid_old = b_grid

    do k = 1, num_iter
        
        if (mod(k,501) == 0) then       ! this drops heat onto the grid
            call random_number(temp); call random_number(temp_2); call random_number(temp_3)
            a_grid(nint(temp*grid_size),nint(temp_2*grid_size)) = temp_3
        end if

        do i = 1, grid_size
            do j = 1, grid_size
                
                left = i-1; right = i+1; top = j+1; down = j-1
                if (i == 1) left = grid_size
                if (i == grid_size) right = 1
                if (j == 1) down = grid_size
                if (j == grid_size) top = 1

                a_grid(i,j) = (1 - 4*dt*diff_a/h**2)*a_grid_old(i,j) &
                              + (diff_a*dt/h**2)*(a_grid_old(left,j) + a_grid_old(right,j) &
                                                  + a_grid_old(i,top) + a_grid_old(i,down)) &
                              + dt*(a_grid_old(i,j) - a_grid_old(i,j)**3 -b_grid_old(i,j) +alpha)

                b_grid(i,j) = (1 - 4*dt*diff_b/h**2)*b_grid_old(i,j) &
                              + (diff_b*dt/h**2)*(b_grid_old(left,j) + b_grid_old(right,j) &
                                                  + b_grid_old(i,top) + b_grid_old(i,down)) &
                              + dt*beta*(a_grid_old(i,j) - b_grid_old(i,j))

            end do

            a_grid_old = a_grid; b_grid_old = b_grid

        end do
    end do
    
    open(file = "E:\computational_physics\Module_7_out\test.csv", unit = 11)
    do i = 1, grid_size
        write (11,*) a_grid(1:grid_size,i)
    end do
    close(11)

    open(file = "E:\computational_physics\Module_7_out\test_2.csv", unit = 12)
    do i = 1, grid_size
        write (12,*) b_grid(1:grid_size,i)
    end do
    close(12)

end program turing_patterns