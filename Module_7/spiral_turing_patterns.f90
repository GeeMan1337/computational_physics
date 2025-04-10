program spiral_patterns
    implicit none 

    integer :: i, j, k
    integer :: left, right, top, down, grid_size, num_spirals
    real*8 :: temp, temp_2
    real*8, allocatable :: a_grid(:,:), b_grid(:,:)
    real*8, allocatable :: a_grid_old(:,:), b_grid_old(:,:)
    real*8 :: dt, h, diff_a, diff_b, alpha, beta   !h= dx =dy  !diff_a and diff_b are diffusion constants
    integer :: num_iter, num_snapshot

    num_iter = 20000; num_snapshot = 100
    grid_size = 100
    num_spirals = 4

    dt = 0.001d0; h = 1.0d0
    diff_a = 0.5d0; diff_b = 100.0d0
    alpha = 0.005d0; beta = 5.0d0

    allocate(a_grid(grid_size,grid_size))
    allocate(a_grid_old(grid_size,grid_size))
    allocate(b_grid(grid_size,grid_size))
    allocate(b_grid_old(grid_size,grid_size))

    a_grid = 0.0d0
    b_grid = 0.0d0

    do i = 1, num_spirals
        call random_number(temp); call random_number(temp_2)
        a_grid(nint(grid_size*temp),nint(grid_size*temp_2)) = 0.4d0
    end do

    a_grid_old = a_grid; b_grid_old = b_grid

    do k = 1, num_iter
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

end program spiral_patterns