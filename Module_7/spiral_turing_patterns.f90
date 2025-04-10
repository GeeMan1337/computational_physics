program spiral_patterns !uses Fitzhugh Nagumo equations
    implicit none 

    integer :: i, j, k
    integer :: left, right, top, down, grid_size, num_spirals
    real*8 :: temp, temp_2
    real*8, allocatable :: a_grid(:,:), b_grid(:,:)
    real*8, allocatable :: a_grid_old(:,:), b_grid_old(:,:)
    real*8 :: dt, h, diff_a, diff_b, alpha, beta   !h= dx =dy  !diff_a and diff_b are diffusion constants
    integer :: num_iter, num_snapshot

    integer :: unit_a, unit_b, num_file
    character(len=100) :: filename_a, filename_b

    num_iter = 40000; num_snapshot = 400
    grid_size = 100
    num_spirals = 10

    dt = 0.002d0; h = 1.0d0
    diff_a = 0.30d0; diff_b = 70.0d0
    alpha = 0.005d0; beta = 3.0d0

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

    ! writing initial conditions
    num_file = 0; unit_a = 10; unit_b = 11
    write(filename_a, '("E:\computational_physics\Module_7_out\spiral_data\spiral_a_", i0 ,".dat")') num_file
    write(filename_b, '("E:\computational_physics\Module_7_out\spiral_data\spiral_b_", i0 ,".dat")') num_file

    open(unit=unit_a, file=filename_a)
    open(newunit=unit_b, file=filename_b)
    do i = 1, grid_size
        do j = 1, grid_size
            write(unit_a,*) i, j, a_grid_old(i,j)
            write(unit_b,*) i, j, b_grid_old(i,j)
        end do
    end do
    close(unit_a)
    close(unit_b)

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
        end do
        a_grid_old = a_grid; b_grid_old = b_grid

        num_file = k

        if (mod(k, num_snapshot) == 0) then

            write(filename_a, '("E:\computational_physics\Module_7_out\spiral_data\spiral_a_", i0 ,".dat")') num_file
            write(filename_b, '("E:\computational_physics\Module_7_out\spiral_data\spiral_b_", i0 ,".dat")') num_file

            open(unit=unit_a, file=filename_a)
            open(unit=unit_b, file=filename_b)

            do i = 1, grid_size
                do j = 1, grid_size
                    write(unit_a,*) i, j, a_grid_old(i,j)
                    write(unit_b,*) i, j, b_grid_old(i,j)
                end do
            end do

            close(unit_a)
            close(unit_b)
        
        end if
    end do

end program spiral_patterns