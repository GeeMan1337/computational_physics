program diffusion_equation
    implicit none 

    integer :: i, j, k
    integer :: left, right, top, down, grid_size
    real*8, allocatable :: a_grid(:,:)
    real*8, allocatable :: a_grid_old(:,:)
    real*8 :: dt, h, diff_a         !h= dx =dy  !diff_a and diff_b are diffusion constants
    real*8 :: temp, temp_2, temp_3
    integer :: num_iter, num_snapshot

    integer :: unit_a, unit_b, num_file
    character(len=100) :: filename_a, filename_b

    num_iter = 20000; num_snapshot = 100; grid_size = 50

    dt = 0.002d0; h = 1.0d0
    diff_a = 0.75d0

    allocate(a_grid(grid_size,grid_size))
    allocate(a_grid_old(grid_size,grid_size))

    a_grid = 0.0d0     
    
    ! setting up initial conditions  
    call random_number(temp); call random_number(temp_2); call random_number(temp_3)
    a_grid(nint(temp*grid_size),nint(temp_2*grid_size)) = 5 + temp_3

    a_grid_old = a_grid
    
    ! writing initial conditions
    num_file = 1; unit_a = 10; unit_b = 11
    write(filename_a, '("E:\computational_physics\Module_7_out\diffusion_data\diffusion_a_", i0 ,".dat")') num_file

    open(unit=unit_a, file=filename_a)
    do i = 1, grid_size
        do j = 1, grid_size
            write(unit_a,*) i, j, a_grid_old(i,j)
        end do
    end do
    close(unit_a)

    unit_a = unit_a + 2; unit_b = unit_b + 2; num_file = num_file + 1

    do k = 1, num_iter
        
        if (mod(k,501) == 0) then       ! this drops heat onto the grid
            call random_number(temp); call random_number(temp_2); call random_number(temp_3)
            a_grid_old(nint(temp*grid_size),nint(temp_2*grid_size)) = 5 + temp_3
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
                                                  + a_grid_old(i,top) + a_grid_old(i,down))
            end do
        end do

        a_grid_old = a_grid

        if (mod(k, num_snapshot) == 0) then

            write(filename_a, '("E:\computational_physics\Module_7_out\diffusion_data\diffusion_a_", i0 ,".dat")') num_file
            open(unit=unit_a, file=filename_a)

            do i = 1, grid_size
                do j = 1, grid_size
                    write(unit_a,*) i, j, a_grid_old(i,j)
                end do
            end do

            close(unit_a)
            unit_a = unit_a + 2; unit_b = unit_b + 2; num_file = num_file + 1
        
        end if
    end do
    
end program diffusion_equation