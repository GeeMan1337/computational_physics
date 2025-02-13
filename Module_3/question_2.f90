program use_mod
    use my_mod
    implicit none
  
    real :: mat(10, 10)
  
    mat(:,:) = public_var
  
    call print_matrix(mat)
  
  end program use_mod