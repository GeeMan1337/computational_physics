program question_1
    use ising
    implicit none

    !real :: j_ising, kt, start_config
    !integer :: len, niter, unit_num
    !character(len=*) :: path

    call ising_model(2,1.0,1.0,10**5,"abc",1,"E:\computational_physics\Module_3_out\question_1_data.dat")
    call ising_model(2,1.0,0.001,10**5,"abc",2,"E:\computational_physics\Module_3_out\question_1_data.dat")


end program question_1