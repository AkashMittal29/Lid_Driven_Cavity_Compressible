
PROGRAM main
    USE mod_parameters
    USE mod_solver
    USE omp_lib

    IMPLICIT NONE
    TYPE(solver_compressible) :: compressible

    !call omp_set_num_threads(3) 
    WRITE(*,'(A," ",I3)') 'Max threads = ', omp_get_max_threads()

    !$OMP PARALLEL
        WRITE(*,'(A," ",I3,"/",I3)') 'Thread id: ', omp_get_thread_num(), omp_get_num_threads()
        !$OMP BARRIER
        !$OMP MASTER
            WRITE(*,'(A," ",I3)') 'from master', omp_get_thread_num()
        !$OMP END MASTER
    !$OMP END PARALLEL

    ! Solution
    CALL compressible%read_input('input.ini')
    CALL compressible%construct()
    CALL compressible%initialize_grid()
    CALL compressible%initialize_field()
    CALL compressible%solve()

    
END PROGRAM main


