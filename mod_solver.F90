
MODULE mod_solver
    USE mod_read_file
    USE mod_grid
    USE mod_field
    USE mod_parameters
    USE mod_utility
    USE omp_lib
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : error_unit 
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: solver_compressible
    
    TYPE :: solver_compressible
        TYPE(input_data_type) :: input_data
        TYPE(grid_type),  POINTER, PUBLIC :: grid => NULL()
        TYPE(field_type), POINTER, PUBLIC :: field => NULL()
        CHARACTER(LEN=:), ALLOCATABLE :: grid_file
        REAL(KIND=rkind)  :: Re_ref                    ! Reference Reynolds number
        REAL(KIND=rkind)  :: M_ref                     ! Reference Mach number
        REAL(KIND=rkind)  :: Pr                        ! Prandtl number
        REAL(KIND=rkind)  :: gamma                     ! Gamma for air (Cp/Cv)
        REAL(KIND=rkind)  :: tolerance                 ! Tolerance for convergence check
        INTEGER           :: niter                     ! Number of time iterations
        REAL(KIND=rkind)  :: eps4                      ! Factor for 4th order numerical diffusion (for interior nodes)
        REAL(KIND=rkind)  :: eps2                      ! Factor for 2nd order numerical diffusion (for boundary nodes)
        REAL(KIND=rkind)  :: dt_factor                 ! Factor to assign dt using dt=dt_factor*min(dx,dy) for numerical stability
        REAL(KIND=rkind)  :: dt                        ! Non-dim time-step (is defined based on grid as dt=dt_factor*min(dx,dy)


        CONTAINS
            PROCEDURE, PASS(self), PUBLIC  :: read_input
            PROCEDURE, PASS(self), PUBLIC  :: construct
            PROCEDURE, PASS(self), PUBLIC  :: initialize_grid
            PROCEDURE, PASS(self), PUBLIC  :: initialize_field
            PROCEDURE, PASS(self), PUBLIC  :: solve
            PROCEDURE, PASS(self), PUBLIC  :: compute_p_mu
            PROCEDURE, PASS(self), PUBLIC  :: conserved_from_primitive_full
            PROCEDURE, PASS(self), PUBLIC  :: conserved_from_primitive_boundary
            PROCEDURE, PASS(self), PUBLIC  :: primitive_from_conserved
            PROCEDURE, PASS(self), PUBLIC  :: solve_continuity
            PROCEDURE, PASS(self), PUBLIC  :: solve_x_momentum
            PROCEDURE, PASS(self), PUBLIC  :: solve_y_momentum
            PROCEDURE, PASS(self), PUBLIC  :: solve_energy
            PROCEDURE, PASS(self), PUBLIC  :: add_numerical_diffusion
            PROCEDURE, PASS(self), PUBLIC  :: conserved_weighted_sum
            GENERIC :: conserved_from_primitive => conserved_from_primitive_full, conserved_from_primitive_boundary
    END TYPE solver_compressible


    CONTAINS
        ! Procedures of user defined type: solver_compressible
        SUBROUTINE read_input(self, file)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: file

            ! Reading input data and storing into input_data object
            CALL self%input_data%read_input(file)
        END SUBROUTINE read_input

        
        SUBROUTINE construct(self)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            
            IF(.NOT. ASSOCIATED(self%grid))  ALLOCATE(self%grid)
            IF(.NOT. ASSOCIATED(self%field)) ALLOCATE(self%field)

            ! Assigning variables from the input file now stored in input_data object
            CALL self%input_data%assign_variable('Grid information','grid_file',self%grid_file)
            CALL self%input_data%assign_variable('Flow parameters','Re_ref',self%Re_ref)
            CALL self%input_data%assign_variable('Flow parameters','M_ref',self%M_ref)
            CALL self%input_data%assign_variable('Fluid properties','Pr',self%Pr)
            CALL self%input_data%assign_variable('Fluid properties','gamma',self%gamma)
            CALL self%input_data%assign_variable('Computational parameters','dt_factor',self%dt_factor)
            CALL self%input_data%assign_variable('Computational parameters','niter',self%niter)
            CALL self%input_data%assign_variable('Computational parameters','tolerance',self%tolerance)
            CALL self%input_data%assign_variable('Computational parameters','eps4',self%eps4)
            CALL self%input_data%assign_variable('Computational parameters','eps2',self%eps2)

            PRINT*, 'Input data:'
            CALL self%input_data%print_input_data()
            PRINT*
        END SUBROUTINE construct


        SUBROUTINE initialize_grid(self)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            
            IF(self%grid_file == '') THEN
                WRITE(error_unit,'("Error: Grid file is not defined in the input file. Aborted.")')
                STOP
            END IF
            CALL self%grid%initialize(self%grid_file)
        END SUBROUTINE initialize_grid


        SUBROUTINE initialize_field(self)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            
            CALL self%field%construct(self%grid) ! BC is also set
            
            ! Assigning p, mu, and Pr
            CALL self%compute_p_mu()
            self%field%Pr = self%Pr ! Pr is assumed constant

            ! Computing conserved variables, q(:,:,1)->rho*u, q(:,:,2)->rho*v, and q(:,:,3)->rho*E
            CALL self%conserved_from_primitive()
        END SUBROUTINE initialize_field


        SUBROUTINE solve(self)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            TYPE(field_type), POINTER :: field_old => NULL(), field_swap => NULL()
            TYPE(field_type), POINTER :: field_inter => NULL() ! For RK3
            INTEGER :: iter, i_rk3

            ! Allocating fields by copying initial field
            CALL copy_field(field_old,self%field)
            ! CALL copy_field(field_swap,self%field) ! Used for swapping between current and old field pointers
            CALL copy_field(field_inter,self%field)

            ! Assigning dt based on minimum grid spacing
            self%dt = self%dt_factor*( MIN( minval(self%grid%dx),minval(self%grid%dy) ) )
            PRINT*, 'dt = ', self%dt

            PRINT*; CALL write_date_time(); PRINT*;  
            WRITE(*,'("iter rho_rms u_rms_diff v_rms_diff p_rms_diff")')   

            iter = 1
            loop_time_iter: DO WHILE(iter<=self%niter)
                WRITE(*,'(I8," ")',ADVANCE='NO') iter
                field_swap => field_inter
                field_inter => field_old

                !#TODO: Applying RK3 
                rk3_loop: DO i_rk3=1,3
                    ! Solving governing equations to obtain new values of the conserved variables
                    !$OMP PARALLEL
                    !$OMP SECTIONS
                        !$OMP SECTION
                        CALL self%solve_continuity(field_inter)

                        !$OMP SECTION
                        CALL self%solve_x_momentum(field_inter)

                        !$OMP SECTION
                        CALL self%solve_y_momentum(field_inter)

                        !$OMP SECTION
                        CALL self%solve_energy(field_inter)
                    !$OMP END SECTIONS
                    !$OMP END PARALLEL

                    CALL self%add_numerical_diffusion(field_inter) ! For energy damping for (rho*E) is with derivative of (rho*E+p)

                    CALL self%primitive_from_conserved()
                    CALL self%field%set_bc()
                    CALL self%compute_p_mu()
                    ! Computing conserved variables, q(:,:,1)->rho*u, q(:,:,2)->rho*v, and q(:,:,3)->rho*E 
                    ! only at boundary and ghost nodes
                    CALL self%conserved_from_primitive('boundary_and_ghosts')

                    ! Weighted sum over conserved variables only.
                    IF(i_rk3==1) THEN
                        ! pointing field_inter to the currently updated field
                        field_inter => self%field
                        self%field => field_swap
                    ELSE
                        ! pointing field_inter to the currently updated field
                        field_swap => field_inter; field_inter => self%field; self%field => field_swap
                        IF(i_rk3==2) THEN ! self%field = 3/4*field_old + 1/4*field_inter
                            CALL self%conserved_weighted_sum(3.0/4.0, field_old, 1.0/4.0, field_inter)
                        ELSE              ! self%field = 1/3*field_old + 2/3*field_inter
                            CALL self%conserved_weighted_sum(1.0/3.0, field_old, 2.0/3.0, field_inter)
                        END IF
                        ! pointing field_inter to the currently updated field
                        field_swap => field_inter; field_inter => self%field; self%field => field_swap
                    END IF
                END DO rk3_loop
                self%field => field_inter
                field_inter => field_swap

                ! Print rms difference
                CALL print_rms_diff_fields(self%field, field_old, FORMAT='E13.6')

                ! Update field
                field_swap => field_old
                field_old  => self%field
                self%field => field_swap

                iter = iter+1
                PRINT*
            END DO loop_time_iter

            PRINT*, 'Solution complete.'
            PRINT*; CALL write_date_time(); PRINT*;
            PRINT*,'New field'
            CALL self%field%printf(FORMAT='E13.6')

            DEALLOCATE(field_old, field_inter)
        END SUBROUTINE solve


        SUBROUTINE compute_p_mu(self)
            CLASS(solver_compressible), INTENT(INOUT) :: self
        
            self%field%p  = self%field%rho*self%field%T/( self%gamma * (self%M_ref)**2 )
            self%field%mu = (self%field%T)**0.7
        END SUBROUTINE compute_p_mu


        SUBROUTINE conserved_from_primitive_full(self)
            ! Calculates conserved variables from primitive for interior nodes
            CLASS(solver_compressible), INTENT(INOUT) :: self

            ASSOCIATE(q => self%field%q, rho => self%field%rho, &
                    & u => self%field%u, v => self%field%v, T => self%field%T, &
                    & M_ref => self%M_ref, gamma => self%gamma)
                q(:,:,1) = rho*u
                q(:,:,2) = rho*v
                q(:,:,3) = rho*( T/( (M_ref**2)*gamma*(gamma-1) ) + 0.5*(u**2+v**2) ) ! non-dim rho*Total_energy(E) 
            END ASSOCIATE
        END SUBROUTINE conserved_from_primitive_full


        SUBROUTINE conserved_from_primitive_boundary(self,str)
            ! Calculates conserved variables from primitive on domain boundaries and ghost nodes.
            CLASS(solver_compressible), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: str ! Any string is accepted. Only to differ from conserved_from_primitive_full
            INTEGER :: ind(4)

            ASSOCIATE(q => self%field%q, rho => self%field%rho, &
                    & u => self%field%u, v => self%field%v, T => self%field%T, &
                    & M_ref => self%M_ref, gamma => self%gamma, nx => self%grid%nx, ny => self%grid%ny)
                ! Bottom and top bounds
                ind = (/0,1,ny,ny+1/)
                q(0:nx+1,ind,1) = rho(0:nx+1,ind)*u(0:nx+1,ind)
                q(0:nx+1,ind,2) = rho(0:nx+1,ind)*v(0:nx+1,ind)
                q(0:nx+1,ind,3) = rho(0:nx+1,ind)                                 &
                                & *( T(0:nx+1,ind)/( (M_ref**2)*gamma*(gamma-1) ) &
                                &   +0.5*( u(0:nx+1,ind)**2+v(0:nx+1,ind)**2 ) ) ! non-dim rho*Total_energy(E) 

                ! Left and right boundaries
                ind = (/0,1,nx,nx+1/)
                q(ind,0:ny+1,1) = rho(ind,0:ny+1)*u(ind,0:ny+1)
                q(ind,0:ny+1,2) = rho(ind,0:ny+1)*v(ind,0:ny+1)
                q(ind,0:ny+1,3) = rho(ind,0:ny+1)                                 &
                                & *( T(ind,0:ny+1)/( (M_ref**2)*gamma*(gamma-1) ) &
                                &   +0.5*( u(ind,0:ny+1)**2+v(ind,0:ny+1)**2 ) ) ! non-dim rho*Total_energy(E) 
            END ASSOCIATE
        END SUBROUTINE conserved_from_primitive_boundary


        SUBROUTINE primitive_from_conserved(self)
            CLASS(solver_compressible), INTENT(INOUT) :: self

            ASSOCIATE(q => self%field%q, rho => self%field%rho, &
                    & u => self%field%u, v => self%field%v, T => self%field%T, &
                    & M_ref => self%M_ref, gamma => self%gamma)
                u = q(:,:,1)/rho
                v = q(:,:,2)/rho
                T = ( q(:,:,3)/rho-0.5*(u**2+v**2) )*(M_ref**2)*gamma*(gamma-1) 
            END ASSOCIATE
        END SUBROUTINE primitive_from_conserved


        SUBROUTINE solve_continuity(self, field_old)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            TYPE(field_type), POINTER :: field_old
            INTEGER :: i, j
            
            ! PRINT*,'in continuity ',omp_get_thread_num()
            ASSOCIATE(rho => field_old%rho, u => field_old%u, v => field_old%v, &
                & T => field_old%T, p => field_old%p, mu => field_old%mu, Pr => field_old%Pr, &
                & q => field_old%q, dt => self%dt, &
                & nx => field_old%grid%nx, ny => field_old%grid%ny, &
                & dx => self%grid%dx, dy => self%grid%dy)
                
                !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
                !!PRINT*,'in continuity ',omp_get_thread_num()
                !!$OMP DO COLLAPSE(2)
                !!PRINT*, 'in continuity ',omp_get_thread_num() ! CAN NOT PUT COMMAND BETWEEN !$OMP DO AND LOOP
                ! Spatial loop
                y_loop: DO j=1,ny
                x_loop: DO i=1,nx
                    self%field%rho(i,j) = rho(i,j) &
                                       & -dt*( ( q(i+1,j,1)-q(i-1,j,1) )/( dx(i-1)+dx(i) ) &
                                       &      +( q(i,j+1,2)-q(i,j-1,2) )/( dy(j-1)+dy(j) )  )
                END DO x_loop
                END DO y_loop
                ! Spatial loop over
                !!$OMP END DO
                !!$OMP END PARALLEL
        END ASSOCIATE
        END SUBROUTINE solve_continuity


        SUBROUTINE solve_x_momentum(self, field_old)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            TYPE(field_type), POINTER :: field_old
            INTEGER :: i, j
            
            ! PRINT*,'in x-momentum ',omp_get_thread_num()
            ASSOCIATE(rho => field_old%rho, u => field_old%u, v => field_old%v, &
                & T => field_old%T, p => field_old%p, mu => field_old%mu, Pr => field_old%Pr, &
                & q => field_old%q, dt => self%dt, &
                & nx => field_old%grid%nx, ny => field_old%grid%ny, &
                & dx => self%grid%dx, dy => self%grid%dy)
                ! Spatial loop
                y_loop: DO j=1,ny
                x_loop: DO i=1,nx
                    ! (rho*u) + dt*( -d(rho*u*u + p)/dx -d(rho*v*u)/dy )
                    self%field%q(i,j,1) = q(i,j,1) + dt*( &
                        & -( (q(i+1,j,1)*u(i+1,j)+p(i+1,j))-(q(i-1,j,1)*u(i-1,j)+p(i-1,j)) )/( dx(i-1)+dx(i) ) &
                        & -( (q(i,j+1,2)*u(i,j+1))-(q(i,j-1,2)*u(i,j-1)) )/( dy(j-1)+dy(j) ) )
                    
                    ! +dt*(1/Re_ref)*d(tau_xx)/dx
                    self%field%q(i,j,1) = self%field%q(i,j,1) + dt*(1/self%Re_ref)*(             & 
                        &    (4.0/3.0)*( 0.5*(mu(i,j)+mu(i+1,j)) * (u(i+1,j)-u(i,j)) * dx(i-1)   &  
                        &               -0.5*(mu(i,j)+mu(i-1,j)) * (u(i,j)-u(i-1,j)) * dx(i)     &
                        &              )/( dx(i)*dx(i-1)*(dx(i)+dx(i-1))*0.5 )                   &
                        &   -(2.0/3.0)*( mu(i+1,j) * (v(i+1,j+1)-v(i+1,j-1))                     &
                        &               -mu(i-1,j) * (v(i-1,j+1)-v(i-1,j-1))                     &
                        &              )/( (dx(i)+dx(i-1))*(dy(j)+dy(j-1)) )                     &
                        & )

                    ! +dt*(1/Re_ref)*d(tau_xy)/dy
                    self%field%q(i,j,1) = self%field%q(i,j,1) + dt*(1/self%Re_ref)*(   &
                        &    ( 0.5*(mu(i,j)+mu(i,j+1)) * (u(i,j+1)-u(i,j)) * dy(j-1)   &  
                        &     -0.5*(mu(i,j)+mu(i,j-1)) * (u(i,j)-u(i,j-1)) * dy(j)     &
                        &    )/( dy(j)*dy(j-1)*(dy(j)+dy(j-1))*0.5 )                   &
                        &   +( mu(i,j+1) * (v(i+1,j+1)-v(i-1,j+1))                     &
                        &     -mu(i,j-1) * (v(i+1,j-1)-v(i-1,j-1))                     &
                        &    )/( (dy(j)+dy(j-1))*(dx(i)+dx(i-1)) )                     &
                        & )         
                END DO x_loop
                END DO y_loop
                ! Spatial loop over
            END ASSOCIATE
        END SUBROUTINE solve_x_momentum


        SUBROUTINE solve_y_momentum(self, field_old)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            TYPE(field_type), POINTER :: field_old
            INTEGER :: i, j
            
            ! PRINT*,'in y momentum ',omp_get_thread_num()
            ASSOCIATE(rho => field_old%rho, u => field_old%u, v => field_old%v, &
                & T => field_old%T, p => field_old%p, mu => field_old%mu, Pr => field_old%Pr, &
                & q => field_old%q, dt => self%dt, &
                & nx => field_old%grid%nx, ny => field_old%grid%ny, &
                & dx => self%grid%dx, dy => self%grid%dy)
                ! Spatial loop
                y_loop: DO j=1,ny
                x_loop: DO i=1,nx
                    ! (rho*v) + dt*( -d(rho*u*v)/dx -d(rho*v*v + p)/dy )
                    self%field%q(i,j,2) = q(i,j,2) + dt*( &
                        & -( (q(i+1,j,1)*v(i+1,j))-(q(i-1,j,1)*v(i-1,j)) )/( dx(i-1)+dx(i) ) &
                        & -( (q(i,j+1,2)*v(i,j+1)+p(i,j+1))-(q(i,j-1,2)*v(i,j-1)+p(i,j-1)) )/( dy(j-1)+dy(j) ) )
                    
                    ! +dt*(1/Re_ref)*d(tau_xy)/dx
                    self%field%q(i,j,2) = self%field%q(i,j,2) + dt*(1/self%Re_ref)*(   &
                        &    ( mu(i+1,j) * (u(i+1,j+1)-u(i+1,j-1))                     &
                        &     -mu(i-1,j) * (u(i-1,j+1)-u(i-1,j-1))                     &
                        &    )/( (dx(i)+dx(i-1))*(dy(j)+dy(j-1)) )                     &
                        &   +( 0.5*(mu(i,j)+mu(i+1,j)) * (v(i+1,j)-v(i,j)) * dx(i-1)   &  
                        &     -0.5*(mu(i,j)+mu(i-1,j)) * (v(i,j)-v(i-1,j)) * dx(i)     &
                        &    )/( dx(i)*dx(i-1)*(dx(i)+dx(i-1))*0.5 )                   &
                        & )         
                    
                    ! +dt*(1/Re_ref)*d(tau_yy)/dy
                    self%field%q(i,j,2) = self%field%q(i,j,2) + dt*(1/self%Re_ref)*(             & 
                        &    (4.0/3.0)*( 0.5*(mu(i,j)+mu(i,j+1)) * (v(i,j+1)-v(i,j)) * dy(j-1)   &  
                        &               -0.5*(mu(i,j)+mu(i,j-1)) * (v(i,j)-v(i,j-1)) * dy(j)     &
                        &              )/( dy(j)*dy(j-1)*(dy(j)+dy(j-1))*0.5 )                   &
                        &   -(2.0/3.0)*( mu(i,j+1) * (u(i+1,j+1)-u(i-1,j+1))                     &
                        &               -mu(i,j-1) * (u(i+1,j-1)-u(i-1,j-1))                     &
                        &              )/( (dx(i)+dx(i-1))*(dy(j)+dy(j-1)) )                     &
                        & )         
                END DO x_loop
                END DO y_loop
                ! Spatial loop over
            END ASSOCIATE
        END SUBROUTINE solve_y_momentum


        SUBROUTINE solve_energy(self, field_old)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            TYPE(field_type), POINTER :: field_old
            INTEGER :: i, j
            
            ! PRINT*,'in energy ',omp_get_thread_num()
            ASSOCIATE(rho => field_old%rho, u => field_old%u, v => field_old%v, &
                & T => field_old%T, p => field_old%p, mu => field_old%mu, Pr => field_old%Pr, &
                & q => field_old%q, dt => self%dt, &
                & nx => field_old%grid%nx, ny => field_old%grid%ny, &
                & dx => self%grid%dx, dy => self%grid%dy)
                ! Spatial loop
                y_loop: DO j=1,ny
                x_loop: DO i=1,nx
                    ! (rho*u*E) + dt*( -d(rho*u*E+p*u)/dx -d(rho*v*E + p*v)/dy )
                    self%field%q(i,j,3) = q(i,j,3) + dt*( &
                        & -( ( (q(i+1,j,3)+p(i+1,j))*u(i+1,j) )-( (q(i-1,j,3)+p(i-1,j))*u(i-1,j) ) )/( dx(i-1)+dx(i) ) &
                        & -( ( (q(i,j+1,3)+p(i,j+1))*v(i,j+1) )-( (q(i,j-1,3)+p(i,j-1))*v(i,j-1) ) )/( dy(j-1)+dy(j) ) )
                    
                    ! +dt*(1/Re_ref)*d(u*tau_xx)/dx
                    self%field%q(i,j,3) = self%field%q(i,j,3) + dt*(1/self%Re_ref)*( &
                        &    (4.0/3.0)*( 0.25*(mu(i,j)+mu(i+1,j)) * (u(i,j)+u(i+1,j)) * (u(i+1,j)-u(i,j)) * dx(i-1) &  
                        &               -0.25*(mu(i,j)+mu(i-1,j)) * (u(i,j)+u(i-1,j)) * (u(i,j)-u(i-1,j)) * dx(i)   &
                        &              )/( dx(i)*dx(i-1)*(dx(i)+dx(i-1))*0.5 )                                      &
                        &   -(2.0/3.0)*( mu(i+1,j) * u(i+1,j) * (v(i+1,j+1)-v(i+1,j-1))                             &
                        &               -mu(i-1,j) * u(i-1,j) * (v(i-1,j+1)-v(i-1,j-1))                             &
                        &              )/( (dx(i)+dx(i-1))*(dy(j)+dy(j-1)) )                                        &
                        & ) 

                    ! +dt*(1/Re_ref)*d(v*tau_xy)/dx
                    self%field%q(i,j,3) = self%field%q(i,j,3) + dt*(1/self%Re_ref)*( &
                        &    ( mu(i+1,j) * v(i+1,j) * (u(i+1,j+1)-u(i+1,j-1))                             &
                        &     -mu(i-1,j) * v(i-1,j) * (u(i-1,j+1)-u(i-1,j-1))                             &
                        &    )/( (dx(i)+dx(i-1))*(dy(j)+dy(j-1)) )                                        &
                        &   +( 0.25*(mu(i,j)+mu(i+1,j)) * (v(i,j)+v(i+1,j)) * (v(i+1,j)-v(i,j)) * dx(i-1) &  
                        &     -0.25*(mu(i,j)+mu(i-1,j)) * (v(i,j)+v(i-1,j)) * (v(i,j)-v(i-1,j)) * dx(i)   &
                        &    )/( dx(i)*dx(i-1)*(dx(i)+dx(i-1))*0.5 )                                      &
                        & )

                    ! +dt*(1/Re_ref)*d(q_x)/dx
                    self%field%q(i,j,3) = self%field%q(i,j,3) + dt*(1/self%Re_ref)*(1/(self%gamma-1)/(self%M_ref**2))*( &
                        &    ( 0.5*(mu(i,j)/Pr(i,j)+mu(i+1,j)/Pr(i+1,j)) * (T(i+1,j)-T(i,j)) * dx(i-1)    &
                        &     -0.5*(mu(i,j)/Pr(i,j)+mu(i-1,j)/Pr(i-1,j)) * (T(i,j)-T(i-1,j)) * dx(i)      &
                        &    )/( dx(i)*dx(i-1)*(dx(i)+dx(i-1))*0.5 )                                      &
                        & )
                    
                    ! +dt*(1/Re_ref)*d(u*tau_xy)/dy
                    self%field%q(i,j,3) = self%field%q(i,j,3) + dt*(1/self%Re_ref)*( &
                        &    ( 0.25*(mu(i,j)+mu(i,j+1)) * (u(i,j)+u(i,j+1)) * (u(i,j+1)-u(i,j)) * dy(j-1) &  
                        &     -0.25*(mu(i,j)+mu(i,j-1)) * (u(i,j)+u(i,j-1)) * (u(i,j)-u(i,j-1)) * dy(j)   &
                        &    )/( dy(j)*dy(j-1)*(dy(j)+dy(j-1))*0.5 )                                      & 
                        &   +( mu(i,j+1) * u(i,j+1) * (v(i+1,j+1)-v(i-1,j+1))                             &
                        &     -mu(i,j-1) * u(i,j-1) * (v(i+1,j-1)-v(i-1,j-1))                             &
                        &    )/( (dx(i)+dx(i-1))*(dy(j)+dy(j-1)) )                                        &
                        & )
                        
                    ! +dt*(1/Re_ref)*d(v*tau_yy)/dy
                    self%field%q(i,j,3) = self%field%q(i,j,3) + dt*(1/self%Re_ref)*( &
                        &    (4.0/3.0)*( 0.25*(mu(i,j)+mu(i,j+1)) * (v(i,j)+v(i,j+1)) * (v(i,j+1)-v(i,j)) * dy(j-1) &  
                        &               -0.25*(mu(i,j)+mu(i,j-1)) * (v(i,j)+v(i,j-1)) * (v(i,j)-v(i,j-1)) * dy(j)   &
                        &              )/( dy(j)*dy(j-1)*(dy(j)+dy(j-1))*0.5 )                                      &
                        &   -(2.0/3.0)*( mu(i,j+1) * v(i,j+1) * (u(i+1,j+1)-u(i-1,j+1))                             &
                        &               -mu(i,j-1) * v(i,j-1) * (u(i+1,j-1)-u(i-1,j-1))                             &
                        &              )/( (dx(i)+dx(i-1))*(dy(j)+dy(j-1)) )                                        &
                        & ) 
                        
                    ! +dt*(1/Re_ref)*d(q_y)/dy
                    self%field%q(i,j,3) = self%field%q(i,j,3) + dt*(1/self%Re_ref)*(1/(self%gamma-1)/(self%M_ref**2))*( &
                        &    ( 0.5*(mu(i,j)/Pr(i,j)+mu(i,j+1)/Pr(i,j+1)) * (T(i,j+1)-T(i,j)) * dy(j-1)    &
                        &     -0.5*(mu(i,j)/Pr(i,j)+mu(i,j-1)/Pr(i,j-1)) * (T(i,j)-T(i,j-1)) * dy(j)      &
                        &    )/( dy(j)*dy(j-1)*(dy(j)+dy(j-1))*0.5 )                                      &
                        & )
                END DO x_loop
                END DO y_loop
                ! Spatial loop over
            END ASSOCIATE
        END SUBROUTINE solve_energy


        SUBROUTINE add_numerical_diffusion(self, field_old)
            CLASS(solver_compressible), INTENT(INOUT) :: self
            TYPE(field_type), POINTER :: field_old
            INTEGER :: i, j
            
            !PRINT*,'in add_numerical_diff ',omp_get_thread_num()
            ASSOCIATE(rho => field_old%rho, u => field_old%u, v => field_old%v, &
                & T => field_old%T, p => field_old%p, mu => field_old%mu, Pr => field_old%Pr, &
                & q => field_old%q, dt => self%dt, &
                & nx => field_old%grid%nx, ny => field_old%grid%ny, &
                & dx => self%grid%dx, dy => self%grid%dy, B => self%grid%B, C => self%grid%C)
                ! Spatial loop: Internal nodes: fourth order diffusion term on the conserved variables
                ! -eps4*dx^4*d4f/dx4-eps4*dy^4*d4f/dy4, f = rho, rho*u, rho*v, rho*E
                y_loop: DO j=2,ny-1
                x_loop: DO i=2,nx-1
                    ! rho
                    self%field%rho(i,j) = self%field%rho(i,j)                                                       &
                        & -(dx(i)**4)*self%eps4*(                                                                   &
                        & rho(i+1,j)*B(i,1)+rho(i-1,j)*B(i,2)+rho(i+2,j)*B(i,3)+rho(i-2,j)*B(i,4)+rho(i,j)*B(i,5) ) &
                        & -(dy(j)**4)*self%eps4*(                                                                   &
                        & rho(i,j+1)*C(j,1)+rho(i,j-1)*C(j,2)+rho(i,j+2)*C(j,3)+rho(i,j-2)*C(j,4)+rho(i,j)*C(j,5) )

                    ! rho*u, rho*v, rho*E
                    self%field%q(i,j,:) = self%field%q(i,j,:)                                                       &
                        & -(dx(i)**4)*self%eps4*(                                                                   &
                        & q(i+1,j,:)*B(i,1)+q(i-1,j,:)*B(i,2)+q(i+2,j,:)*B(i,3)+q(i-2,j,:)*B(i,4)+q(i,j,:)*B(i,5) ) &
                        & -(dy(j)**4)*self%eps4*(                                                                   &
                        & q(i,j+1,:)*C(j,1)+q(i,j-1,:)*C(j,2)+q(i,j+2,:)*C(j,3)+q(i,j-2,:)*C(j,4)+q(i,j,:)*C(j,5) )
						
					! adding damping on p into rho*E
                    self%field%q(i,j,3) = self%field%q(i,j,3)                                                       &
                        & -(dx(i)**4)*self%eps4*(                                                                   &
                        & p(i+1,j)*B(i,1)+p(i-1,j)*B(i,2)+p(i+2,j)*B(i,3)+p(i-2,j)*B(i,4)+p(i,j)*B(i,5) ) &
                        & -(dy(j)**4)*self%eps4*(                                                                   &
                        & p(i,j+1)*C(j,1)+p(i,j-1)*C(j,2)+p(i,j+2)*C(j,3)+p(i,j-2)*C(j,4)+p(i,j)*C(j,5) )
                END DO x_loop
                END DO y_loop
                ! Spatial loop over

                ! Spatial loop: Domain boundary nodes: second order diffusion term on the conserved variables (only rho and rho*E)
                ! +eps2*dx^2*d2f/dx2+eps2*dy^2*d2f/dy2, f = rho, rho*E
                DO j=1,ny ! Left and right domain boundaries
                    ! rho
                    self%field%rho((/1,nx/),j) = self%field%rho((/1,nx/),j)                                              &
                    & +self%eps2*(dx((/1,nx/))**2)*( ( rho((/1,nx/)+1,j)-rho((/1,nx/),j) )*dx((/1,nx/)-1)                &
                    &                               -( rho((/1,nx/),j)-rho((/1,nx/)-1,j) )*dx((/1,nx/))                  &
                    &                              )/( dx((/1,nx/))*dx((/1,nx/)-1)*0.5*( dx((/1,nx/))+dx((/1,nx/)-1) ) ) &
                    & +self%eps2*(dy(j)**2)*( ( rho((/1,nx/),j+1)-rho((/1,nx/),j) )*dy(j-1)                              &
                    &                        -( rho((/1,nx/),j)-rho((/1,nx/),j-1) )*dy(j)                                &
                    &                       )/( dy(j)*dy(j-1)*0.5*( dy(j)+dy(j-1) ) ) 

                    ! rho*E
                    self%field%q((/1,nx/),j,3) = self%field%q((/1,nx/),j,3)                                              &
                    & +self%eps2*(dx((/1,nx/))**2)*( ( q((/1,nx/)+1,j,3)-q((/1,nx/),j,3) )*dx((/1,nx/)-1)                &
                    &                               -( q((/1,nx/),j,3)-q((/1,nx/)-1,j,3) )*dx((/1,nx/))                  &
                    &                              )/( dx((/1,nx/))*dx((/1,nx/)-1)*0.5*( dx((/1,nx/))+dx((/1,nx/)-1) ) ) &
                    & +self%eps2*(dy(j)**2)*( ( q((/1,nx/),j+1,3)-q((/1,nx/),j,3) )*dy(j-1)                              &
                    &                        -( q((/1,nx/),j,3)-q((/1,nx/),j-1,3) )*dy(j)                                &
                    &                       )/( dy(j)*dy(j-1)*0.5*( dy(j)+dy(j-1) ) ) 

                    ! adding damping on p into rho*E
                    self%field%q((/1,nx/),j,3) = self%field%q((/1,nx/),j,3)                                              &
                    & +self%eps2*(dx((/1,nx/))**2)*( ( p((/1,nx/)+1,j)-p((/1,nx/),j) )*dx((/1,nx/)-1)                    &
                    &                               -( p((/1,nx/),j)-p((/1,nx/)-1,j) )*dx((/1,nx/))                      &
                    &                              )/( dx((/1,nx/))*dx((/1,nx/)-1)*0.5*( dx((/1,nx/))+dx((/1,nx/)-1) ) ) &
                    & +self%eps2*(dy(j)**2)*( ( p((/1,nx/),j+1)-p((/1,nx/),j) )*dy(j-1)                                  &
                    &                        -( p((/1,nx/),j)-p((/1,nx/),j-1) )*dy(j)                                    &
                    &                       )/( dy(j)*dy(j-1)*0.5*( dy(j)+dy(j-1) ) ) 
                END DO

                ! Top and bottom domain boundaries
                DO i=1,nx ! Left and right domain boundaries
                    ! rho
                    self%field%rho(i,(/1,ny/)) = self%field%rho(i,(/1,ny/))                                              &
                    & +self%eps2*(dx(i)**2)*( ( rho(i+1,(/1,ny/))-rho(i,(/1,ny/)) )*dx(i-1)                              &
                    &                        -( rho(i,(/1,ny/))-rho(i-1,(/1,ny/)) )*dx(i)                                &
                    &                       )/( dx(i)*dx(i-1)*0.5*( dx(i)+dx(i-1) ) )                                    &
                    & +self%eps2*(dy((/1,ny/))**2)*( ( rho(i,(/1,ny/)+1)-rho(i,(/1,ny/)) )*dy((/1,ny/)-1)                &
                    &                               -( rho(i,(/1,ny/))-rho(i,(/1,ny/)-1) )*dy((/1,ny/))                  &
                    &                              )/( dy((/1,ny/))*dy((/1,ny/)-1)*0.5*( dy((/1,ny/))+dy((/1,ny/)-1) ) ) 

                    ! rho*E
                    self%field%q(i,(/1,ny/),3) = self%field%q(i,(/1,ny/),3)                                              &
                    & +self%eps2*(dx(i)**2)*( ( q(i+1,(/1,ny/),3)-q(i,(/1,ny/),3) )*dx(i-1)                              &
                    &                        -( q(i,(/1,ny/),3)-q(i-1,(/1,ny/),3) )*dx(i)                                &
                    &                       )/( dx(i)*dx(i-1)*0.5*( dx(i)+dx(i-1) ) )                                    &
                    & +self%eps2*(dy((/1,ny/))**2)*( ( q(i,(/1,ny/)+1,3)-q(i,(/1,ny/),3) )*dy((/1,ny/)-1)                &
                    &                               -( q(i,(/1,ny/),3)-q(i,(/1,ny/)-1,3) )*dy((/1,ny/))                  &
                    &                              )/( dy((/1,ny/))*dy((/1,ny/)-1)*0.5*( dy((/1,ny/))+dy((/1,ny/)-1) ) )

                    ! adding damping on p into rho*E
                    self%field%q(i,(/1,ny/),3) = self%field%q(i,(/1,ny/),3)                                              &
                    & +self%eps2*(dx(i)**2)*( ( p(i+1,(/1,ny/))-p(i,(/1,ny/)) )*dx(i-1)                                  &
                    &                        -( p(i,(/1,ny/))-p(i-1,(/1,ny/)) )*dx(i)                                    &
                    &                       )/( dx(i)*dx(i-1)*0.5*( dx(i)+dx(i-1) ) )                                    &
                    & +self%eps2*(dy((/1,ny/))**2)*( ( p(i,(/1,ny/)+1)-p(i,(/1,ny/)) )*dy((/1,ny/)-1)                    &
                    &                               -( p(i,(/1,ny/))-p(i,(/1,ny/)-1) )*dy((/1,ny/))                      &
                    &                              )/( dy((/1,ny/))*dy((/1,ny/)-1)*0.5*( dy((/1,ny/))+dy((/1,ny/)-1) ) )
                END DO
            END ASSOCIATE
        END SUBROUTINE add_numerical_diffusion


        SUBROUTINE conserved_weighted_sum(self, w1, field1, w2, field2)
            ! This subroutine is to implement RK3 time integration
            CLASS(solver_compressible), INTENT(INOUT) :: self
            TYPE(field_type), POINTER, INTENT(IN) :: field1, field2
            REAL, INTENT(IN) :: w1, w2

            ! Assiigning conserved field variables rho, rho*u, rho*v, rho*E
            self%field%rho = w1*field1%rho + w2*field2%rho
            self%field%q = w1*field1%q + w2*field2%q

            ! Obtaining primitive variables u, v, T, p, mu 
            CALL self%primitive_from_conserved()
            CALL self%field%set_bc()
            CALL self%compute_p_mu()

            ! Obtaining conserved variables at domain boundaries and ghost nodes
            CALL self%conserved_from_primitive('boundary_and_ghosts')
        END SUBROUTINE conserved_weighted_sum

END MODULE mod_solver