
MODULE mod_field
    USE mod_parameters
    USE mod_utility
    USE mod_grid
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : error_unit
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: field_type
    PUBLIC :: print_rms_diff_fields
    public :: copy_field
    
    ! PUBLIC :: ASSIGNMENT(=) ! With ifort (intel) compiler, assignment operator can not be defined on allocatable or pointers. gfortran accepts it. 
    ! ! Interface block for assignment operator overloading
    ! INTERFACE ASSIGNMENT(=)
    !   MODULE PROCEDURE copy_field
    ! END INTERFACE

    TYPE :: field_type
        REAL(KIND=rkind), ALLOCATABLE, DIMENSION(:,:)   :: rho, u, v, p, T, mu, Pr ! non-dim 
        REAL(KIND=rkind), ALLOCATABLE, DIMENSION(:,:,:) :: q ! Non-dim conserved variables [rho*u, rho*v, rho*E], E = total energy                                                  
        CLASS(grid_type), POINTER :: grid => NULL()
        
        CONTAINS
            PROCEDURE, PUBLIC,  PASS(self) :: construct
            PROCEDURE, PUBLIC,  PASS(self) :: printf
            PROCEDURE, PUBLIC,  PASS(self) :: set_bc
    END TYPE field_type


    CONTAINS
        ! Procedures of user defined type: field_type
        SUBROUTINE construct(self, grid)
            CLASS(field_type), INTENT(INOUT) :: self
            TYPE(grid_type), POINTER :: grid

            self%grid => grid ! Pointing to grid (associating to grid)

            ! Allocating field variables with collocated grid and ghost nodes
            IF(.NOT. ALLOCATED(self%rho)) ALLOCATE( self%rho(0:self%grid%nx+1, 0:self%grid%ny+1) )
            IF(.NOT. ALLOCATED(self%u))   ALLOCATE( self%u(0:self%grid%nx+1,   0:self%grid%ny+1) )
            IF(.NOT. ALLOCATED(self%v))   ALLOCATE( self%v(0:self%grid%nx+1,   0:self%grid%ny+1) )
            IF(.NOT. ALLOCATED(self%p))   ALLOCATE( self%p(0:self%grid%nx+1,   0:self%grid%ny+1) )
            IF(.NOT. ALLOCATED(self%T))   ALLOCATE( self%T(0:self%grid%nx+1,   0:self%grid%ny+1) )
            IF(.NOT. ALLOCATED(self%mu))  ALLOCATE( self%mu(0:self%grid%nx+1,  0:self%grid%ny+1) )
            IF(.NOT. ALLOCATED(self%Pr))  ALLOCATE( self%Pr(0:self%grid%nx+1,  0:self%grid%ny+1) )
            IF(.NOT. ALLOCATED(self%q))   ALLOCATE( self%q(0:self%grid%nx+1,   0:self%grid%ny+1, 1:3) )

            ! Initializing field variables with zeros
            self%rho = 1.0
            self%u   = 0.0
            self%v   = 0.0
            self%T   = 1.0
            self%p   = 0.0 ! Will be assigned in solver after setting BC
            self%mu  = 0.0 ! Will be assigned in solver after setting BC
            self%Pr  = 0.0 ! Will be assigned in solver after setting BC
            self%q   = 0.0 ! Will be assigned in solver
            
            ! setting BC
            CALL self%set_bc()

            CALL self%printf()
        END SUBROUTINE construct


        SUBROUTINE set_bc(self)
            ! BCs are applied to the field variables: rho, u, v, and T
            CLASS(field_type), INTENT(INOUT) :: self
            ASSOCIATE(rho => self%rho, u => self%u, v => self%v, T => self%T, &
                    & nx => self%grid%nx, ny => self%grid%ny,&
                    & dx => self%grid%dx, dy => self%grid%dy)
                ! No slip condition on u and v
                    u(1,1:ny)  = 0.0         ! left boundary
                    u(nx,1:ny) = 0.0         ! right boundary
                    u(1:nx,1)  = 0.0         ! bottom boundary
                    u(1:nx,ny) = 1.0         ! top boundary (lid velocity (non-dim)
                    v(1,1:ny)  = 0.0         ! left boundary
                    v(nx,1:ny) = 0.0         ! right boundary
                    v(1:nx,1)  = 0.0         ! bottom boundary
                    v(1:nx,ny) = 0.0         ! top boundary (lid velocity (non-dim)

                ! Linear extrapolatoiion for u and v on Ghost node (u,v at corner ghost nodes are not needed)
                    u(0,1:ny)     = 2*u(1,1:ny)-u(2,1:ny)     ! left ghost nodes
                    v(0,1:ny)     = 2*v(1,1:ny)-v(2,1:ny)
                    u(nx+1,1:ny)  = 2*u(nx,1:ny)-u(nx-1,1:ny) ! right ghost nodes
                    v(nx+1,1:ny)  = 2*v(nx,1:ny)-v(nx-1,1:ny)
                    u(1:nx,0)     = 2*u(1:nx,1)-u(1:nx,2)     ! bottom ghost nodes
                    v(1:nx,0)     = 2*v(1:nx,1)-v(1:nx,2)
                    u(1:nx,ny+1)  = 2*u(1:nx,ny)-u(1:nx,ny-1) ! top ghost nodes
                    v(1:nx,ny+1)  = 2*v(1:nx,ny)-v(1:nx,ny-1)

                ! BC on rho and T using ghost nodes (zero normal gradient) (p is calculated in the solver using eqn. of state)
                    rho(0,1:ny)    = rho(2,1:ny)    ! left boundary ghost nodes
                    T(0,1:ny)      = T(2,1:ny)
                    rho(nx+1,1:ny) = rho(nx-1,1:ny) ! right boundary ghost nodes
                    T(nx+1,1:ny)   = T(nx-1,1:ny)
                    rho(1:nx,0)    = rho(1:nx,2)    ! bottom boundary ghost nodes
                    T(1:nx,0)      = T(1:nx,2)
                    rho(1:nx,ny+1) = rho(1:nx,ny-1) ! top boundary ghost nodes
                    T(1:nx,ny+1)   = T(1:nx,ny-1)

                ! Averaging rho and T for corner nodes of the domain
                    rho(1,1)   = ( rho(2,1)*dy(1)+rho(1,2)*dx(1) )/(dy(1)+dx(1))
                    rho(nx,1)  = ( rho(nx-1,1)*dy(1)+rho(nx,2)*dx(nx-1) )/(dy(1)+dx(nx-1))
                    rho(nx,ny) = ( rho(nx-1,ny)*dy(ny-1)+rho(nx,ny-1)*dx(nx-1) )/(dy(ny-1)+dx(nx-1))
                    rho(1,ny)  = ( rho(2,ny)*dy(ny-1)+rho(1,ny-1)*dx(1) )/(dy(ny-1)+dx(1))
                    T(1,1)     = ( T(2,1)*dy(1)+T(1,2)*dx(1) )/(dy(1)+dx(1))
                    T(nx,1)    = ( T(nx-1,1)*dy(1)+T(nx,2)*dx(nx-1) )/(dy(1)+dx(nx-1))
                    T(nx,ny)   = ( T(nx-1,ny)*dy(ny-1)+T(nx,ny-1)*dx(nx-1) )/(dy(ny-1)+dx(nx-1))
                    T(1,ny)    = ( T(2,ny)*dy(ny-1)+T(1,ny-1)*dx(1) )/(dy(ny-1)+dx(1))
            END ASSOCIATE
        END SUBROUTINE set_bc


        SUBROUTINE printf(self, format)
            ! Prints field variables value of a 2D grid in matrix form
            CLASS(field_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: format ! eg. format='F7.3'
            INTEGER :: ix_low, ix_up, iy_low, iy_up

            ASSOCIATE(rho => self%rho, u => self%u, v => self%v, p => self%p, T => self%T)
                ! Lower and upper bounds for each field variable are identical due to collocated grid.
                ix_low = LBOUND(rho,1)+1; ix_up  = UBOUND(rho,1)-1;
                iy_low = LBOUND(rho,2)+1; iy_up  = UBOUND(rho,2)-1;
                
                ! rho
                PRINT*, 'field: rho, shape: (dim1, dim2) = (',SHAPE(rho(ix_low:ix_up,iy_low:iy_up)),')'
                IF(PRESENT(format)) THEN
                    CALL print_data(data=rho(ix_low:ix_up,iy_low:iy_up), format=format, dim_along_row = 1)
                ELSE
                    CALL print_data(data=rho(ix_low:ix_up,iy_low:iy_up), dim_along_row = 1)
                END IF

                ! u velocity
                PRINT*, 'field: u, shape: (dim1, dim2) = (',SHAPE(u(ix_low:ix_up,iy_low:iy_up)),')'
                IF(PRESENT(format)) THEN
                    CALL print_data(data=u(ix_low:ix_up,iy_low:iy_up), format=format, dim_along_row = 1)
                ELSE
                    CALL print_data(data=u(ix_low:ix_up,iy_low:iy_up), dim_along_row = 1)
                END IF

                ! v velocity
                PRINT*, 'field: v, shape: (dim1, dim2) = (',SHAPE(v(ix_low:ix_up,iy_low:iy_up)),')'
                IF(PRESENT(format)) THEN
                    CALL print_data(data=v(ix_low:ix_up,iy_low:iy_up), format=format, dim_along_row = 1)
                ELSE
                    CALL print_data(data=v(ix_low:ix_up,iy_low:iy_up), dim_along_row = 1)
                END IF

                ! p pressure
                PRINT*, 'field: p, shape: (dim1, dim2) = (',SHAPE(p(ix_low:ix_up,iy_low:iy_up)),')'
                IF(PRESENT(format)) THEN
                    CALL print_data(data=p(ix_low:ix_up,iy_low:iy_up), format=format, dim_along_row = 1)
                ELSE
                    CALL print_data(data=p(ix_low:ix_up,iy_low:iy_up), dim_along_row = 1)
                END IF

                ! T pressure
                PRINT*, 'field: T, shape: (dim1, dim2) = (',SHAPE(T(ix_low:ix_up,iy_low:iy_up)),')'
                IF(PRESENT(format)) THEN
                    CALL print_data(data=T(ix_low:ix_up,iy_low:iy_up), format=format, dim_along_row = 1)
                ELSE
                    CALL print_data(data=T(ix_low:ix_up,iy_low:iy_up), dim_along_row = 1)
                END IF
            END ASSOCIATE
        END SUBROUTINE printf


        ! Procedures of module mod_field
        SUBROUTINE copy_field(left, right) ! Instead of copying, point to the old field pointer, if possible, for efficiency.
            TYPE(field_type), INTENT(OUT), POINTER  :: left
            TYPE(field_type), INTENT(IN),  POINTER  :: right
            INTEGER :: lx, ux, ly, uy
    
            IF(ASSOCIATED(left)) DEALLOCATE(left) ! Necessary to avoid memory leak
            ALLOCATE(left) ! Allocating afresh

            ! Lower and upper bounds; For collocated grid, the bounds are identical for all the field variables.
            lx = LBOUND(right%rho,1); ux = UBOUND(right%rho,1);
            ly = LBOUND(right%rho,2); uy = UBOUND(right%rho,2);
    
            ALLOCATE( left%rho(lx:ux,ly:uy) )
            ALLOCATE( left%u(lx:ux,ly:uy) )
            ALLOCATE( left%v(lx:ux,ly:uy) )
            ALLOCATE( left%T(lx:ux,ly:uy) )
            ALLOCATE( left%p(lx:ux,ly:uy) )
            ALLOCATE( left%mu(lx:ux,ly:uy) )
            ALLOCATE( left%Pr(lx:ux,ly:uy) )
            ALLOCATE( left%q(lx:ux,ly:uy,SIZE(right%q,3)) )
            left%rho  = right%rho
            left%u    = right%u
            left%v    = right%v
            left%T    = right%T
            left%p    = right%p
            left%mu   = right%mu
            left%Pr   = right%Pr
            left%q    = right%q
            left%grid => right%grid ! Associating to the same grid memory
        END SUBROUTINE copy_field


        SUBROUTINE print_rms_diff_fields(field1,field2,format)
            TYPE(field_type), INTENT(OUT), POINTER  :: field1
            TYPE(field_type), INTENT(IN),  POINTER  :: field2
            CHARACTER(LEN=*), INTENT(IN) :: format ! eg. format='F7.3'
            INTEGER :: ix_low, ix_up, iy_low, iy_up
            REAL(KIND=rkind) :: rho_rms, u_rms , v_rms, p_rms

            ! Lower and upper bounds; For collocated grid, the bounds are identical for all the field variables.
            ix_low = LBOUND(field1%u,1)+1; ix_up  = UBOUND(field1%u,1)-1;
            iy_low = LBOUND(field1%u,2)+1; iy_up  = UBOUND(field1%u,2)-1;
            
            rho_rms = SQRT( SUM( ( field1%rho(ix_low:ix_up, iy_low:iy_up) &
                               &-field2%rho(ix_low:ix_up, iy_low:iy_up)   &
                               )**2                                       &
                             )/SIZE(field1%rho(ix_low:ix_up, iy_low:iy_up)) )


            u_rms = SQRT( SUM( ( field1%u(ix_low:ix_up, iy_low:iy_up) &
                               &-field2%u(ix_low:ix_up, iy_low:iy_up) &
                               )**2                                   &
                             )/SIZE(field1%u(ix_low:ix_up, iy_low:iy_up)) )

            
            v_rms = SQRT( SUM( ( field1%v(ix_low:ix_up, iy_low:iy_up) &
                               &-field2%v(ix_low:ix_up, iy_low:iy_up) &
                               )**2                                   &
                             )/SIZE(field1%v(ix_low:ix_up, iy_low:iy_up)) )

            
            p_rms = SQRT( SUM( ( field1%p(ix_low:ix_up, iy_low:iy_up) &
                               &-field2%p(ix_low:ix_up, iy_low:iy_up) &
                               )**2                                   &
                             )/SIZE(field1%p(ix_low:ix_up, iy_low:iy_up)) )

            WRITE(*, FMT='(*('//format//' '//'))', ADVANCE='NO') rho_rms, u_rms, v_rms, p_rms
        END SUBROUTINE print_rms_diff_fields

END MODULE mod_field