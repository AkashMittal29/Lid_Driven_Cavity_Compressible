
MODULE mod_grid
    USE mod_parameters
    USE mod_utility
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : error_unit
    IMPLICIT NONE
    PRIVATE ! by default: private
    PUBLIC :: grid_type


    TYPE :: grid_type
        REAL(KIND=rkind), ALLOCATABLE, DIMENSION(:,:)   :: x, y   ! node location matrices
        REAL(KIND=rkind), ALLOCATABLE, DIMENSION(:)     :: dx, dy ! grid spacing matrices
        INTEGER :: nx, ny ! Number of nodes in x and y direction (without ghost nodes)
        REAL(KIND=rkind), ALLOCATABLE, DIMENSION(:,:)   :: B, C ! For numerical diffusion 
        CONTAINS   
            PROCEDURE, PUBLIC,  PASS(self) :: initialize
            PROCEDURE, PRIVATE, PASS(self) :: read_grid_file
            PROCEDURE, PUBLIC,  PASS(self) :: printf
            PROCEDURE, PUBLIC,  PASS(self) :: numerical_diff_const
    END TYPE grid_type


    CONTAINS
        ! Procedures of user defined type: grid_type
        SUBROUTINE initialize(self,file)
            CLASS(grid_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: file

            CALL read_grid_file(self,file)
            CALL self%printf(format='F13.8')

            ! Allocating and defining dx dy with ghost nodes 
            ALLOCATE(self%dx(0:self%nx), self%dy(0:self%ny))
            ASSOCIATE(nx => self%nx, ny => self%ny, x => self%x, y => self%y, &
                      dx => self%dx, dy => self%dy)
                dx(1:nx-1) = x(2:nx,1)-x(1:nx-1,1)
                dx((/0,nx/)) = (/dx(1), dx(nx-1)/)
                dy(1:ny-1) = y(1,2:ny)-y(1,1:ny-1)
                dy((/0,ny/)) = (/dy(1), dy(ny-1)/)

                ! Assigning B1, B2, ..., B5, C1, C2, ..., C5 for numerical diffusion
                  CALL self%numerical_diff_const()
            END ASSOCIATE
            PRINT*; WRITE(*,'("dx: ")')
            CALL print_data(RESHAPE(self%dx, shape=(/1,SIZE(self%dx)/)), dim_along_row = 2, format='F13.6')
            PRINT*; WRITE(*,'("dy: ")')
            CALL print_data(RESHAPE(self%dy, shape=(/1,SIZE(self%dy)/)), dim_along_row = 2, format='F13.6')
        END SUBROUTINE initialize


        SUBROUTINE read_grid_file(self,file)
            CLASS(grid_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: file
            INTEGER :: funit, io_status, i, j, id
            CHARACTER(LEN=100) :: text, text1, text2
            INTEGER :: nz_dummy
            REAL ::  z_dummy ! for 2d z_dummy is not used

            ! opening file in reading mode
            PRINT*, TRIM(ADJUSTL(file))
            OPEN(NEWUNIT=funit, FILE=TRIM(ADJUSTL(file)), ACTION='read', &
                IOSTAT=io_status, STATUS='old', POSITION='REWIND')
            IF(io_status/=0) THEN
                WRITE(error_unit,*) 'Unable to open file.'
                STOP 'stopping the execution' ! stops the code execution
            END IF

            READ(funit,*)
            DO i=1,3,1 ! Reading nx, ny, and nz
                READ(funit,'(A)',IOSTAT=io_status) text
                text1 = TRIM(ADJUSTL(text(1:index(text,'=')-1)))
                text2 = TRIM(ADJUSTL(text(index(text,'=')+1:))) 
                SELECT CASE (trim(adjustl(text1)))
                    CASE ('nx'); READ(text2,*) self%nx
                    CASE ('ny'); READ(text2,*) self%ny
                    CASE ('nz'); READ(text2,*) nz_dummy     
                END SELECT
            END DO

            ALLOCATE(self%x(self%nx,self%ny), self%y(self%nx,self%ny))
            self%x = 0.0; self%y=0.0;

            DO j=1,self%ny,1
                DO i=1,self%nx,1
                    READ(funit,*) id, self%x(i,j), self%y(i,j), z_dummy
                END DO
            END DO
        END SUBROUTINE read_grid_file


        SUBROUTINE printf(self, format)
            ! Prints x and y values of a 2D grid in matrix form
            CLASS(grid_type), INTENT(INOUT) :: self ! Here self object is a dummy argument (i.e. passed by reference), hence CLASS can be used.
            CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: format ! eg. format='F7.3'

            ASSOCIATE(x => self%x, y => self%y)
                IF(PRESENT(format)) THEN
                    PRINT*, 'grid: x, shape: (nx, ny) = (',SHAPE(x),')'
                    CALL print_data(data=x, format=format, dim_along_row = 1)
                    PRINT*, 'grid: y, shape: (nx, ny) = (',SHAPE(y),')'
                    CALL print_data(data=y, format=format, dim_along_row = 1)
                ELSE
                    PRINT*, 'grid: x, shape: (nx, ny) = (',SHAPE(x),')'
                    CALL print_data(data=x, dim_along_row = 1)
                    PRINT*, 'grid: y, shape: (nx, ny) = (',SHAPE(y),')'
                    CALL print_data(data=y, dim_along_row = 1)
                END IF
            END ASSOCIATE
        END SUBROUTINE printf


        SUBROUTINE numerical_diff_const(self)
            ! Assigns B1, B2, ..., B5, C1, C2, ..., C5 for numerical diffusion for structured grid
            CLASS(grid_type), INTENT(INOUT) :: self
            REAL(KIND=rkind) :: a, b, c, d, D_temp
            INTEGER :: i

            ALLOCATE(self%B(0:self%nx-1,5), &
                   & self%C(0:self%ny-1,5)) ! Ghost nodes included only to keep indexing consistent.
            DO i=2,self%nx-1
                a = self%dx(i);
                b = self%dx(i-1);
                c = self%dx(i)+self%dx(i+1);
                d = self%dx(i-1)+self%dx(i-2);
                D_temp = determinant(RESHAPE((/a, (a**2)/2,  (a**3)/6, (a**4)/24, &
                                            & -b, (b**2)/2, -(b**3)/6, (b**4)/24, &
                                            &  c, (c**2)/2,  (c**3)/6, (c**4)/24, &
                                            & -d, (d**2)/2, -(d**3)/6, (d**4)/24 /), SHAPE=(/4,4/), ORDER=(/2,1/)));
                self%B(i,1) =-determinant(RESHAPE((/-b, (b**2)/2, -(b**3)/6,  c, (c**2)/2,  (c**3)/6,  -d, (d**2)/2, -(d**3)/6/), SHAPE=(/3,3/), ORDER=(/2,1/)));
                self%B(i,2) = determinant(RESHAPE((/ a, (a**2)/2,  (a**3)/6,  c, (c**2)/2,  (c**3)/6,  -d, (d**2)/2, -(d**3)/6/), SHAPE=(/3,3/), ORDER=(/2,1/)));
                self%B(i,3) =-determinant(RESHAPE((/ a, (a**2)/2,  (a**3)/6, -b, (b**2)/2, -(b**3)/6,  -d, (d**2)/2, -(d**3)/6/), SHAPE=(/3,3/), ORDER=(/2,1/)));
                self%B(i,4) = determinant(RESHAPE((/ a, (a**2)/2,  (a**3)/6, -b, (b**2)/2, -(b**3)/6,   c, (c**2)/2,  (c**3)/6/), SHAPE=(/3,3/), ORDER=(/2,1/)));
                self%B(i,5) =-SUM(self%B(i,1:4))
                self%B(i,:) = self%B(i,:)/D_temp 
            END DO
            DO i=2,self%ny-1
                a = self%dy(i);
                b = self%dy(i-1);
                c = self%dy(i)+self%dy(i+1);
                d = self%dy(i-1)+self%dy(i-2);
                D_temp = determinant(RESHAPE((/a, (a**2)/2,  (a**3)/6, (a**4)/24, &
                                            & -b, (b**2)/2, -(b**3)/6, (b**4)/24, &
                                            &  c, (c**2)/2,  (c**3)/6, (c**4)/24, &
                                            & -d, (d**2)/2, -(d**3)/6, (d**4)/24 /), SHAPE=(/4,4/), ORDER=(/2,1/)));
                self%C(i,1) =-determinant(RESHAPE((/-b, (b**2)/2, -(b**3)/6,  c, (c**2)/2,  (c**3)/6,  -d, (d**2)/2, -(d**3)/6/), SHAPE=(/3,3/), ORDER=(/2,1/)));
                self%C(i,2) = determinant(RESHAPE((/ a, (a**2)/2,  (a**3)/6,  c, (c**2)/2,  (c**3)/6,  -d, (d**2)/2, -(d**3)/6/), SHAPE=(/3,3/), ORDER=(/2,1/)));
                self%C(i,3) =-determinant(RESHAPE((/ a, (a**2)/2,  (a**3)/6, -b, (b**2)/2, -(b**3)/6,  -d, (d**2)/2, -(d**3)/6/), SHAPE=(/3,3/), ORDER=(/2,1/)));
                self%C(i,4) = determinant(RESHAPE((/ a, (a**2)/2,  (a**3)/6, -b, (b**2)/2, -(b**3)/6,   c, (c**2)/2,  (c**3)/6/), SHAPE=(/3,3/), ORDER=(/2,1/)));
                self%C(i,5) =-SUM(self%C(i,1:4))
                self%C(i,:) = self%C(i,:)/D_temp 
            END DO
            ! PRINT*,'B:'
            ! DO i=1,5
            !     WRITE(*,'(*(E10.3, " "))') self%B(2:,i)
            ! END DO
            ! PRINT*,'C:'
            ! DO i=1,5
            !     WRITE(*,'(*(E10.3, " "))') self%C(2:,i)
            ! END DO
        END SUBROUTINE numerical_diff_const
END MODULE mod_grid
