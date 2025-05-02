MODULE mod_read_file
    !  ! Usage example
    !  USE mod_read_file
    !  TYPE(input_data_type) :: input_data
    !  REAL(KIND=rkind) :: a ! Will have to update rkind consistent with the KIND of variable
    !  
    !  CALL input_data%read_input('input_file_name.ini')
    !  CALL input_data%assign_variable('section name','variable_name',a)

    USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
    USE mod_parameters, ONLY : rkind 
    IMPLICIT NONE
    !INTEGER, PARAMETER, PUBLIC :: rkind = SELECTED_REAL_KIND(P=15, R=300) ! double


    PRIVATE
    PUBLIC :: input_data_type


    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE assign_section, assign_var
    END INTERFACE
    

    TYPE :: input_data_type
        TYPE(section), ALLOCATABLE :: sections(:)

        CONTAINS
            PROCEDURE, PASS(self), PUBLIC :: add_section
            PROCEDURE, PASS(self), PUBLIC :: add_variable
            PROCEDURE, PASS(self), PUBLIC :: print_input_data
            PROCEDURE, PASS(self), PUBLIC :: get_index
            PROCEDURE, PASS(self), PUBLIC :: read_input
            PROCEDURE, PASS(self), PUBLIC :: assign_variable_real
            PROCEDURE, PASS(self), PUBLIC :: assign_variable_int
            PROCEDURE, PASS(self), PUBLIC :: assign_variable_char
            GENERIC :: assign_variable => assign_variable_real, assign_variable_int, assign_variable_char
    END TYPE


    TYPE :: section
        CHARACTER(LEN=:), ALLOCATABLE :: sect
        TYPE(variable),   ALLOCATABLE :: vars(:)

        CONTAINS
            PROCEDURE, PASS(self), PUBLIC :: print_section
            PROCEDURE, PASS(self), PUBLIC :: adding_variable
    END TYPE section


    TYPE :: variable
        CHARACTER(LEN=:), ALLOCATABLE :: var
        CHARACTER(LEN=:), ALLOCATABLE :: value
    END TYPE variable


    CONTAINS
        SUBROUTINE add_section(self,  sect_name)
            CLASS(input_data_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: sect_name
            TYPE(section), ALLOCATABLE :: sections_temp(:)
            INTEGER :: i

            IF(.NOT. ALLOCATED(self%sections)) THEN
                ALLOCATE(self%sections(1))
                self%sections(1)%sect = sect_name
                !PRINT*, self%sections(1)%sect, LEN(self%sections(1)%sect)
                RETURN
            END IF

            ALLOCATE(sections_temp(SIZE(self%sections)))
            sections_temp=self%sections

            ! Adding new section
            DEALLOCATE(self%sections)
            ALLOCATE(self%sections(SIZE(sections_temp)+1))
            DO i=1,SIZE(sections_temp)
                self%sections(i)=sections_temp(i)
            END DO
            i = SIZE(sections_temp)+1
            self%sections(i)%sect=sect_name

            DEALLOCATE(sections_temp)
        END SUBROUTINE add_section


        SUBROUTINE add_variable(self,  sect_name, var_name, var_value)
            CLASS(input_data_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: sect_name, var_name, var_value
            INTEGER :: ind

            DO ind=1,SIZE(self%sections)
                IF(self%sections(ind)%sect==sect_name) EXIT
            END DO
            CALL self%sections(ind)%adding_variable(var_name, var_value)
        END SUBROUTINE add_variable


        SUBROUTINE adding_variable(self, var_name, var_value)
            CLASS(section), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: var_name, var_value
            TYPE(variable), ALLOCATABLE :: vars_temp(:)
            INTEGER :: i

            IF(.NOT. ALLOCATED(self%vars)) THEN
                ALLOCATE(self%vars(1))
                self%vars(1)%var = var_name
                self%vars(1)%value = var_value
                RETURN
            END IF

            ALLOCATE(vars_temp(SIZE(self%vars)))
            vars_temp=self%vars

            ! Adding new variable
            DEALLOCATE(self%vars)
            ALLOCATE(self%vars(SIZE(vars_temp)+1))
            DO i=1,SIZE(vars_temp)
                self%vars(i)=vars_temp(i)
            END DO
            i = SIZE(vars_temp)+1
            self%vars(i)%var=var_name
            self%vars(i)%value=var_value

            DEALLOCATE(vars_temp)
        END SUBROUTINE adding_variable


        SUBROUTINE assign_section(left, right)
            TYPE(section), INTENT(OUT), ALLOCATABLE  :: left(:)
            TYPE(section), INTENT(IN)    :: right(:)
            INTEGER :: i
            
            IF(ALLOCATED(left)) DEALLOCATE(left)
            ALLOCATE(left(SIZE(right)))

            DO i=1,SIZE(right)
                left(i)%sect = right(i)%sect
                ALLOCATE(left(i)%vars(SIZE(right(i)%vars)))
                left(i)%vars=right(i)%vars
            END DO
        END SUBROUTINE assign_section


        SUBROUTINE assign_var(left, right)
            TYPE(variable), INTENT(OUT), ALLOCATABLE :: left(:)
            TYPE(variable), INTENT(IN)    :: right(:)
            INTEGER :: i
            
            IF(ALLOCATED(left)) DEALLOCATE(left)
            ALLOCATE(left(SIZE(right)))

            DO i=1,SIZE(right)
                left(i)%var = right(i)%var
                left(i)%value = right(i)%value
            END DO
        END SUBROUTINE assign_var


        SUBROUTINE print_section(self)
            CLASS(section), INTENT(INOUT) :: self
            INTEGER :: i
            PRINT*, self%sect
            DO i=1,SIZE(self%vars)
                PRINT*, '    ',self%vars(i)%var,' = ',self%vars(i)%value
            END DO
        END SUBROUTINE print_section


        SUBROUTINE print_input_data(self)
            CLASS(input_data_type), INTENT(INOUT) :: self
            INTEGER :: i

            DO i=1,SIZE(self%sections)
                CALL self%sections(i)%print_section()
                PRINT*
            END DO
        END SUBROUTINE print_input_data


        SUBROUTINE read_input(self, file)
            CLASS(input_data_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: file ! file name
            INTEGER :: funit, io_status, ind, ind1
            CHARACTER(LEN=100) :: text
            CHARACTER(LEN=:), ALLOCATABLE :: sect_name, text1, text2
            
            ! opening file in reading mode
            OPEN(NEWUNIT=funit, FILE=file, ACTION='read', &
                IOSTAT=io_status, STATUS='old', POSITION='REWIND')
            IF(io_status/=0) THEN
                WRITE(error_unit,*) 'Unable to open file: '//file
                STOP 'stopping the execution' ! stops the code execution
            END IF

            DO WHILE(.TRUE.)
                READ(funit,'(A)',IOSTAT=io_status) text
                IF (io_status/=0) EXIT ! End of file
                text = TRIM(ADJUSTL(text))

                ! Discarding commment
                IF( text(1:1)==';' .OR. text(1:1)=='#' &
                & .OR. text(1:1)=='!' &
                & .OR. ( text(1:1)/='[' .AND. .NOT.(index(text,'=') > 0) ) ) THEN 
                    CYCLE
                END IF

                ! Sections and Variables
                IF(text(1:1)=='[') THEN ! section
                    ind = INDEX(text,']')
                    text = text(2:ind-1)
                    sect_name = TRIM(ADJUSTL(text))
                    !PRINT*, sect_name
                    CALL self%add_section(sect_name)
                    
                ELSE IF(INDEX(text,'=')>0) THEN ! variable
                    ! Splitting text
                    text1 = TRIM(ADJUSTL( tabs_to_spaces( text(1:INDEX(text,'=')-1) ) )) ! TRIM and ADJUSTL only remove spaces, not tabs.
                    text2 = TRIM(ADJUSTL( tabs_to_spaces( text(INDEX(text,'=')+1:) ) ))
                    IF (LEN_TRIM(text1)==0 .OR. LEN_TRIM(text2)==0) THEN
                        WRITE(error_unit,*) TRIM(text),' Empty variable name or value';
                        CYCLE
                    END IF
                    IF(INDEX(text2, ';')>0 .OR. INDEX(text2,'#')>0 .OR. INDEX(text2,'!')>0) THEN
                        ind1 = INDEX(text2, ';')
                        IF(ind1/=0) ind=ind1
                        ind1 = INDEX(text2, '#')
                        IF(ind1/=0 .AND. ind1<ind) ind=ind1
                        ind1 = INDEX(text2, '!')
                        IF(ind1/=0 .AND. ind1<ind) ind=ind1
                       
                        text2 = text2(1:ind-1)
                    END IF

                    CALL self%add_variable(sect_name, text1, text2)
                END IF

                IF(ALLOCATED(text1)) DEALLOCATE(text1)
                IF(ALLOCATED(text2)) DEALLOCATE(text2)
            END DO
            DEALLOCATE(sect_name)

            !CALL self%print_input_data()
        END SUBROUTINE read_input


        FUNCTION tabs_to_spaces(str) RESULT(str2)
            ! Replaces tabs in a string with spaces.
            CHARACTER(LEN=*), INTENT(IN):: str
            CHARACTER(LEN=:), ALLOCATABLE :: str2
            INTEGER :: i

            str2 = str
            DO i=1,LEN(str)
                IF(str2(i:i)==CHAR(9)) THEN ! CHAR(9) -> tab
                    str2(i:i)=' '
                END IF
            END DO
        END FUNCTION tabs_to_spaces


        FUNCTION get_index(self, sec_name, var_name) RESULT(ind)
            CLASS(input_data_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: sec_name, var_name
            INTEGER :: ind(2)
            INTEGER :: i

            DO i=1,SIZE(self%sections)
                IF(self%sections(i)%sect==sec_name) THEN
                    ind(1) = i
                    EXIT
                END IF
            END DO

            DO i=1,SIZE(self%sections(ind(1))%vars)
                IF(self%sections(ind(1))%vars(i)%var==var_name) THEN
                    ind(2) = i
                    EXIT
                END IF
            END DO
        END FUNCTION get_index


        SUBROUTINE assign_variable_real(self, sect_name, var_name, var)
            CLASS(input_data_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: sect_name, var_name
            REAL(KIND=rkind), INTENT(OUT) :: var
            INTEGER :: index(2), sec_i, var_i

            index = self%get_index(sect_name, var_name); 
            sec_i = index(1); var_i = index(2);
            READ(self%sections(sec_i)%vars(var_i)%value,*) var
        END SUBROUTINE assign_variable_real


        SUBROUTINE assign_variable_int(self, sect_name, var_name, var)
            CLASS(input_data_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: sect_name, var_name
            INTEGER, INTENT(OUT) :: var
            INTEGER :: index(2), sec_i, var_i

            index = self%get_index(sect_name, var_name); 
            sec_i = index(1); var_i = index(2);
            READ(self%sections(sec_i)%vars(var_i)%value,*) var
        END SUBROUTINE assign_variable_int


        SUBROUTINE assign_variable_char(self, sect_name, var_name, var)
            CLASS(input_data_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), INTENT(IN) :: sect_name, var_name
            CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: var
            INTEGER :: index(2), sec_i, var_i

            index = self%get_index(sect_name, var_name); 
            sec_i = index(1); var_i = index(2);
            var = self%sections(sec_i)%vars(var_i)%value
        END SUBROUTINE assign_variable_char

END MODULE mod_read_file
