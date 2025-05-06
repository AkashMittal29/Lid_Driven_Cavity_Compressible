MODULE mod_restart
    USE mod_grid
    USE mod_field
    USE mod_parameters
    USE mod_utility
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : error_unit 
    IMPLICIT NONE

    !#TODO : write subroutine to binary open file, write file(full grid and field data from grid and field), close file.
    !#TODO : write subroutine to open binary file, read file(to field), close file. After reading, have to call primitive to conserved.
    !#TODO : write subroutine to write readable data from field into a file after reading restart file (make a call in main to do just this task when a flag is given in input file).



END MODULE mod_restart