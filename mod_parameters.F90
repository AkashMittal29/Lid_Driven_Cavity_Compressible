
MODULE mod_parameters
    IMPLICIT NONE
    PUBLIC

#ifdef SINGLE_PRECISION
    ! SINGLE_PRCISION is a Preprocessor Directive
    INTEGER, PARAMETER :: rkind = SELECTED_REAL_KIND(P=6, R=30) ! single precision
#else
    INTEGER, PARAMETER :: rkind = SELECTED_REAL_KIND(P=15, R=300) ! double
#endif

    REAL(rkind), PARAMETER :: pi = acos(-1._rkind)

END MODULE mod_parameters

