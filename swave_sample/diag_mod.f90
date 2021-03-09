MODULE diag_mod
  USE nrtype
  IMPLICIT NONE



CONTAINS

  SUBROUTINE eigen(JOBZ,A,W)
    USE nrtype
    ! This is a wrapper for the lapack zheevd subroutine, designed to make it easier for use. It computers the eigenvalues and eigenvectors of a Hermitian matrix 
    ! A: Matrix, Eigenvectors
    ! W: eigenvalues
    IMPLICIT NONE
!    INTEGER,   INTENT(out),OPTIONAL  :: INFO
    CHARACTER(len=1), INTENT(in)     :: JOBZ ! calculate eigenvalues only = 'N', calculate eigenvalues and vectors = 'V'
    character(1) :: UPLO
    REAL(DP), DIMENSION(:),  INTENT(inout) :: W
    COMPLEX(DP), DIMENSION(:,:),  INTENT(inout) :: A

    INTEGER                  :: N, LDA, info_inner, lwork, lrwork, liwork
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rwork
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: work
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwork

    UPLO='U'
    
    N= SIZE(A,1)
    IF (SIZE(W).ne.N) STOP 'ZHEEVD95: dimensions of matrix and array for EV do not match!'
    lda=N
    IF (.NOT. ((jobz.eq.'N') .OR. (JOBZ.eq.'V'))) STOP 'ZHEEVD95: Invalid option, must be "N" or "V"!'

    lwork  =  2*N + N**2
    lrwork = 1 + 5*N + 2*N**2
    liwork = 3 + 5*N
    ALLOCATE(work(lwork),rwork(lrwork),iwork(liwork))
    CALL ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, -1, &
&      RWORK, -1, IWORK, -1, INFO_inner ) 
    lwork = work(1)
    liwork = iwork(1)
    lrwork = rwork(1)
    DEALLOCATE(work,rwork,iwork)
    ALLOCATE(work(lwork),rwork(lrwork),iwork(liwork))

    CALL ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, &
&      RWORK, LRWORK, IWORK, LIWORK, INFO_inner ) 
    DEALLOCATE(work, rwork, iwork)

!    IF (PRESENT(info)) info = info_inner
  END SUBROUTINE eigen


END MODULE diag_mod
