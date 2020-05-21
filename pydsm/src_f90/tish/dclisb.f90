SUBROUTINE dclisb0_pretreatment(A, N, NUD, N1, B, EPS, DR, Z, IER)
!************************************************************************
!*  SIMULTANEOUS LINEAR EQUATIONS WITH REAL SYMMETRIC POSITIVE DEFINITE *
!*      BAND MATRIX BY CHOLESKY METHOD.                                 *
!*  PARAMETERS                                                          *
!*    (1) A : 2-DIM. ARRAY CONTAINING THE MATRIX.                       *
!*    (2) N : ORDER OF THE MATRIX.                                      *
!*    (3) NUD : SIZE OF BAND'S HALF WIDTH.                              *
!*    (4) N1 : ROW SIZE OF THE ARRAY A IN THE 'DIMENSION' STATEMENT.    *
!*    (5) B : 1-DIM. ARRAY CONTAINING THE RIGHT HAND SIDE VECTOR.       *
!*    (6) EPS : PARAMETER TO CHECK SINGURARITY OFF THE MATRIX           *
!*              STANDARD VALUE = 1.0D-14                                *
!*    (7) DR : 1-DIM. WORKING ARRAY.                                    *
!*    (8) Z : 1-DIM. WORKING ARRAY.                                     *
!*    (9) IER : ERROR CODE.                                             *
!*  COPY RIGHT   T. OGUNI   JULY 30 1989   VERSION 1.0                  *
!* modified Kensuke Konishi 2018 in Paris
!************************************************************************
    use parameters
    implicit none
    COMPLEX(dp):: A(N1,N), B(N), DR(N), Z(N)
    double precision:: EPS
    INTEGER:: N, NUD, N1 ,IER
    COMPLEX(dp)::XX, S, SUM, AU, T
    double precision:: EPS1
    INTEGER:: I ,M, J, K1, MJ, K
    ! CHECK THE INPUT DATA
    IER = 0
    EPS1 = 1.0D-14
    M = NUD + 1
    IF ( N <= 0 .OR. NUD <= 0 .OR. N1 < M ) THEN
        IER = 2
        WRITE(*,*) '(SUBR. LISB) INVALID ARGUMENT. ', N, NUD, N1
        RETURN
    ENDIF
    IF (EPS <= 0.0) EPS = EPS1
    ! MODIFIED CHOLESKY DECOMPOSITION
    J = 1
    IF (CDABS(A(M,1)) <= EPS) THEN
        IER = 1
        WRITE(*,*) '(SUBR. LISB) SINGULAR AT STEP # ', J
        RETURN
    ENDIF
    DR(1) = DCMPLX(1.0D0) / A(M,1)
    XX = A(M-1,2)
    A(M-1,2) = A(M-1,2) * DR(1)
    S = A(M,2) - XX * A(M-1,2)
    J = 2
    IF (CDABS(S) <= EPS) THEN
        IER = 1
        WRITE(*,*) '(SUBR. LISB) SINGULAR AT STEP # ', J
        RETURN
    ENDIF
    DR(2) = 1 / S
    IF (M < 3) THEN
        DO J=3,N
            XX = A(1,J)
            A(1,J) = XX * DR(J-1)
            S = A(2,J) - XX * A(1,J)
            IF (CDABS(S) <= EPS) THEN
                IER = 1
                WRITE(*,*) ' (SUBR. LISB) SINGULAR AT STEP # ', J
                RETURN
            ENDIF
            DR(J) = 1 / S
        enddo
    ELSE
        DO J=3,N
            K1 = 1
            IF (J >= M) K1 = J - M + 1
            MJ = M - J
            DO I=K1+1,J-1
                SUM = 0
                DO K=K1,I-1
                    SUM = SUM + A(M-I+K,I) * A(MJ+K,J)
                enddo
                A(MJ+I,J) = A(MJ+I,J) - SUM
            enddo
            SUM = 0
            DO I=K1,J-1
                XX = A(MJ+I,J)
                AU = XX * DR(I)
                SUM = SUM + XX *AU
                A(MJ+I,J) = AU
            enddo
            T = A(M,J) - SUM
            IF (CDABS(T) <= EPS) THEN
                IER = 1
                WRITE(*,*) ' (SUBR. LISB) SINGULAR AT STEP # ', J
                RETURN
            ENDIF
            DR(J) = 1 / T
        enddo
    ENDIF
    RETURN
END

SUBROUTINE DCLISB0(A, N, NUD, N1, B, EPS, DR, Z, IER)
!************************************************************************
!*  SIMULTANEOUS LINEAR EQUATIONS WITH REAL SYMMETRIC POSITIVE DEFINITE *
!*      BAND MATRIX BY CHOLESKY METHOD.                                 *
!*  PARAMETERS                                                          *
!*    (1) A : 2-DIM. ARRAY CONTAINING THE MATRIX.                       *
!*    (2) N : ORDER OF THE MATRIX.                                      *
!*    (3) NUD : SIZE OF BAND'S HALF WIDTH.                              *
!*    (4) N1 : ROW SIZE OF THE ARRAY A IN THE 'DIMENSION' STATEMENT.    *
!*    (5) B : 1-DIM. ARRAY CONTAINING THE RIGHT HAND SIDE VECTOR.       *
!*    (6) EPS : PARAMETER TO CHECK SINGURARITY OFF THE MATRIX           *
!*              STANDARD VALUE = 1.0D-14                                *
!*    (7) DR : 1-DIM. WORKING ARRAY.                                    *
!*    (8) Z : 1-DIM. WORKING ARRAY.                                     *
!*    (9) IER : ERROR CODE.                                             *
!*  COPY RIGHT   T. OGUNI   JULY 30 1989   VERSION 1.0                  *
!* modified Kensuke Konishi 2018 in Paris
!************************************************************************
    use parameters
    implicit none
    COMPLEX(dp):: A(N1,N), B(N), DR(N), Z(N)
    double precision:: EPS
    INTEGER:: N, NUD, N1 ,IER

    call dclisb0_pretreatment(A, N, NUD, N1, B, EPS, DR, Z, IER)
    ENTRY DCSBSUB0(A, N, NUD, N1, B, EPS, DR, Z, IER)
    !  FORWARD SUBSTITUTION
    call dclisb0_kenja(A, N, NUD, N1, B, EPS, DR, Z, IER)
    RETURN
END

SUBROUTINE dclisb0_kenja(A, N, NUD, N1, B, EPS, DR, Z, IER)
!************************************************************************
!*  SIMULTANEOUS LINEAR EQUATIONS WITH REAL SYMMETRIC POSITIVE DEFINITE *
!*      BAND MATRIX BY CHOLESKY METHOD.                                 *
!*  PARAMETERS                                                          *
!*    (1) A : 2-DIM. ARRAY CONTAINING THE MATRIX.                       *
!*    (2) N : ORDER OF THE MATRIX.                                      *
!*    (3) NUD : SIZE OF BAND'S HALF WIDTH.                              *
!*    (4) N1 : ROW SIZE OF THE ARRAY A IN THE 'DIMENSION' STATEMENT.    *
!*    (5) B : 1-DIM. ARRAY CONTAINING THE RIGHT HAND SIDE VECTOR.       *
!*    (6) EPS : PARAMETER TO CHECK SINGURARITY OFF THE MATRIX           *
!*              STANDARD VALUE = 1.0D-14                                *
!*    (7) DR : 1-DIM. WORKING ARRAY.                                    *
!*    (8) Z : 1-DIM. WORKING ARRAY.                                     *
!*    (9) IER : ERROR CODE.                                             *
!*  COPY RIGHT   T. OGUNI   JULY 30 1989   VERSION 1.0                  *
!* modified Kensuke Konishi 2018 in Paris
!************************************************************************
    use parameters
    implicit none
    COMPLEX(dp):: A(N1,N), B(N), DR(N), Z(N)
    double precision:: EPS
    INTEGER:: N, NUD, N1 ,IER
    COMPLEX (dp)::SUM
    INTEGER:: M, J, I1, K, J1
    !  FORWARD SUBSTITUTION
    M = NUD + 1
    IF (M < 3) THEN
        Z(1) = B(1)
        DO J=2,N
            Z(J) = B(J) - A(1,J) * Z(J-1)
        enddo
        DO J=1,N
            Z(J) = Z(J) * DR(J)
        enddo
        B(N) = Z(N)
        DO J=1,N-1
            B(N-J) = Z(N-J) - A(1,N-J+1) * B(N-J+1)
        enddo
    ELSE
        Z(1) = B(1)
        Z(2) = B(2) - A(M-1,2) * Z(1)
        DO J=3,N
            IF (J > M) THEN
                I1 = 1
            ELSE
                I1 = M - J + 1
            ENDIF
            SUM = 0
            DO K=I1,M-1
                SUM = SUM + A(K,J) * Z(J-M+K)
            enddo
            Z(J) = B(J) - SUM
        enddo
        DO  J=1,N
            Z(J) = Z(J) * DR(J)
        enddo

        B(N) = Z(N)
        B(N-1) = Z(N-1) - A(M-1,N) * Z(N)
        DO J=3,N
            J1 = N - J + 1
            I1 = 1
            IF (J < M) I1 = M - J + 1
            SUM = 0
            DO K=I1,M-1
                SUM = SUM + A(K,M-K+J1) * B(M-K+J1)
            enddo
            B(J1) = Z(J1) - SUM
        enddo
    ENDIF
    RETURN
END


SUBROUTINE dclisb_pretreatment(A, N, NUD, N1, NP, B, EPS, DR, Z, IER)
!************************************************************************
!*  SIMULTANEOUS LINEAR EQUATIONS WITH REAL SYMMETRIC POSITIVE DEFINITE *
!*      BAND MATRIX BY CHOLESKY METHOD.                                 *
!*  PARAMETERS                                                          *
!*    (1) A : 2-DIM. ARRAY CONTAINING THE MATRIX.                       *
!*    (2) N : ORDER OF THE MATRIX.                                      *
!*    (3) NUD : SIZE OF BAND'S HALF WIDTH.                              *
!*    (4) N1 : ROW SIZE OF THE ARRAY A IN THE 'DIMENSION' STATEMENT.    *
!*    (5) B : 1-DIM. ARRAY CONTAINING THE RIGHT HAND SIDE VECTOR.       *
!*    (6) EPS : PARAMETER TO CHECK SINGURARITY OFF THE MATRIX           *
!*              STANDARD VALUE = 1.0D-14                                *
!*    (7) DR : 1-DIM. WORKING ARRAY.                                    *
!*    (8) Z : 1-DIM. WORKING ARRAY.                                     *
!*    (9) IER : ERROR CODE.                                             *
!*  COPY RIGHT   T. OGUNI   JULY 30 1989   VERSION 1.0                  *
!* modified Kensuke Konishi 2018 in Paris
!************************************************************************
    use parameters
    implicit none
    COMPLEX(dp):: A(N1,N), B(N), DR(N), Z(N)
    double precision:: EPS, EPS1
    INTEGER:: N, NUD, N1, NP, IER
    DOUBLE COMPLEX:: XX, S, SUM, AU, T
    INTEGER:: I ,M, J, K1, MJ, K
    !  CHECK THE INPUT DATA
    IER = 0
    EPS1 = 1.0D-14
    M = NUD + 1
    IF ( N <= 0 .OR. NUD <= 0 .OR. N1 < M)  THEN
        IER = 2
        WRITE(*,*) '(SUBR. LISB) INVALID ARGUMENT. ', N, NUD, N1
        RETURN
    ENDIF
    IF (EPS <= 0.0) EPS = EPS1
    ! MODIFIED CHOLESKY DECOMPOSITION
    J = 1
    IF (CDABS(A(M,1)) <= EPS) THEN
        IER = 1
        WRITE(*,*) '(SUBR. LISB) SINGULAR AT STEP # ', J
        RETURN
    ENDIF
    DR(1) = 1 / A(M,1)
    XX = A(M-1,2)
    A(M-1,2) = A(M-1,2) * DR(1)
    S = A(M,2) - XX * A(M-1,2)
    J = 2
    IF (CDABS(S) <= EPS) THEN
        IER = 1
        WRITE(*,*) '(SUBR. LISB) SINGULAR AT STEP # ', J
        RETURN
    ENDIF
    DR(2) = 1 / S
    IF (M < 3) THEN
        DO J=3,N
            XX = A(1,J)
            A(1,J) = XX * DR(J-1)
            S = A(2,J) - XX * A(1,J)
            IF (CDABS(S) <= EPS) THEN
                IER = 1
                WRITE(*,*) ' (SUBR. LISB) SINGULAR AT STEP # ', J
                RETURN
            ENDIF
            DR(J) = 1 / S
        enddo
    ELSE
        DO J=3,N
            K1 = 1
            IF (J >= M) K1 = J - M + 1
            MJ = M - J
            DO  I=K1+1,J-1
                SUM = 0
                DO K=K1,I-1
                    SUM = SUM + A(M-I+K,I) * A(MJ+K,J)
                enddo
                A(MJ+I,J) = A(MJ+I,J) - SUM
            enddo
            SUM = 0
            DO  I=K1,J-1
                XX = A(MJ+I,J)
                AU = XX * DR(I)
                SUM = SUM + XX *AU
                A(MJ+I,J) = AU
            enddo
            T = A(M,J) - SUM
            IF (CDABS(T) <= EPS) THEN
                IER = 1
                WRITE(*,*) ' (SUBR. LISB) SINGULAR AT STEP # ', J
                RETURN
            ENDIF
            DR(J) = 1 / T
        enddo
    ENDIF

    RETURN
END

SUBROUTINE DCLISB(A, N, NUD, N1, NP, B, EPS, DR, Z, IER)
!************************************************************************
!*  SIMULTANEOUS LINEAR EQUATIONS WITH REAL SYMMETRIC POSITIVE DEFINITE *
!*      BAND MATRIX BY CHOLESKY METHOD.                                 *
!*  PARAMETERS                                                          *
!*    (1) A : 2-DIM. ARRAY CONTAINING THE MATRIX.                       *
!*    (2) N : ORDER OF THE MATRIX.                                      *
!*    (3) NUD : SIZE OF BAND'S HALF WIDTH.                              *
!*    (4) N1 : ROW SIZE OF THE ARRAY A IN THE 'DIMENSION' STATEMENT.    *
!*    (5) B : 1-DIM. ARRAY CONTAINING THE RIGHT HAND SIDE VECTOR.       *
!*    (6) EPS : PARAMETER TO CHECK SINGURARITY OFF THE MATRIX           *
!*              STANDARD VALUE = 1.0D-14                                *
!*    (7) DR : 1-DIM. WORKING ARRAY.                                    *
!*    (8) Z : 1-DIM. WORKING ARRAY.                                     *
!*    (9) IER : ERROR CODE.                                             *
!*  COPY RIGHT   T. OGUNI   JULY 30 1989   VERSION 1.0                  *
!* modified Kensuke Konishi 2018 in Paris
!************************************************************************
    use parameters
    COMPLEX(dp):: A(N1,N), B(N), DR(N), Z(N)
    double precision:: EPS
    INTEGER:: N, NUD, N1, NP, IER

    call dclisb_pretreatment(A, N, NUD, N1, NP, B, EPS, DR, Z, IER)
    ! SUBTITUTION
    ENTRY DCSBSUB(A, N, NUD, N1, NP, B, EPS, DR, Z, IER)
    call dclisb_kenja(A, N, NUD, N1, NP, B, EPS, DR, Z, IER)

    RETURN
END

SUBROUTINE dclisb_kenja(A, N, NUD, N1, NP, B, EPS, DR, Z, IER)
!************************************************************************
!*  SIMULTANEOUS LINEAR EQUATIONS WITH REAL SYMMETRIC POSITIVE DEFINITE *
!*      BAND MATRIX BY CHOLESKY METHOD.                                 *
!*  PARAMETERS                                                          *
!*    (1) A : 2-DIM. ARRAY CONTAINING THE MATRIX.                       *
!*    (2) N : ORDER OF THE MATRIX.                                      *
!*    (3) NUD : SIZE OF BAND'S HALF WIDTH.                              *
!*    (4) N1 : ROW SIZE OF THE ARRAY A IN THE 'DIMENSION' STATEMENT.    *
!*    (5) B : 1-DIM. ARRAY CONTAINING THE RIGHT HAND SIDE VECTOR.       *
!*    (6) EPS : PARAMETER TO CHECK SINGURARITY OFF THE MATRIX           *
!*              STANDARD VALUE = 1.0D-14                                *
!*    (7) DR : 1-DIM. WORKING ARRAY.                                    *
!*    (8) Z : 1-DIM. WORKING ARRAY.                                     *
!*    (9) IER : ERROR CODE.                                             *
!*  COPY RIGHT   T. OGUNI   JULY 30 1989   VERSION 1.0                  *
!* modified Kensuke Konishi 2018 in Paris
!************************************************************************
    use parameters
    COMPLEX(dp):: A(N1,N), B(N), DR(N), Z(N)
    double precision:: EPS
    INTEGER:: N, NUD, N1, NP, IER
    DOUBLE COMPLEX:: SUM
    INTEGER:: M, J, I1, K
    ! FORWARD SUBSTITUTION
    M = NUD + 1
    IF (M < 3) THEN
        Z(NP) = B(NP)
        DO J=NP+1,N
            Z(J) = B(J) - A(1,J) * Z(J-1)
        enddo
        B(N) = Z(N) * DR(N)
    ELSE
        Z(NP) = B(NP)
        Z(NP+1) = B(NP+1) - A(M-1,NP+1) * Z(NP)
        DO J=NP+2,N
            IF (J > NP-1+M) THEN
                I1 = 1
            ELSE
                I1 = NP-1+M - J + 1
            ENDIF
            SUM = 0
            DO K=I1,M-1
                SUM = SUM + A(K,J) * Z(J-M+K)
            enddo
            Z(J) = B(J) - SUM
        enddo
        DO J=N-1,N
            Z(J) = Z(J) * DR(J)
        enddo

        B(N) = Z(N)
        B(N-1) = Z(N-1) - A(M-1,N) * Z(N)
    ENDIF

    RETURN
END


