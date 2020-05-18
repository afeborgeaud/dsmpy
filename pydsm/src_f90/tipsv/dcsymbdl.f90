SUBROUTINE DCSYMBDL(A,M,N,NN,EPS,Z,W,L,LI,LJ,IER)
!**********************************************************************
!  GAUSS METHOD FOR A SYMMETRIC BAND MATRIX WITH A SHORT WIDTH.       *
!  THIS ROUTINE IS FOR A VECTOR COMPUTER.                             *
!                                                                     *
!  PARAMETERS:                                                        *
!   ON ENTRY:                                                         *
!     A      THE ARRAY WITH DIMENSION (M+1)*N WHICH CONTAINS          *
!            THE LOWER BAND MATRIX IN ROW-WISE.                       *
!     M      THE HALF BAND WIDTH OF THE MATRIX A, EXCLUDING DIADONALS.*
!     N      THE ORDER OF THE MATRIX A.                               *
!     NN     THE ORDER OF WORKING ARRAYS L, LI, AND LJ.               *
!     EPS    THE TOLERANCE FOR PIVOTAL ELEMENTS.                      *
!   ON RETURN:                                                        *
!     A      THE LOWER TRIANGULAR MATRIX L AFTER DECOMPOSITION.       *
!     IER    ERROR CODE. IF IER=0, NORMAL RETURN.                     *
!   OTHERS: WORKING PARAMETERS.                                       *
!                                                                     *
!  COPYRIGHT:       FUMIKO NAGAHORI    1 SEP. 1991      VER. 1        *
!    modified: Kensuke Konishi 2018 in Paris
!**********************************************************************
    use parameters
    implicit none
    INTEGER:: N,M,MM,NN
    INTEGER:: L(NN),LI(NN),LJ(NN),IER
    double precision::EPS
    COMPLEX(dp):: A((M+1)*N),Z(M+1),W(M+1)
    INTEGER:: I,J,K,IJ,KK,NK,NKK,NKI
    COMPLEX(dp):: PIV
    !
    IER = 0
    IJ = 0
    DO  I=1,M
        IJ = IJ + 1
        LI(IJ) = I + 1
        LJ(IJ) = 2
        L(IJ) = M * (I-1) + 1
        DO  J=I+1,M
            IJ = IJ + 1
            L(IJ) = L(IJ-1) + M + 1
            LI(IJ) = LI(IJ-1) + 1
            LJ(IJ) = LJ(IJ-1) + 1
        enddo
    enddo
    MM = (M+1) * M / 2
    DO  K=1,N-M
        NK = (M+1) * (K-1) + 1
        IF (cdabs(A(NK+M)) < EPS) THEN
            WRITE(*,*) '(SUBR. SYMBDL) SINGULAR AT STEP = ', K
            IER = 1
            RETURN
        ENDIF
        PIV = 1.0D0 / A(NK+M)
        DO  J=2,M+1
            Z(J) = -A(NK+M*J)
            W(J) = A(NK+M*J) * PIV
            A(NK+M*J) = W(J)
        enddo

        KK = NK + M + M
        !VORTION VEC
        A(KK+L(1:mm)) = A(KK+L(1:mm)) + W(LJ(1:mm)) * Z(LI(1:mm))

    !
    enddo
    !
    DO K=N-M+1,N-1
        NK = (M+1) * (K-1) + 1
        NKK = (M+1) * K - 1
        IF (cdabs(A(NK+M)) < EPS) THEN
            WRITE(*,*) '(SUBR. SYMBDL) SINGULAR AT STEP = ', K
            IER = 1
            RETURN
        ENDIF
        PIV = 1.0D0 / A(NK+M)
        DO  J=2,N-K+1
            Z(J) = -A(NK+M*J)
            W(J) = A(NK+M*J) * PIV
            A(NK+M*J) = W(J)

            NKI = NKK + M * (J-1)

            A(NKI+2:nki+j) = A(NKI+2:nki+j) + W(J) * Z(2:j)

        enddo
    enddo
    !
    RETURN
!  END OF SYMBDL
END SUBROUTINE DCSYMBDL

SUBROUTINE DCSBDLV(A,B,M,N,NP,EPS,Z,IER)
!**********************************************************************
!  GAUSS METHOD FOR A SYMMETRIC BAND MATRIX WITH A SHORT WIDTH.       *
!  THIS ROUTINE IS FOR A VECTOR COMPUTER.                             *
!                                                                     *
!  PARAMETERS:                                                        *
!   ON ENTRY:                                                         *
!     A      THE ARRAY WITH DIMENSION (M+1),N WHICH CONTAINS          *
!            THE LEFT HAND SIDE LOWER BAND MATRIX.                    *
!     M      THE HALF BAND WIDTH OF THE MATRIX A, EXCLUDING DIADONALS.*
!     N      THE ORDER OF THE MATRIX A.                               *
!     EPS    THE TOLERANCE FOR PIVOTAL ELEMENTS.                      *
!     B      RIGHT HAND SIDE VECTOR                                   *
!   ON RETURN:                                                        *
!     B      SOLUTION.                                                *
!     IER    ERROR CODE. IF IER=0, NORMAL RETURN.                     *
!   OTHERS: WORKING PARAMETERS.                                       *
!                                                                     *
!  COPYRIGHT:       FUMIKO NAGAHORI    1 SEP. 1991      VER. 1        *
! modified: Kensuke Konishi 2018 in Paris
!**********************************************************************
    use parameters
    implicit none
    INTEGER:: M,N,NP,IER
    double precision::EPS
    COMPLEX(dp):: A(M+1,N),B(N),Z(N)
    INTEGER:: MM,J,K,I1
    COMPLEX(dp):: SUM
    !
    !  FORWARD SUBSTITUTION
    MM = M + 1
    IF (MM < 3) THEN
        Z(NP) = B(NP)
        DO J=NP+1,N
            Z(J) = B(J) - A(1,J) * Z(J-1)
        enddo
        B(N) = Z(N) / A(M+1,N)
    ELSE
        Z(NP) = B(NP)
        Z(NP+1) = B(NP+1) - A(MM-1,NP+1) * Z(NP)
        DO  J=NP+2,N
            IF (J > NP-1+MM) THEN
                I1 = 1
            ELSE
                I1 = NP-1+MM - J + 1
            ENDIF
            SUM = 0.0D0
            DO  K=I1,MM-1
                SUM = SUM + A(K,J) * Z(J-MM+K)
            enddo
            Z(J) = B(J) - SUM
        enddo

        Z(n-1:n) = Z(n-1:n) / A(M+1,n-1:n)

        !
        B(N) = Z(N)
        B(N-1) = Z(N-1) - A(MM-1,N) * Z(N)
    ENDIF
    !
    RETURN
END SUBROUTINE DCSBDLV

SUBROUTINE DCSYMBDL0(A,M,N,NN,EPS,Z,W,L,LI,LJ,IER)
!**********************************************************************
!  GAUSS METHOD FOR A SYMMETRIC BAND MATRIX WITH A SHORT WIDTH.       *
!  THIS ROUTINE IS FOR A VECTOR COMPUTER.                             *
!                                                                     *
!  PARAMETERS:                                                        *
!   ON ENTRY:                                                         *
!     A      THE ARRAY WITH DIMENSION (M+1)*N WHICH CONTAINS          *
!            THE LOWER BAND MATRIX IN ROW-WISE.                       *
!     M      THE HALF BAND WIDTH OF THE MATRIX A, EXCLUDING DIADONALS.*
!     N      THE ORDER OF THE MATRIX A.                               *
!     NN     THE ORDER OF WORKING ARRAYS L, LI, AND LJ.               *
!     EPS    THE TOLERANCE FOR PIVOTAL ELEMENTS.                      *
!   ON RETURN:                                                        *
!     A      THE LOWER TRIANGULAR MATRIX L AFTER DECOMPOSITION.       *
!     IER    ERROR CODE. IF IER=0, NORMAL RETURN.                     *
!   OTHERS: WORKING PARAMETERS.                                       *
!                                                                     *
!  COPYRIGHT:       FUMIKO NAGAHORI    1 SEP. 1991      VER. 1        *
!modified: Kensuke Konishi 2018 in Paris
!**********************************************************************
    use parameters
    implicit none
    integer:: N,M,MM,NN
    integer:: L(NN),LI(NN),LJ(NN),IER
    double precision:: EPS
    complex(dp):: A((M+1)*N),Z(M+1),W(M+1)
    integer:: I,J,K,IJ,KK,NK,NKK,NKI
    complex(dp):: PIV
    !
    IER = 0
    IJ = 0
    DO I=1,M
        IJ = IJ + 1
        LI(IJ) = I + 1
        LJ(IJ) = 2
        L(IJ) = M * (I-1) + 1
        DO  J=I+1,M
            IJ = IJ + 1
            L(IJ) = L(IJ-1) + M + 1
            LI(IJ) = LI(IJ-1) + 1
            LJ(IJ) = LJ(IJ-1) + 1
        enddo
    enddo
    MM = (M+1) * M / 2
    DO K=1,N-M
        NK = (M+1) * (K-1) + 1
        IF (cdabs(A(NK+M)) < EPS) THEN
            WRITE(*,*) '(SUBR. SYMBDL) SINGULAR AT STEP = ', K
            IER = 1
            RETURN
        ENDIF
        PIV = 1D0 / A(NK+M)
        DO J=2,M+1
            Z(J) = -A(NK+M*J)
            W(J) = A(NK+M*J) * PIV
            A(NK+M*J) = W(J)
        enddo
        !
        KK = NK + M + M
        !VORTION VEC
        DO  I=1,MM
            A(KK+L(I)) = A(KK+L(I)) + W(LJ(I)) * Z(LI(I))
        enddo
    enddo
    !
    !
    DO  K=N-M+1,N-1
        NK = (M+1) * (K-1) + 1
        NKK = (M+1) * K - 1
        IF (cdabs(A(NK+M)) < EPS) THEN
            WRITE(*,*) '(SUBR. SYMBDL) SINGULAR AT STEP = ', K
            IER = 1
            RETURN
        ENDIF
        PIV = 1D0 / A(NK+M)
        DO  J=2,N-K+1
            Z(J) = - A(NK+M*J)
            W(J) = A(NK+M*J) * PIV
            A(NK+M*J) = W(J)
            !
            NKI = NKK + M * (J-1)
            DO I=2,J
                A(NKI+I) = A(NKI+I) + W(J) * Z(I)
            enddo
        enddo
    enddo
    !
    RETURN
!  END OF SYMBDL
END

SUBROUTINE DCSBDLV0(A,B,M,N,EPS,Z,IER)
!**********************************************************************
!  GAUSS METHOD FOR A SYMMETRIC BAND MATRIX WITH A SHORT WIDTH.       *
!  THIS ROUTINE IS FOR A VECTOR COMPUTER.                             *
!                                                                     *
!  PARAMETERS:                                                        *
!   ON ENTRY:                                                         *
!     A      THE ARRAY WITH DIMENSION (M+1),N WHICH CONTAINS          *
!            THE LEFT HAND SIDE LOWER BAND MATRIX.                    *
!     M      THE HALF BAND WIDTH OF THE MATRIX A, EXCLUDING DIADONALS.*
!     N      THE ORDER OF THE MATRIX A.                               *
!     EPS    THE TOLERANCE FOR PIVOTAL ELEMENTS.                      *
!     B      RIGHT HAND SIDE VECTOR                                   *
!   ON RETURN:                                                        *
!     B      SOLUTION.                                                *
!     IER    ERROR CODE. IF IER=0, NORMAL RETURN.                     *
!   OTHERS: WORKING PARAMETERS.                                       *
!                                                                     *
!  COPYRIGHT:       FUMIKO NAGAHORI    1 SEP. 1991      VER. 1        *
!**********************************************************************
    use parameters
    IMPLICIT NONE
    integer:: M,N,IER
    double precision:: EPS
    complex(dp):: A(M+1,N),B(N),Z(N)
    integer:: MM,J,K,I1,J1
    complex(dp):: SUM
    !
    !  FORWARD SUBSTITUTION
    MM = M + 1
    IF (MM < 3) THEN
        Z(1) = B(1)
        DO J=2,N
            Z(J) = B(J) - A(1,J) * Z(J-1)
        enddo
        DO J=1,N
            Z(J) = Z(J) / A(M+1,J)
        enddo
        B(N) = Z(N)
        DO J=1,N-1
            B(N-J) = Z(N-J) - A(1,N-J+1) * B(N-J+1)
        enddo
    ELSE
        Z(1) = B(1)
        Z(2) = B(2) - A(M-1,2) * Z(1)
        DO J=3,N
            IF (J > MM) THEN
                I1 = 1
            ELSE
                I1 = MM - J + 1
            ENDIF
            SUM = 0
            DO K=I1,MM-1
                SUM = SUM + A(K,J) * Z(J-MM+K)
            enddo
            Z(J) = B(J) - SUM
        enddo

        Z(1:n) = Z(1:n) / A(M+1,1:n)

        B(N) = Z(N)
        B(N-1) = Z(N-1) - A(MM-1,N) * Z(N)
        DO J=3,N
            J1 = N - J + 1
            I1 = 1
            IF (J < MM) I1 = MM - J + 1
            SUM = DCMPLX(0.0D0)
            DO  K=I1,MM-1
                SUM = SUM + A(K,MM-K+J1) * B(MM-K+J1)
            enddo
            B(J1) = Z(J1) - SUM
        enddo
    ENDIF
    RETURN
END
