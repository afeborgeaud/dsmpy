SUBROUTINE GLU( A, N, N1, B, EPS, WK, IP, IER )
!        GLU, GLUSUB
!               COPYRIGHT : H.HASEGAWA, AUG. 26 1989 V.1
!
!               SOLVES SIMULTANEOUS LINEAR EQUATIONS
!               BY GAUSSIAN ELIMINATION METHOD.
!
!        INPUT - -
!             A(N1,N)  R *8  : 2-DIM. ARRAY CONTAINING THE COEFFICIENTS.
!             N        I *4  : ORDER OF MATRIX.
!             N1       I *4  : SIZE OF ARRAY A.
!             B(N)     R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT-HAND
!                              SIDE VECTOR.
!             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE
!                              MATRIX. ( STANDARD VALUE 3.52D-15 )
!        OUTPUT - -
!             A(N1,N)        : RESULT OF GAUSSIAN ELIMINATION.
!             B(N)           : SOLUTION.
!             IP(N)    I *4  : PIVOT NUMBER.
!             IER      I *4  : = 0,  FOR NORMAL EXECUTION.
!                              = 1,  FOR SINGULAR MATRIX.
!                              = 2,  FOR SINGULAR ORIGINAL MATRIX.
!                              = 3,  FOR INVALID ARGUEMENT.
!        WORKING  -
!             WK(N)    R *8  : 1-DIM. ARRAY.
! modified : Kensuke Konishi 2018 in Paris
    use parameters
    implicit none
    INTEGER:: N,N1
    INTEGER:: IP(N),IER
    double precision:: EPS
    COMPLEX(dp)::A(N1,N),B(N),WK(N)
    INTEGER:: I,J,K,IPK
    COMPLEX(dp):: AMAX,AIK,W,T
    !             LEFT-HAND SIDE
    IF( EPS<0.0D0 )  EPS = 3.52D-15
    IF( ( N1<N ).OR.( N<=0 ) )  THEN
        IER = 3
        WRITE(*,*) '  (SUBR. GLU)  INVALID ARGUMENT.  N1, N =', N1, N
        RETURN
    END IF
    !             CHECK ORIGINAL MATRIX.
    WK(1:n) = cdabs(A(1:n,1))

    DO J = 2, N
        WK(1:n) = DMAX1( cdabs( WK(1:n) ), cdabs(A(1:n,J)) )
    enddo
    DO I = 1, N
        IF( cdabs( WK(I) )<EPS )  THEN
            IER = 2
            WRITE(*,*) '  (SUBR. GLU)  ORIGINAL MATRIX IS SINGULAR.'
            RETURN
        END IF
    enddo
    !
    IER = 0
    DO K = 1, N
        !             FIND MAXIMUM ELEMENT IN THE K-TH COLUMN.
        AMAX = cdabs(A(K,K))
        IPK = K
        DO I = K+1, N
            AIK = cdabs(A(I,K))
            IF( cdabs( AIK )>cdabs( AMAX ) )  THEN
                IPK = I
                AMAX = AIK
            END IF
        enddo
        IP(K) = IPK
        !
        IF( cdabs( AMAX )>EPS )  THEN
            IF( IPK/=K )  THEN
                W = A(IPK,K)
                A(IPK,K) = A(K,K)
                A(K,K) = W
            END IF
            !             COMPUTE ALFA
            DO I = K+1, N
                A(I,K) = -A(I,K)/A(K,K)
                WK(I) = A(I,K)
            enddo

            DO J = K+1, N
                IF( IPK/=K )  THEN
                    W = A(IPK,J)
                    A(IPK,J) = A(K,J)
                    A(K,J) = W
                END IF
                !             GAUSSIAN ELIMINATION
                T = A(K,J)

                A(K+1:N,J) = A(K+1:N,J) + WK(K+1:N)*T

            enddo
        !             MATRIX IS SINGULAR.
        ELSE
            IER = 1
            IP(K) = K

            A(k+1:n,K) = 0.0D0

            WRITE(*,*) '  (SUBR. GLU)  MATRIX IS SINGULAR AT K =', K
            RETURN
        END IF
    enddo
    !             RIGHT-HAND SIDE
    !	ENTRY GLUSUB( A, B )
    ENTRY GLUSUB( A, N, N1, B, EPS, WK, IP, IER )
    !             FORWARD ELIMINATION PROCESS
    DO K = 1, N
        IF( IP(K)/=K ) THEN
            W = B(IP(K))
            B(IP(K)) = B(K)
            B(K) = W
        END IF
        !
        T = B(K)
        B(K+1:N) = B(K+1:N) + A(K+1:N,K)*T

    enddo
    !             BACKWARD SUBSTITUTION PROCESS
    B(N) = B(N)/A(N,N)
    DO K = N-1, 1, -1
        T = B(K+1)
        B(1:k) = B(1:k) - A(1:k,K+1)*T
        B(K) = B(K)/A(K,K)
    enddo
    RETURN
END
