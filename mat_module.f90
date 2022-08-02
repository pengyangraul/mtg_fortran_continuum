MODULE mat_module
IMPLICIT NONE
INTEGER, PARAMETER :: sp = selected_real_kind(6, 37)
INTEGER, PARAMETER :: dp = selected_real_kind(15, 307) 
INTEGER, PARAMETER :: qp = selected_real_kind(33, 4931) 
! Define Pauli matrices
INTEGER, PARAMETER :: np = dp
REAL(kind=np), DIMENSION (2,2) :: sigma_x = RESHAPE([0.,1.,1.,0.], [2,2])
COMPLEX(kind=np), DIMENSION (2,2) :: sigma_y = RESHAPE([(0.,0.), (0.,1.), (0.,-1.), (0.,0.)], [2,2])
REAL(kind=np), DIMENSION (2,2) :: sigma_z = RESHAPE([1.,0.,0.,-1.], [2,2])
REAL(kind=np), DIMENSION (2,2) :: sigma_0 = RESHAPE([1.,0.,0.,1.], [2,2])

public :: kron
interface kron
  procedure kron_real, kron_complex, kron_int
end interface kron
SAVE
CONTAINS
FUNCTION linspace(ni,nf,N)  ! Linearly spaced  points
    IMPLICIT NONE
    REAL(kind = np), INTENT(IN) :: ni, nf    ! initial point and final point
    INTEGER, INTENT(IN) :: N ! N points between the above two points (included)
    REAL(kind = np), DIMENSION(N) :: linspace
    INTEGER :: i
    if( N == 1 ) then
      linspace(1) = ni
    else
    FORALL (i = 1:N) linspace(i) = (i - 1)*1.0
    linspace = ni + (nf - ni)/(N-1)*linspace
    end if
END FUNCTION linspace

subroutine cumsum(b,a)
  implicit none
  real(kind = np), dimension(:), allocatable, intent(in) :: a
  real(kind = np), dimension(:), allocatable, intent(out) :: b
  integer :: ii, Ndim
  Ndim = size(a)
  allocate(b(Ndim))
  b(1) = a(1)
  do ii = 2, Ndim
  b(ii) = b(ii-1) + a(ii)
  end do
end subroutine

FUNCTION eye(N) ! idenitty of N x N
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N ! size of matrix
    REAL(kind=np), DIMENSION(N,N) :: eye
    INTEGER :: i,j ! loop index
    FORALL(i=1:N,j=1:N) eye(i,j) = (i/j)*(j/i)
END FUNCTION eye

FUNCTION diag_vec_mat(vec, ith)  ! create a diagonal matrix from a vector 
    IMPLICIT NONE
    REAL(kind=np), DIMENSION(:), INTENT(IN) :: vec
    INTEGER, INTENT(IN), OPTIONAL :: ith
    REAL(kind=np), DIMENSION(:,:), ALLOCATABLE :: diag_vec_mat
    INTEGER :: nth_diag, mat_dim, i
    IF (PRESENT(ith)) THEN
    nth_diag = ith
    ELSE
    nth_diag = 0
    END IF
    mat_dim = SIZE(vec) + ABS(nth_diag)

    IF (ALLOCATED(diag_vec_mat)) DEALLOCATE(diag_vec_mat) !Discard any lingering storage.
    ALLOCATE(diag_vec_mat(mat_dim,mat_dim))

    diag_vec_mat(:,:) = 0
    IF (nth_diag>=0) THEN
        FORALL (i = 1:SIZE(vec)) diag_vec_mat(i,i+nth_diag) = vec(i)
    ELSE
        FORALL (i = 1:SIZE(vec)) diag_vec_mat(i-nth_diag,i) = vec(i)
    END IF
END FUNCTION diag_vec_mat

FUNCTION diag_mat_vec(mat,ith) ! get the ith diagonal elements of a matrix
    IMPLICIT NONE
    REAL(kind = np), DIMENSION(:,:), INTENT(IN) :: mat
    INTEGER, OPTIONAL, INTENT(IN) :: ith
    REAL, DIMENSION(:), ALLOCATABLE :: diag_mat_vec
    INTEGER :: nth_diag, vec_length, i
    IF (PRESENT(ith)) THEN
    nth_diag = ith
    ELSE
    nth_diag = 0
    END IF

    IF (SIZE(mat,1) == SIZE(mat,2)) THEN
        vec_length = SIZE(mat,1) - ABS(nth_diag)
    ELSE
        WRITE(*,*) 'Non-square matrix was inputted!'
        STOP
    END IF

    IF (ALLOCATED(diag_mat_vec)) DEALLOCATE(diag_mat_vec) !Discard any lingering storage.
    ALLOCATE(diag_mat_vec(vec_length))

    IF (nth_diag >= 0) THEN
        FORALL (i = 1:vec_length) diag_mat_vec(i) = mat(i,i+nth_diag)
    ELSE
        FORALL (i = 1:vec_length) diag_mat_vec(i) = mat(i-nth_diag,i)
    END IF
  END FUNCTION diag_mat_vec

  pure function kron_real(A,B)
    !kron = Kronecker product of A and B, both two-dimensional arrays.
    !Considers the arrays to be addressed as A(row,column), despite any storage order arrangements. 
    IMPLICIT NONE
    REAL(kind=np), DIMENSION(:,:), INTENT(IN) :: A,B !Two-dimensional arrays, lower bound one.
    REAL(kind=np), DIMENSION(:,:), ALLOCATABLE :: kron_real !To be created to fit.
    INTEGER R,RA,RB,C,CA,CB,I,J ! assistant integers
    RA = SIZE(A,DIM = 1)  !Ascertain the upper bounds of the incoming arrays.
    CA = SIZE(A,DIM = 2)  !Their lower bounds will be deemed one,
    RB = SIZE(B,DIM = 1)  !And the upper bound as reported will correspond.
    CB = SIZE(B,DIM = 2)  !SIZE(A) would give an array of two values, RA and CA, more for higher dimensionality.

    !WRITE (*,100) "A",RA,CA,"B",RB,CB,"A.k.B",RA*RB,CA*CB  !Announce.
    !100 FORMAT (3(A," is ",I0,"x",I0,1X))  !Three sets of sizes.
    IF (ALLOCATED(kron_real)) DEALLOCATE(kron_real) !Discard any lingering storage.
    ALLOCATE (kron_real(RA*RB,CA*CB))    !Obtain the exact desired size.
    R = 0   !Syncopation: start the row offset.
    DO I = 1,RA !Step down the rows of A.
      C = 0   !For each row, start the column offset.
      DO J = 1,CA   !Step along the columns of A.
        kron_real(R + 1:R + RB,C + 1:C + CB) = A(I,J)*B  !Place a block of B values.
        C = C + CB    !Advance a block of columns.
      END DO    !On to the next column of A.
      R = R + RB    !Advance a block of rows.
    END DO  !On to the next row of A.
  end function  kron_real

  pure function kron_complex(A,B)
    !kron = Kronecker product of A and B, both two-dimensional arrays.
    !Considers the arrays to be addressed as A(row,column), despite any storage order arrangements. 
    IMPLICIT NONE
    complex(kind=np), DIMENSION(:,:), INTENT(IN) :: A,B !Two-dimensional arrays, lower bound one.
    complex(kind=np), DIMENSION(:,:), ALLOCATABLE :: kron_complex !To be created to fit.
    INTEGER R,RA,RB,C,CA,CB,I,J ! assistant integers
    RA = SIZE(A,DIM = 1)  !Ascertain the upper bounds of the incoming arrays.
    CA = SIZE(A,DIM = 2)  !Their lower bounds will be deemed one,
    RB = SIZE(B,DIM = 1)  !And the upper bound as reported will correspond.
    CB = SIZE(B,DIM = 2)  !SIZE(A) would give an array of two values, RA and CA, more for higher dimensionality.

    !WRITE (*,100) "A",RA,CA,"B",RB,CB,"A.k.B",RA*RB,CA*CB  !Announce.
    !100 FORMAT (3(A," is ",I0,"x",I0,1X))  !Three sets of sizes.
    IF (ALLOCATED(kron_complex)) DEALLOCATE(kron_complex) !Discard any lingering storage.
    ALLOCATE (kron_complex(RA*RB,CA*CB))    !Obtain the exact desired size.
    R = 0   !Syncopation: start the row offset.
    DO I = 1,RA !Step down the rows of A.
      C = 0   !For each row, start the column offset.
      DO J = 1,CA   !Step along the columns of A.
        kron_complex(R + 1:R + RB,C + 1:C + CB) = A(I,J)*B  !Place a block of B values.
        C = C + CB    !Advance a block of columns.
      END DO    !On to the next column of A.
      R = R + RB    !Advance a block of rows.
    END DO  !On to the next row of A.
  end function  kron_complex

  pure function kron_int(A,B)
    !kron = Kronecker product of A and B, both two-dimensional arrays.
    !Considers the arrays to be addressed as A(row,column), despite any storage order arrangements. 
    IMPLICIT NONE
    integer, DIMENSION(:,:), INTENT(IN) :: A,B !Two-dimensional arrays, lower bound one.
    integer, DIMENSION(:,:), ALLOCATABLE :: kron_int !To be created to fit.
    INTEGER R,RA,RB,C,CA,CB,I,J ! assistant integers
    RA = SIZE(A,DIM = 1)  !Ascertain the upper bounds of the incoming arrays.
    CA = SIZE(A,DIM = 2)  !Their lower bounds will be deemed one,
    RB = SIZE(B,DIM = 1)  !And the upper bound as reported will correspond.
    CB = SIZE(B,DIM = 2)  !SIZE(A) would give an array of two values, RA and CA, more for higher dimensionality.

    !WRITE (*,100) "A",RA,CA,"B",RB,CB,"A.k.B",RA*RB,CA*CB  !Announce.
    !100 FORMAT (3(A," is ",I0,"x",I0,1X))  !Three sets of sizes.
    IF (ALLOCATED(kron_int)) DEALLOCATE(kron_int) !Discard any lingering storage.
    ALLOCATE (kron_int(RA*RB,CA*CB))    !Obtain the exact desired size.
    R = 0   !Syncopation: start the row offset.
    DO I = 1,RA !Step down the rows of A.
      C = 0   !For each row, start the column offset.
      DO J = 1,CA   !Step along the columns of A.
        kron_int(R + 1:R + RB,C + 1:C + CB) = A(I,J)*B  !Place a block of B values.
        C = C + CB    !Advance a block of columns.
      END DO    !On to the next column of A.
      R = R + RB    !Advance a block of rows.
    END DO  !On to the next row of A.
  end function  kron_int

FUNCTION dirsum(A,B)  ! direct sum of A and B
    IMPLICIT NONE
    REAL(kind=np), DIMENSION (:,:), INTENT(IN)  :: A, B
    REAL(kind=np), DIMENSION (:,:), ALLOCATABLE :: dirsum
    INTEGER :: p,  q, RA, RB, CA, CB
    RA = SIZE(A,DIM=1)
    RB = SIZE(B,DIM=1)
    CA = SIZE(A,DIM=2)
    CB = SIZE(B,DIM=2)
    p = RA + RB
    q = CA + CB

    IF (ALLOCATED(dirsum)) DEALLOCATE(dirsum) !Discard any lingering storage.
    ALLOCATE(dirsum(p,p))
    dirsum = 0
    dirsum(1:RA,1:CA) = A
    dirsum(RA+1:p,CA+1:q) = B
    return
END FUNCTION dirsum

FUNCTION histogram(x, rmin, rmax, Nbins)
!  given a 1D array x, count the frequency of appearance within intervals [ranges(j),ranges(j+1)),
! ranges = linspace(rmin,rmax, Nbins +1 )
IMPLICIT NONE
REAL(kind = np), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: x
REAL(kind = np), INTENT(IN) :: rmin, rmax
INTEGER, INTENT(IN) :: Nbins
INTEGER, DIMENSION(:), ALLOCATABLE :: histogram
INTEGER :: jj, kk, info
REAL(kind = np) :: delta_bin
INTEGER, DIMENSION(:), ALLOCATABLE :: y
REAL(kind = np), DIMENSION(:), ALLOCATABLE :: xtemp
xtemp = pack(x,(x >= rmin))
xtemp = pack(xtemp,(xtemp < rmax)) - rmin
delta_bin = (rmax - rmin)/Nbins
y = floor(xtemp/delta_bin)
ALLOCATE(histogram(Nbins))
FORALL(jj=1:Nbins)
    histogram(jj) = COUNT(y == (jj-1))
END FORALL
END FUNCTION histogram

function gcd_bin(u_aux, v_aux)
  integer :: gcd_bin
  integer, intent(in) :: u_aux, v_aux
  integer :: k, t, u, v 
 
  u = abs(u_aux)
  v = abs(v_aux)
  if( u < v ) then
     t = u
     u = v
     v = t
  endif
  if( v == 0 ) then
     gcd_bin = u
     return
  endif
  k = 1
  do while( (mod(u, 2) == 0).and.(mod(v, 2) == 0) )
     u = u / 2
     v = v / 2
     k = k * 2
  enddo
  if( (mod(u, 2) == 0) ) then
     t = u
  else
     t = -v
  endif
  do while( t /= 0 )
     do while( (mod(t, 2) == 0) )
        t = t / 2
     enddo
     if( t > 0 ) then
        u = t
     else
        v = -t
     endif
     t = u - v
  enddo
  gcd_bin = u * k
end function gcd_bin

END MODULE mat_module
