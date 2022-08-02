PROGRAM main
USE f95_precision, ONLY: wp => dp
USE blas95
USE lapack95
USE mat_module
USE indnth
USE io_module
IMPLICIT NONE
! data dictionary
! ======================================================================
real(kind = np), parameter :: w = 102._np, w0 = 60._np, vF = 0.873534e6_np, hbar = 1.0546e-34_np, &
  pi = 3.1415926535897932385,  J_to_meV = 6.242e21, a0 = 2.4595e-10, &
  kd = 4*pi/(3*a0), echarge = 1.6e-19, epsilon_0 = 11.15*8.854e-12, d_separation = 3.0e-10, threshold = 1e-3
real(kind = np), dimension(2), parameter :: QQ1=[-sqrt(3._np)/2._np,3._np/2._np], QQ2 = [sqrt(3._np)/2._np,3._np/2._np]
complex(kind = np), parameter :: i = (0._np, 1._np) 
integer, parameter :: Nk = 3, bzN = 40, nL=5, num_of_Nk = 2*Nk+1
integer :: hamiltonian_dim, num_of_states, CN_num_states
real(kind = np), dimension(:), allocatable ::  Ne_layer_from_CN
real(kind = np), dimension(:,:,:), allocatable :: evsq_list
real(kind = np) :: ktheta, alpha, energy_scale, energy_scale_int_scr, LM, &
scaled_Dfield, field_strength_0, ren_int_scr_scale, filling, mu
real(kind = np), dimension(:), allocatable :: fillinglist
real(kind = np), dimension(:), allocatable :: Nklist, Nkoneslist, Dfield_list_0, nnlist, mmlist, Qxlist, Qylist, Dfield_list
real(kind = np), dimension(:), allocatable :: odd_layer_list, even_layer_list 
complex(kind = np), dimension(:,:), allocatable :: interlayer_ham
real(kind = np), dimension(bzN**2,2) :: klist
integer :: ii, jj, Ne, nfilling
CHARACTER(len=80) :: filename
! ======================================================================
! input parameters
real(kind = np) :: theta, eta, dfield
! ======================================================================
! Usage: ./main theta dfield
! Input from command lines: theta = read_param_real(1) eta = w0/w
theta = read_param_real(1)
dfield = read_param_real(2)
eta = w0/w
write(*,*) '================================================================'
write(*,*) '================PRINT PARAMETERS================================'
write(*,*) "Theta = ", theta, "   |   ", "eta = ", eta, "  |  ", "dfield = ", dfield
write(*,*) '================================================================'
! ======================================================================
theta = theta/180._np * pi
ktheta = sin(theta/2._np)*2._np*kD
energy_scale = hbar * J_to_meV * vF * ktheta
alpha = w/energy_scale
hamiltonian_dim = 2*nL*num_of_Nk**2
num_of_states = hamiltonian_dim * bzN**2
CN_num_states = num_of_states/2
LM = a0/(2*sin(theta/2))
energy_scale_int_scr =2._np* echarge**2*d_separation/(LM**2 *sqrt(3._np)*epsilon_0) *J_to_meV *4._np
ren_int_scr_scale = energy_scale_int_scr/energy_scale
! ======================================================================
! initialize parametersr
scaled_Dfield = dfield/energy_scale
Dfield_list_0 =linspace(-scaled_Dfield,scaled_Dfield,nL)/2._np
field_strength_0 = scaled_Dfield / (nL-1)
allocate(odd_layer_list(nL))
allocate(even_layer_list(nL))
odd_layer_list = 0.0_np
even_layer_list = 0.0_np
odd_layer_list(1:nL:2) = 1.0_np
even_layer_list(2:nL:2) = 1.0_np
! ======================================================================
allocate(interlayer_ham(hamiltonian_dim, hamiltonian_dim))
klist = bz_sampling(bzn, qq1, qq2)
interlayer_ham = tcoupling(alpha,eta)
nfilling = 16
fillinglist = linspace(0.5_np,6._np,nfilling)
open(unit=2,file="potentials.csv")
write(2,*) "filling ", "v1 ", "v2 ", "v3 ", "v4 ", "v5 "
do ii = 1,nfilling
  allocate(Dfield_list(nL))
  Dfield_list = Dfield_list_0
  filling = fillinglist(ii)
  call self_consistent_loop(Dfield_list,threshold, filling,mu)
  write(2,'(7e16.8)') filling, Dfield_list, mu
  deallocate(Dfield_list)
end do
close(2)
!
! ======================================================================
 CONTAINS
! INTERNAL PROCEDURES
  ! ======================================================================

! Define interlayer tunneling
function Tcoupling(alpha, eta)
  implicit none
  real(kind = np), intent(in) :: eta, alpha
  complex(kind = np), parameter :: xi = exp(i*pi*2/3._np), xi_star = exp(-i*pi*2/3._np)
  complex(kind = np), dimension(2,2) :: TT0, TT1, TT2
  real(kind = np), dimension(:,:), allocatable :: diagp, diag_Nk, delta_nn_mp1m, delta_np1n_mm, &
  delta_nn_mm, delta_nn_mm1m, delta_nm1n_mm
  complex(kind = np), dimension(:,:), allocatable :: Tcoupling, Toe, Teo

  TT0 = reshape([eta, 1._np, 1._np, eta]*(1._np,0._np),[2,2])
  TT1 = reshape([eta*(1._np,0._np), xi_star, xi, eta*(1._np,0._np)],[2,2])
  TT2 = conjg(TT1)

  allocate(diagp(num_of_Nk,num_of_Nk))
  allocate(diag_Nk(num_of_Nk,num_of_Nk))
  diagp = 0._np
  diag_Nk = 0._np
  do concurrent (jj=1:(num_of_Nk-1))
    diagp(jj,jj+1) = 1.0_np
    diag_Nk(jj,jj) = 1.0_np
  end do
  diag_Nk(num_of_Nk,num_of_Nk) = 1._np
  delta_nn_mp1m = kron(diag_Nk,diagp)
  delta_np1n_mm = kron(diagp,diag_Nk)
  delta_nn_mm = kron(diag_Nk,diag_Nk)
  allocate(delta_nn_mm1m(num_of_Nk**2,num_of_Nk**2))
  allocate(delta_nm1n_mm(num_of_Nk**2,num_of_Nk**2))
  delta_nn_mm1m = transpose(delta_nn_mp1m)
  delta_nm1n_mm = transpose(delta_np1n_mm)
  allocate(Toe(2*num_of_Nk**2,2*num_of_Nk**2))
  allocate(Teo(2*num_of_Nk**2,2*num_of_Nk**2))
  Toe = kron(delta_nn_mm*(1._np,0._np),TT0) + kron(delta_np1n_mm*(1._np,0._np),TT1) + kron(delta_nn_mp1m*(1._np,0._np), TT2)
  Teo = kron(delta_nn_mm*(1._np,0._np),TT0) + kron(delta_nm1n_mm*(1._np,0._np),TT1) + kron(delta_nn_mm1m*(1._np,0._np), TT2)
  allocate(Tcoupling(2*nL*num_of_Nk**2,2*nL*num_of_Nk**2))
  Tcoupling = kron(diag_vec_mat(odd_layer_list(1:(num_of_Nk-1)),1)*(1._np,0._np), Toe) + kron(diag_vec_mat(even_layer_list(1:(num_of_Nk-1)),1)*(1._np,0._np), Teo)
  Tcoupling = (Tcoupling + conjg(transpose(Tcoupling)))*alpha
end function

! Generating odd/even layer hamiltonian
subroutine gen_odd_even_ham(hodd,heven,k)
  implicit none
  real (kind = np), dimension(2), parameter :: q0 = [0._np,-1._np], q1=[-sqrt(3._np)/2._np,1._np/2._np], q2 = [sqrt(3._np)/2._np,1._np/2._np], qh = [sqrt(3._np)/2._np, 0._np]
  real(kind = np), dimension(:), allocatable :: Nklist, Nkoneslist, nnlist, mmlist, Qxlist, Qylist, kxlist_odd, kylist_odd, kxlist_even, kylist_even
  complex(kind = np), dimension(2,2) :: sigma_x_mtheta, sigma_y_mtheta, sigma_x_ptheta, sigma_y_ptheta, exp_itheta_z, exp_itheta_z_conj
  real(kind = np), dimension(2), intent(in) :: k
  complex(kind = np), dimension(:,:), intent(out), allocatable :: hodd, heven
  complex(kind = np), dimension(:,:), allocatable :: kxlist_odd_mat, kylist_odd_mat, kxlist_even_mat, kylist_even_mat
  integer :: num_of_Nk_2 = num_of_Nk**2

  exp_itheta_z = reshape([exp(i*theta), (0._np,0._np), (0._np,0._np), exp(-i*theta)],[2,2])
  exp_itheta_z_conj = conjg(exp_itheta_z)
  sigma_x_ptheta = matmul(exp_itheta_z,sigma_x)
  sigma_y_ptheta = matmul(exp_itheta_z,sigma_y)
  sigma_x_mtheta = matmul(exp_itheta_z_conj,sigma_x)
  sigma_y_mtheta = matmul(exp_itheta_z_conj,sigma_y)

  Nklist = linspace(-Nk*1._np, Nk*1._np, num_of_Nk)
  allocate(Nkoneslist(num_of_Nk))
  Nkoneslist = 1._np
  allocate(nnlist(num_of_Nk_2))
  allocate(mmlist(num_of_Nk_2))
  do concurrent (jj=1:num_of_Nk, ii = 1:num_of_Nk)
      nnlist(jj+(ii-1)*num_of_Nk) = Nklist(ii) * Nkoneslist(jj)
      mmlist(jj+(ii-1)*num_of_Nk) = Nklist(jj) * Nkoneslist(ii)
  end do
  allocate(Qxlist(num_of_Nk_2))
  allocate(Qylist(num_of_Nk_2))
  allocate(kxlist_odd(num_of_Nk_2))
  allocate(kxlist_even(num_of_Nk_2))
  allocate(kylist_odd(num_of_Nk_2))
  allocate(kylist_even(num_of_Nk_2))
  Qxlist = nnlist*QQ1(1) + mmlist*QQ2(1) + qh(1)
  Qylist = nnlist*QQ1(2) + mmlist*QQ2(2) + qh(2)
  kxlist_odd = Qxlist + k(1) - q0(1)/2
  kylist_odd = Qylist + k(2) - q0(2)/2
  kxlist_even = Qxlist + k(1) + q0(1)/2
  kylist_even = Qylist + k(2) + q0(2)/2
  allocate(hodd(2*num_of_Nk_2,2*num_of_Nk_2))
  allocate(heven(2*num_of_Nk_2,2*num_of_Nk_2))
  allocate(kxlist_odd_mat(num_of_Nk_2,num_of_Nk_2))
  allocate(kylist_odd_mat(num_of_Nk_2,num_of_Nk_2))
  allocate(kxlist_even_mat(num_of_Nk_2,num_of_Nk_2))
  allocate(kylist_even_mat(num_of_Nk_2,num_of_Nk_2))
  kxlist_odd_mat = cmplx(diag_vec_mat(kxlist_odd,0))
  kylist_odd_mat = cmplx(diag_vec_mat(kylist_odd,0))
  kxlist_even_mat = cmplx(diag_vec_mat(kxlist_even,0))
  kylist_even_mat = cmplx(diag_vec_mat(kylist_even,0))
  
  hodd = kron(kxlist_odd_mat, sigma_x_mtheta) + kron(kylist_odd_mat, sigma_y_mtheta)
  heven = kron(kxlist_even_mat, sigma_x_ptheta) + kron(kylist_even_mat, sigma_y_ptheta)
end subroutine

function gen_ham0(k, Dfield_list)
  implicit none
  real(kind = np), dimension(2), intent(in) :: k
  complex(kind = np), dimension(:,:), allocatable :: hodd, heven, gen_ham0
  real(kind = np), dimension(:), allocatable, intent(in) :: Dfield_list
  complex(kind = np), dimension(:,:), allocatable :: odd_layer_mat,   even_layer_mat
  call gen_odd_even_ham(hodd,heven,k)
  allocate(gen_ham0(hamiltonian_dim, hamiltonian_dim))
  allocate(odd_layer_mat(nL,nL))
  allocate(even_layer_mat(nL,nL))
  odd_layer_mat = cmplx(diag_vec_mat(odd_layer_list,0))
  even_layer_mat = cmplx(diag_vec_mat(even_layer_list,0))
  gen_ham0 = kron(odd_layer_mat,hodd)+ kron(even_layer_mat,heven) + kron(diag_vec_mat(Dfield_list,0), eye(2*num_of_Nk**2))
end function

function construct_ham(k,Dfield_list,interlayer_ham)
  implicit none
  real(kind = np), dimension(2), intent(in) :: k
  real(kind = np), dimension(:), allocatable, intent(in) :: Dfield_list
  complex(kind = np), dimension(:,:), allocatable, intent(in) :: interlayer_ham
  complex(kind = np), dimension(:,:), allocatable :: construct_ham
  allocate(construct_ham(hamiltonian_dim,hamiltonian_dim))
  construct_ham = gen_ham0(k,Dfield_list) + interlayer_ham
end function

subroutine spectrum(eigenenergies, hamiltonian)
implicit none
real(kind = np), dimension(:), allocatable, intent(out) :: eigenenergies
complex(kind = np), dimension(:,:), allocatable, intent(inout) :: hamiltonian 
allocate(eigenenergies(hamiltonian_dim))
call heevd(hamiltonian,eigenenergies, jobz='V') ! hamiltonian is overwritten by eigenvectors (columns)
end subroutine


! Sample points in BZ
FUNCTION bz_sampling(n,Q1,Q2)
! Sample k points in the BZ (k1,k2) = [-pi to pi, -pi to pi)/*(2*pi)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  real(kind = np), dimension(2), intent(in) :: Q1, Q2
  REAL(kind = np), DIMENSION(n**2,2) :: bz_sampling
  INTEGER :: ii,jj ! loop indices
  REAL(kind = np), DIMENSION(n+1) :: klist

  klist = linspace(-0.5_np,0.5_np,n+1)
  do concurrent (ii = 1:n, jj = 1:n)
    bz_sampling((ii-1)*n+jj,1) =Q1(1)*klist(ii) + Q2(1)*klist(jj)
    bz_sampling((ii-1)*n+jj,2) =Q1(2)*klist(ii) + Q2(2)*klist(jj)
  end do
END FUNCTION bz_sampling

subroutine  EE_EV_list(eigenenergy_list,evsq_list,Dfield_list,interlayer_ham)
  implicit none
  real(kind = np), dimension(:,:,:), allocatable, intent(out) :: evsq_list
  real(kind = np), dimension(:,:), allocatable, intent(out) :: eigenenergy_list
  real(kind = np), dimension(:), allocatable, intent(in) :: Dfield_list
  real(kind = np), dimension(:), allocatable  :: eigenenergies
  complex(kind = np), dimension(:,:), allocatable :: full_ham
  complex(kind = np), dimension(:,:), allocatable, intent(in) :: interlayer_ham
  integer :: ii

  allocate(evsq_list(hamiltonian_dim,hamiltonian_dim,bzN**2))
  allocate(eigenenergy_list(hamiltonian_dim,bzN**2))
  do ii = 1, bzN**2
    full_ham = construct_ham(klist(ii,:),Dfield_list,interlayer_ham)
    call spectrum(eigenenergies, full_ham)
    eigenenergy_list(:,ii) = eigenenergies
    evsq_list(:,:,ii) = abs(full_ham)**2
  end do
end subroutine

subroutine layerNe(Nelist,mu,ee,evsq,Ne)
  implicit none 
  real(kind = np), dimension(:), allocatable, intent(out) :: Nelist
  real(kind = np), dimension(:,:), allocatable, intent(in) :: ee
  real(kind = np), dimension(:,:,:), allocatable, intent(in) :: evsq
  real(kind = np), dimension(:,:,:), allocatable :: tmparray
  integer, intent(inout) :: Ne
  real(kind = np), intent(out) :: mu
  REAL(kind = np), DIMENSION(:), ALLOCATABLE :: eeflat  
  REAL(kind = np), DIMENSION(:,:), ALLOCATABLE :: evsqflat, evsqflatnew
  LOGICAL, DIMENSION(:), ALLOCATABLE:: newidx
  integer :: Nlargest, jj
  eeflat = pack(ee,.True.)
  Nlargest = d_indnth(eeflat,Ne)
  mu = eeflat(Nlargest)
  newidx = PACK(ee<=mu+1e-10,.TRUE.)
  Ne = count(newidx)
  evsqflat = reshape(evsq,[hamiltonian_dim,hamiltonian_dim*bzN**2])
  allocate(evsqflatnew(hamiltonian_dim, Ne))
  !do concurrent(jj = 1:490)  
  forall(jj = 1:hamiltonian_dim)  
  !do jj = 1,hamiltonian_dim
  !  write(*,*) jj
    evsqflatnew(jj,:) = pack(evsqflat(jj,:),newidx)
  end forall 
  !end do
  deallocate(evsqflat)
  !allocate(evsq(2*num_of_Nk**2,nL,Ne))
  tmparray = reshape(evsqflatnew,[2*num_of_Nk**2,nL,Ne])
  !write(*,*) shape(tmparray)
  deallocate(evsqflatnew)
  evsqflat = sum(tmparray,1)
  deallocate(tmparray)
  Nelist = sum(evsqflat,2)
end subroutine

subroutine update_Ne_layer(Ne_layer_from_CN, mu, Dfield_list, interlayer_ham, Ne)
  implicit none
  integer, intent(inout) :: Ne
  real(kind = np), dimension(:), allocatable, intent(in) :: Dfield_list
  real(kind = np), dimension(:), allocatable, intent(out) :: Ne_layer_from_CN
  real(kind = np), intent(out) :: mu
  real(kind = np), dimension(:,:), allocatable :: eigenenergy_list
  real(kind = np), dimension(:,:,:), allocatable :: evsq_list
  complex(kind = np), dimension(:,:), allocatable, intent(in) :: interlayer_ham
  real(kind = np), dimension(:), allocatable :: Nelist
  call EE_EV_list(eigenenergy_list,evsq_list,Dfield_list,interlayer_ham)
  call layerNe(Nelist,mu,eigenenergy_list,evsq_list,Ne)
  allocate(Ne_layer_from_CN(nL))
  Ne_layer_from_CN = Nelist - CN_num_states*1._np/nL
end subroutine

function field_list(Ne_layer_from_CN)
  implicit none
  real(kind = np), dimension(:), allocatable, intent(in) :: Ne_layer_from_CN
  real(kind = np), dimension(:), allocatable :: field_pos_cumsum, field_neg_cumsum, field_list
  integer :: Ndim, jj, ii
  Ndim = size(Ne_layer_from_CN)
  allocate(field_pos_cumsum(Ndim-1))
  allocate(field_neg_cumsum(Ndim-1))
  field_pos_cumsum(1) = Ne_layer_from_CN(1)
  field_neg_cumsum(Ndim-1) = Ne_layer_from_CN(Ndim)
  do jj=2,(Ndim-1)
  ii = Ndim - jj
  field_pos_cumsum(jj) = field_pos_cumsum(jj-1) + Ne_layer_from_CN(jj)
  field_neg_cumsum(ii) = field_neg_cumsum(ii+1) + Ne_layer_from_CN(ii+1)
  end do
  field_list = (field_pos_cumsum - field_neg_cumsum)/2._np
end function

function update_Dfields(Ne_layer_from_CN)
  implicit none
  real(kind = np), dimension(:), allocatable, intent(in) :: Ne_layer_from_CN
  real(kind = np), dimension(:), allocatable :: F_list, F_list_cumsum
  real(kind = np), dimension(:), allocatable :: update_Dfields

  F_list = field_list(Ne_layer_from_CN)/(bzN**2) * ren_int_scr_scale + field_strength_0
  call cumsum(F_list_cumsum,F_list)
  allocate(update_Dfields(nL))
  update_Dfields = F_list_cumsum(nL-1)/2._np
  update_Dfields(2:nL) =  update_Dfields(2:nL) - F_list_cumsum
end function

function NeCN_from_filling(filling)
  implicit none
  real(kind = np), intent(in) :: filling
  real(kind = np) :: filling_single
  integer :: NeCN_from_filling
  filling_single = filling / 4._np
  NeCN_from_filling = int(filling_single*bzN**2)
end function

subroutine self_consistent_loop(Dfield_list,threshold, filling, mu)
  implicit none
  real(kind = np), dimension(:), allocatable, intent(inout) :: Dfield_list
  real(kind = np), intent(in) :: threshold
  real(kind = np), intent(inout) :: filling
  real(kind = np), intent(out) :: mu
  integer :: Ne_from_CN, Ne, jj
  real(kind = np), dimension(:), allocatable :: Ne_layer_from_CN_old, Ne_layer_from_CN
  real(kind = np) :: error
  error = 10._np
  Ne_from_CN = NeCN_from_filling(filling)
  allocate(Ne_layer_from_CN_old(nL))
  allocate(Ne_layer_from_CN(nL))
  Ne = Ne_from_CN + CN_num_states
  Ne_layer_from_CN_old = Ne_from_CN*1._np/nL
  jj = 0
  do while (error > threshold)
      jj = jj + 1
      if (jj >= 100) then
        exit
      end if
      call update_Ne_layer(Ne_layer_from_CN, mu, Dfield_list, interlayer_ham, Ne)
      Dfield_list = update_Dfields(Ne_layer_from_CN)
      error = norm2(Ne_layer_from_CN - Ne_layer_from_CN_old)/(nL*bzN**2)
      write(*,*) "......................................................."
      write(*,*) "Loop ", jj
      write(*,*) "Error = ", error
      write(*,*) "Nu_layers: ", Ne_layer_from_CN/(bzN**2)
      write(*,*) "......................................................."
      Ne_layer_from_CN_old = Ne_layer_from_CN
  end do
  filling = (Ne - CN_num_states)*1._np/(bzN**2)*4._np
  write(*,*) "======================================="
  write(*,*) "     Filling = ", filling  
  write(*,*) "||     Self-consistency is met!      ||"
  write(*,*) "======================================="
end subroutine

END PROGRAM main


