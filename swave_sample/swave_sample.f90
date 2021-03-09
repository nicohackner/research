!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! swave_sample.f90
! ================
!
! This program solves the self-consistent gap equations for an impurity in
! an s-wave singlet superconductor on a finite lattice with periodic
! boundary conditions.
!
! PMRB, 3.3.2021
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM swave_sample
  USE nrtype
  USE diag_mod

  IMPLICIT NONE

  integer :: ix, iy, m1, m2
  integer :: counter
  integer, parameter :: Nlat=31 ! side length of square lattice
  integer :: iter
  integer, parameter :: itermax=500
  integer :: middle 
  complex(dp), dimension(Nlat,Nlat) :: Delta, Deltanew
  complex(dp), dimension(2*Nlat**2,2*Nlat**2) :: Ham
  real(dp) :: fermi
  real(dp) :: mu, t, Vimp, temp
  real(dp) :: Uinter !interaction potential
  real(dp), dimension(2*Nlat**2) :: evals 
  real(dp), parameter :: tol=1.0d-6 ! tolerance for convergence
  real(dp), parameter :: alpha = 0.0d0 ! update controller
  character(9), parameter :: fmt='(31G16.8)' ! output - make sure that the first two numbers match Nlat!!

  t = 1.0d0 ! hopping defines energy scale
  temp = 0.001d0 ! really k_BT, in units of hopping t
  mu = -1.0d0 ! chemical potential
  Vimp = 5.0d0 ! impurity potential
  Uinter = 1.8d0 ! attractive interaction responsible for the pairing
  
  middle = (Nlat-1)/2 + 1 ! middle point of lattice

  do iter = 1, itermax ! important to have upper iteration bound

     ! must initialize Hamiltonian
     Ham(:,:) = (0.0d0,0.0d0)
  
! set up Hamiltonian
     do ix = 1,Nlat
        do iy = 1,Nlat

           ! initialize pairing potential to random value
           if (iter.eq.1) then
              Delta(ix,iy) = cmplx(ran((ix-1)*Nlat+iy),0.0d0)
           end if

        
           ! chemical potential
           ! electron-like part
           Ham((ix-1)*Nlat+iy,(ix-1)*Nlat+iy) = -mu

           ! hole-like part
           Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2 + (ix-1)*Nlat+iy) = mu

           ! hopping along y direction
           if (iy.eq.Nlat) then ! periodic boundary conditions
              ! electron-like part
              Ham((ix-1)*Nlat+iy,(ix-1)*Nlat+1) = -t 
              
              ! hole-like part
              Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2+(ix-1)*Nlat+1) = t 
           else
              ! electron-like part
              Ham((ix-1)*Nlat+iy,(ix-1)*Nlat+iy+1) = -t
              
              ! hole-like part
              Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2+(ix-1)*Nlat+iy+1) = t
           end if
           
           if (iy.eq.1) then ! periodic boundary conditions
              ! electron-like part
              Ham((ix-1)*Nlat+iy,(ix-1)*Nlat+Nlat) = -t
              
              ! hole-like part
              Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2+(ix-1)*Nlat+Nlat) = t
              
           else
              ! electron-like part
              Ham((ix-1)*Nlat+iy,(ix-1)*Nlat+iy-1) = -t
              
              ! hole-like part
              Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2+(ix-1)*Nlat+iy-1) = t
           end if
           
           ! hopping along x direction
           if (ix.eq.Nlat) then ! periodic boundary conditions
              ! electron-like part
              Ham((ix-1)*Nlat+iy,iy) = -t

              ! hole-like part
              Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2+iy) = t
           else
              ! electron-like part
              Ham((ix-1)*Nlat+iy,(ix)*Nlat+iy) = -t
              
              ! hole-like part
              Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2+(ix)*Nlat+iy) = t
           end if
           
           if (ix.eq.1) then ! periodic boundary conditions
              ! electron-like part
              Ham((ix-1)*Nlat+iy,(Nlat-1)*Nlat+iy) = -t

              ! hole-like part
              Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2+(Nlat-1)*Nlat+iy) = t 
           else
              ! electron-like part
              Ham((ix-1)*Nlat+iy,(ix-2)*Nlat+iy) = -t

              ! hole-like part
              Ham(Nlat**2+(ix-1)*Nlat+iy,Nlat**2+(ix-2)*Nlat+iy) = t 
           end if


           ! pairing potential, upper right-hand block

           Ham((ix-1)*Nlat+iy,Nlat**2+(ix-1)*Nlat+iy) = Delta(ix,iy)

           ! pairing potential, lower left-hand block

           Ham(Nlat**2+(ix-1)*Nlat+iy,(ix-1)*Nlat+iy) = conjg(Delta(ix,iy))
        
        end do
     end do

     ! place impurity at middle of lattice

     Ham((middle-1)*Nlat+middle, (middle-1)*Nlat+middle) = Vimp
     Ham(Nlat**2+(middle-1)*Nlat+middle, Nlat**2+(middle-1)*Nlat+middle) = -Vimp
      
     ! calculate eigenvalues and eigenvectors
     ! NB: Ham is replaced by the matrix of eigenvectors.
     call eigen('V',Ham,evals)

     Deltanew(:,:) = (0.0d0,0.0d0)
     
  ! iterate the gap equation
     do ix=1,Nlat
        do iy=1,Nlat
           do m1=1,2*Nlat**2
              m2 = (ix-1)*Nlat + iy 
              Deltanew(ix,iy) = Deltanew(ix,iy) - Uinter*Ham(m2,m1)*conjg(Ham(m2+Nlat**2,m1))*fermi(evals(m1),temp) ! eq. 61 of notes_1.pdf
           end do
           ! controlled update
           Deltanew(ix,iy) = alpha*Delta(ix,iy) + (1.0d0 - alpha)*Deltanew(ix,iy)

           if (abs(Deltanew(ix,iy) - Delta(ix,iy)).gt.tol) then
              counter = counter + 1
           end if
        end do
     end do

     if ((counter.eq.0).and.(iter.gt.5)) then
        ! self-consistent solution obtained for all sites, exit iteration loop
        write(*,*) "solution converged"
        exit
     else
        write(*,*)  iter, counter, Nlat**2
        write(*,*) Deltanew(1,1)
        counter = 0
        ! update Delta
        Delta = Deltanew
     end if

  end do

  ! output data
  open(unit=1,file='./potential.dat',status='replace')

  do ix = 1, Nlat
     write(1,fmt) real(Deltanew(ix,:))
  end do
  
  close(1)
      
end program swave_sample

function fermi(evals,temp)
  USE nrtype
  IMPLICIT NONE

  real(dp), intent(in) :: evals, temp
  real(dp) :: fermi

  ! we do this odd thing because we want to avoid any kind of overflow error, i.e. the exponential might give a number which is too large!
  if (evals.gt.0.0d0) then
     fermi = exp(-evals/temp)/(1.0d0 + exp(-evals/temp))
  else
     fermi = 1.0d0/(1.0d0 + exp(evals/temp))
  end if
     
end function fermi

! random number generator, used to initialize Delta
FUNCTION ran(idum) 
  use nrtype
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  INTEGER(K4B), INTENT(INOUT) :: idum
  REAL(dp) :: ran
  INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  REAL(dp), SAVE :: am
  INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
  if (idum <= 0 .or. iy < 0) then
     am=nearest(1.0d0,-1.0d0)/IM
     iy=ior(ieor(888889999,abs(idum)),1)
     ix=ieor(777755555,abs(idum))
     idum=abs(idum)+1
  end if
  ix=ieor(ix,ishft(ix,13))
  ix=ieor(ix,ishft(ix,-17))
  ix=ieor(ix,ishft(ix,5))
  k=iy/IQ
  iy=IA*(iy-k*IQ)-IR*k
  if (iy < 0) iy=iy+IM
  ran=am*ior(iand(IM,ieor(ix,iy)),1)
END FUNCTION ran
