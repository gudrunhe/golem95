! computes the Susy gg to Hbb example to scan Landau singularity
program main
  !
  use precision_golem
  use matrice_s
  use form_factor_type
  use form_factor_4p 
  use parametre
  use constante
  !
  implicit none
  !
  type(form_factor) :: res6
  real(ki) :: ti1,ti2
  real(ki) :: p1sq,p2sq,p3sq,p4sq
  real(ki) :: s_var,t_var,s2,s2min,s2max,s,sqrts,s1
  integer :: Nsteps,istep
  real(ki) :: mu02,Gamsquark,Gammaxi
  complex(ki) :: m1sq,m2sq,m3sq,m4sq
  real(ki) :: msquark,mxi,mH
  !
  ! width:
   Gamsquark=3.5_ki
   Gammaxi=1.5_ki
  !Gamsquark=0._ki
  !Gammaxi=0._ki
  !
  msquark=800._ki
  mxi=200._ki
  !
  mH=450._ki
  ! 
  sqrts = 1700._ki
  s=sqrts**2
  s1 = 2.0_ki * (msquark*msquark + mxi*mxi)
  !
  !s2min = mH**2*s/s1+1000._ki
  !s2max = mH**2+s-s1
  s2min = 900._ki**2
  s2max=1200._ki**2
  Nsteps = 200
  !
  !
  open(unit=19,file='errorLandau.txt',status='unknown')
  !
  ! Opening of the files containing the results
  !
  !open(unit=16,file='SusyLandau.dat',status='unknown')
    open(unit=16,file='SusyLandauC.dat',status='unknown')
  open(unit=17,file='SusyLandau.txt',status='unknown')
  !
  t_var = s1
  !
  p1sq = s
  p2sq = 0.0_ki
  p3sq = mH**2
  p4sq = 0.0_ki
  !    m1sq=msquark**2
  m1sq=cmplx(msquark**2,-msquark*Gamsquark,ki)
  m2sq=cmplx(mxi**2,-mxi*Gammaxi,ki)
  !    m2sq=mxi**2
  m3sq=m2sq
  m4sq=m1sq
  !
  call cpu_time(ti1)
  !
  call initgolem95(4)
  !
  do istep = 0, Nsteps
     !  
     s2 = s2min + istep * (s2max - s2min) / real(Nsteps, ki)
     s_var = s2
     !
     ! Definition of the S matrix
     !
     s_mat(1,1) = -2._ki*m1sq
     s_mat(1,2) = p2sq-m1sq-m2sq
     s_mat(1,3) = t_var-m1sq-m3sq
     s_mat(1,4) = p1sq-m1sq-m4sq
     !
     s_mat(2,1) = s_mat(1,2)
     s_mat(2,2) = -2._ki*m2sq
     s_mat(2,3) = p3sq-m2sq-m3sq
     s_mat(2,4) = s_var-m2sq-m4sq
     !
     s_mat(3,1) = s_mat(1,3)
     s_mat(3,2) = s_mat(2,3)
     s_mat(3,3) = -2._ki*m3sq
     s_mat(3,4) = p4sq-m3sq-m4sq
     !
     s_mat(4,1) = s_mat(1,4)
     s_mat(4,2) = s_mat(2,4)
     s_mat(4,3) = s_mat(3,4)
     s_mat(4,4) = -2._ki*m4sq
     !
     !
     call preparesmatrix()
     !
     ! To change the value of mu^2 (in GeV) (set to 1. by default)
     ! uncomment this line
     !
     !mu2_scale_par = 12._ki
     mu02 = mu2_scale_par
     !
     ! Result for the scalar box integral in n dimension
     !
     res6 =  a40(s_null)
     !
     !
     write (16,'(F16.10,6(1x,F26.16))') sqrt(s2),  real(res6%c,ki),aimag(res6%c), &
          & real(res6%b,ki),aimag(res6%b), &
          & real(res6%a,ki),aimag(res6%a) 
     ! The results are written to the file with unit 17
     !
     write (17,*) ''
     write (17,'("1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
     write (17,'("+ 1/epsilon * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki),aimag(res6%b)
     write (17,'("+ (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%c,ki),aimag(res6%c)
     !
     write (6,'("(",e16.10,1x,"+ I*",1x,e16.10,")")"') real(res6%c,ki),aimag(res6%c)
     write (17,*) ''
     !
     !
  end do  ! end loop over s2
  !
  call cpu_time(ti2)
  !
  write (17,*) 'CPU time=',ti2-ti1
  !
  close(16)
  close(17)
  close(19)
  !
  call exitgolem95()
  !
end program main
