! program to test five point tensor coefficients with respect to their
! behaviour when detS and sub-Gram determinant goes to zero.
!
! pentagon (1,2,3,4,5) with s5 > 0, else sj=0
!
!                      has detS = 2*s12*s23*s34*( s15*s45 - s5*s23 )
!
! box      (1,23,4,5)  has detG = 2*s14*( s15*s45 - s5*s23 )
!                               = 2*(s23-s15-s45+s5)*( s15*s45 - s5*s23 )
!
! using momentum parametrisation for 1+4 -> 2+3+5
!
!                          detS = 2*s12*s23*s34*s14*pt5^2
!                          detG = 2*s14^2*pt5^2
!
! where pt5 is the transverse momentum of particle 5 or the system 34
! relative to the beam axis (=z-axis)

! 
! a rotation of 2,3,5 around the z-axis is evaluated to check for
! stability in the limit pt5 -> 0.
!
! output files: demo_a55_dets_sing.txt, demo_a55_dets_sing.dat   
! file to plot demo_a55_dets_sing.dat: plot_demo_A55.gp
!
program main
  !
  use precision_golem
  use matrice_s
  use parametre
  use form_factor_type, only: form_factor
  use form_factor_5p, only: a55
  use cache, only: allocate_cache, clear_cache
  use constante, only: s_null, pi
  use spinor, only: scalar
  !
  implicit none
  !
  type(form_factor) :: res6
  real(ki), dimension(4) :: p1,p2,p3,p4,p5,p2i,p3i,p5i,ptest
  real(ki), dimension(4,4) :: boost,rotation
  real(ki) :: s12,s23,s34,s45,s15,s5,rs,e3,e5,s14
  real(ki) :: t1,t2,theta,sh,ch,phi3,theta3,lambda,s,c,pt,p5z,detg
  integer :: nstep,n
  !
  !
  !open(unit=19,file='error.txt',status='unknown')
  !
  ! Opening of the file containing the results
  !
  open(unit=16,file='demo_a55_dets_sing.dat',status='unknown')
  open(unit=17,file='demo_a55_dets_sing.txt',status='unknown')
  !
  ! These are the entries of the S matrix
  ! They are related to the cuts of the following diagram
  !
  ! p1            p5
  !  \           /
  !   \   (5)   /
  !    |---<---\ (4)
  !    |        \
  !    |         \____ p4
  ! (1)|         /
  !    |        / (3)
  !    |--->---/
  !   /   (2)   \
  !  /           \
  ! p2            p3
  !
  ! S(1,3) = (p2+p3)^2 = s23
  ! S(2,4) = (p3+p4)^2 = s34
  ! S(2,5) = (p1+p2)^2 = s12
  ! S(3,5) = (p4+p5)^2 = s45
  ! S(1,4) = (p1+p5)^2 = s15
  ! S(1,2) = p2^2 = 0
  ! S(2,3) = p3^2 = 0
  ! S(3,4) = p4^2 = 0
  ! S(4,5) = p5^2 = s5
  ! S(1,5) = p1^2 = 0
  !
  !
  call initgolem95(5)
  !  
  !  rs =  sqrt(s14)              ! kinematics normalized such that s14 = 1.0_ki 
  rs = 1.0_ki                               
  theta3 = 0.2_ki*pi                            ! must be in range [0,  pi]
  phi3   = 0.6_ki*2.0_ki*pi                       ! must be in range [0,2*pi]
  s14    = rs**2                               !
  s5     = 0.1_ki*s14                           ! 0 < s5  < s12
  s23    = 0.00003_ki*( sqrt(s14) - sqrt(s5) )**2   ! 0 < s23 < sqrt(s14) - sqrt(s5)
  ! 
  !  
  p1(1) = rs/2.0_ki
  p1(2) =   0.0_ki
  p1(3) =   0.0_ki
  p1(4) = rs/2.0_ki
  !
  p4(1) = rs/2.0_ki
  p4(2) =   0.0_ki
  p4(3) =   0.0_ki
  p4(4) =-rs/2.0_ki
  !
  ! center of mass frame (i) of massless particles 2 and 3
  !  
  !   p2i+p3i = ( sqrt(s23),0,0,0 )
  !
  e3    = sqrt(s23)/2.0_ki   
  p3i(1) =-e3
  p3i(2) =-e3*sin(phi3)*sin(theta3)
  p3i(3) =-e3*cos(phi3)*sin(theta3)
  p3i(4) =-e3*cos(theta3)
  !
  p2i(1) =  p3i(1)
  p2i(2) = -p3i(2)
  p2i(3) = -p3i(3)
  p2i(4) = -p3i(4)
  !  
  ! boost from 2+3 rest frame to pt=0 region
  !
  lambda = s14**2 + s5**2 + s23**2 - 2.0_ki*( s14*s5 + s14*s23 + s23*s5 )
  sh     = sqrt( lambda/4.0_ki/s14/s23 )
  ch     = sqrt( 1.0_ki + sh**2 )
  !
  boost(1,:) = (/   ch, 0.0_ki, 0.0_ki,  -sh /)
  boost(2,:) = (/ 0.0_ki, 1.0_ki, 0.0_ki, 0.0_ki /)
  boost(3,:) = (/ 0.0_ki, 0.0_ki, 1.0_ki, 0.0_ki /)
  boost(4,:) = (/  -sh, 0.0_ki, 0.0_ki,   ch /)
  !
  p2i = matmul( boost,p2i )
  p3i = matmul( boost,p3i )
  !
  p5z   = sqrt( lambda/s14 )/2.0_ki
  e5    = sqrt( s5 + p5z**2 )  
  !
  p5i(1) =-e5
  p5i(2) = 0.0_ki
  p5i(3) = 0.0_ki
  p5i(4) =-p5z
  !
  ! given kinematics leads to scattering singularity,
  ! now rotate 2,3,5 around x-axis.   
  !
  write (17,*)'  n  theta           pt              detG             Re(A55)         Im(A55)' 
  !
  ! 
  ! number of evaluated points:
  !
  nstep = 100
  !
  do n=1,nstep-1   
  ! loop over theta to rotate kinematics around x-axis into a scattering singularity
  ! 
  theta = pi/2.0_ki*( 1.0_ki - real( n,ki )/real(nstep,ki) )
  s = sin(theta)
  c = cos(theta) 
  !
  rotation(1,:) = (/ 1.0_ki, 0.0_ki, 0.0_ki, 0.0_ki /)
  rotation(2,:) = (/ 0.0_ki, 1.0_ki, 0.0_ki, 0.0_ki /)
  rotation(3,:) = (/ 0.0_ki, 0.0_ki,    c,   -s /)
  rotation(4,:) = (/ 0.0_ki, 0.0_ki,    s,    c /)
  !
  p2 = matmul( rotation,p2i )
  p3 = matmul( rotation,p3i )
  p5 = matmul( rotation,p5i ) 
  !
  ! p1+p2+p3+p4+p5=0  
  ! 
  ptest = p1+p2+p3+p4+p5
  !
  !write(*,*) 'Sum of momenta ='
  !write(*,*)ptest
  !
  ! remaining scalar products
  ! 
  s12 = scalar( p1+p2, p1+p2 )
  s34 = scalar( p3+p4, p3+p4 )
  s45 = scalar( p4+p5, p4+p5 )
  s15 = scalar( p1+p5, p1+p5 )
  !
  ! gram determinant of box( 1,2,34,5 ) 
  !
  ! note that: detg = 2.0_ki*s14**2*pt**2
  !
  pt   = p5(3)
  detg = 2.0_ki*s14*( s15*s45 - s23*s5 ) 
  ! write(*,*) detg/2.0_ki/s14**2/pt**2 ! should be equal to 1.0_ki
  !
  !
  ! Definition of the S matrix
  !
  s_mat(1,:) = (/ 0.0_ki, 0.0_ki,s23,s15,0.0_ki /)
  s_mat(2,:) = (/ 0.0_ki, 0.0_ki,0.0_ki,s34,s12 /)
  s_mat(3,:) = (/ s23,0.0_ki,0.0_ki,0.0_ki, s45 /)
  s_mat(4,:) = (/ s15,s34,0.0_ki,0.0_ki,s5 /)
  s_mat(5,:) = (/ 0.0_ki,s12,s45,s5,0.0_ki /)
  !
  ! The inverse of the S matrix is computed numerically, and also all
  ! the related quantities: b coefficients, ......
  call preparesmatrix()
  !
  !
  call cpu_time(t1)
  !
  !
  res6 = a55(1,1,1,1,1,s_null)
  ! 
  call cpu_time(t2)
  !
  write (17,'(i4,6(1x, d15.8))') n,theta,pt,detg,real(res6%c,ki),aimag(res6%c)
  write (16,'(4(1x,f25.5))') pt,detg,real(res6%c,ki),aimag(res6%c)
  !
  !
  end do
  !
  call exitgolem95()
  !
  close(16)
  close(17)
  !close(19)
  !
  write(6,*) 'The result has been written to demo_a55_dets_sing.txt, demo_a55_dets_sing.dat'
  !
end program main
!
