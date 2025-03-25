!
! This program computes the scattering of two massive particles 
! into two light-like particles 
! in a region where the velocity of the ingoing particles
! is small.
!
program main
  !
  use precision_golem
  use matrice_s
  use parametre
  use function_4p2m_adj, only: f4p2m_adj
  use constante, only: pi ! to get the value of pi
  use spinor, only: scalar
  !
  implicit none
  !
  real(ki), dimension(4) :: temp,rest
  real(ki) :: m_sq,pp,ener,theta,s12,s13,s23,s34,s1,s2,x,gr_b
  real(ki), dimension(4) :: p1,p2,p3,p4
  real(ki), dimension(4) :: p23,p34,p12
  integer :: n
  integer :: choix
  !
  choix = 2
  !
  open(unit=8,file='demo_detG.txt',status='unknown')
  !
  open(unit=9,file='demo_detG.dat',status='unknown')
  !
  !
  !
  m_sq = 7.0_ki
  theta = 35.0_ki/180.0_ki*pi
 !
  ! tolerance is the relative precision for the Gaussian integration
  ! it is declared in the file parametre.f90 
  !
  tolerance = 1.d-6
  coupure_4p2m_adj = 5.d-3
  !
  ! Iteration on the three momentum of the ingoing particles
  !
    write (6,*) 'B gets smaller in each iteration'
  !   
    call initgolem95(4)
  !  
  do n=0,30
    !
    write (6,*) 'iteration:',n
    ! pp=M*x
    pp = 0.5d0**n
    ener = sqrt(m_sq + pp*pp)
    x=pp/sqrt(m_sq)
   ! 
    p1 = (/ener,0.0_ki,0.0_ki,pp/)
    p2 = (/ener,0.0_ki,0.0_ki,-pp/)
    p3 = (/ener,0.0_ki,ener*sin(theta),ener*cos(theta)/)
    p4 = (/ener,0.0_ki,-ener*sin(theta),-ener*cos(theta)/)
    !   
    p12 = p1 + p2
    p23 = p2 - p3
    p34 = p3 + p4
    ! 
    !  The S matrix elements
    !
    s34 = scalar(p34,p34)
    s23 = scalar(p23,p23)
    s12 = scalar(p12,p12)
    s13 = -s12-s23 + 2.0_ki*m_sq
    s1  = scalar(p1,p1)
    s2  = scalar(p2,p2)
    !
    ! Definition of the S matrix
    !
    !
    s_mat(1,1) = 0.0_ki
    s_mat(1,2) = s2
    s_mat(1,3) = scalar(p23,p23)
    s_mat(1,4) = s1
    !
    s_mat(2,1) = s_mat(1,2)
    s_mat(2,2) = 0.0_ki
    s_mat(2,3) = 0.0_ki
    s_mat(2,4) = scalar(p12,p12)
    !
    s_mat(3,1) = s_mat(1,3)
    s_mat(3,2) = s_mat(2,3)
    s_mat(3,3) = 0.0_ki
    s_mat(3,4) = 0.0_ki
    !
    s_mat(4,1) = s_mat(1,4)
    s_mat(4,2) = s_mat(2,4)
    s_mat(4,3) = s_mat(3,4)
    s_mat(4,4) = 0.0_ki
    !
    call preparesmatrix
    !
    !  This is minus the coefficient B
    !
    gr_b=(2.0_ki*(sin(theta))**2*x**2)/(m_sq*(1.0_ki+2.0_ki*x**2+2.0_ki*x*sqrt(1.0_ki+x**2)*cos(theta))**2)
    !
    if (choix == 1) then
      !
      temp = f4p2m_adj('n+2',s34,s23,s1,s2,0,0,0,0)
      !
    else if (choix == 2) then
      !
      temp = f4p2m_adj('n+2',s34,s23,s1,s2,0,1,2,2)
      !
    end if
    !
    rest = temp
    !
    write (8,*) 'n=',n
    write (8,*) 'x=',x
    write (6,*) 'x=',x
    write (8,*) '|B|=',gr_b
    write (6,*) '|B|=',gr_b
    write (8,*) 's12=',s12
    write (8,*) 's23=',s23
    write (8,*) 'real part:',rest(3)
    write (8,*) 'imaginary part:',rest(4)
    write (6,*) 'real part:',rest(3)
    write (6,*) 'imaginary part:',rest(4)
    !
    write (9,'(3(2x,f16.14))') x,rest(3),rest(4)
    !
  end do
  !
  call exitgolem95()
  !
  close(8)
  !
end program main
