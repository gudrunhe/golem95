!
! This program computes form factors for five-point functions and related
! (by pinches) diagrams with less external legs.
! The normalisation is as follows:
! We define I_N^n= mu^(4-n) \int d^n k/(i*Pi^(n/2))*func(k,p_i)
! = r_Gam *(P2/eps^2+P1/eps+P0)
! n=4-2*eps
! r_Gam= Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2eps)
! the program gives numbers for P2,P1,P0
!
program main
  !
  use precision_golem
  use matrice_s
  use form_factor_type
  use form_factor_3p ! export all
  use form_factor_4p ! export all
  use form_factor_5p ! export all
  use constante, only: s_null
  use parametre, only: mu2_scale_par
  use form_factor_higher_ranks ! includes form factors for pentagons: A56, ...
  !
  implicit none
  !
  type(form_factor) :: res6
  real(ki) :: t1,t2
  integer :: choix
  !
  ! Opening of the error files
  !
  open(unit=19,file='error_5point.txt',status='unknown')
  !
  ! Opening of the files containing the results
  !
  open(unit=17,file='test5point.txt',status='unknown')
  !
  ! These are the entries of the S matrix
  ! They are related to the cuts of the following diagram
  ! All the momenta are incoming : p1+p2+p3+p4+p5 = 0
  !
  !   p1            p5
  !    \             /
  !     \   (5)    /
  !      |---<---\ (4)
  !      |           \
  !      |            \____ p4
  ! (1) |            /
  !      |           / (3)
  !      |--->---/
  !     /   (2)    \
  !    /             \
  !   p2            p3
  !
  ! S(1,3) = (p2+p3)^2
  ! S(2,4) = (p3+p4)^2
  ! S(2,5) = (p1+p2)^2
  ! S(3,5) = (p4+p5)^2
  ! S(1,4) = (p1+p5)^2
  ! S(1,2) = p2^2
  ! S(2,3) = p3^2
  ! S(3,4) = p4^2
  ! S(4,5) = p5^2
  ! S(1,5) = p1^2
  !
  ! Allocates memory to store the set of initial propagators, the S matrix,
  ! its inverse and the b coefficients.
  ! This call will allocate a derived type s_mat_p object.
  ! Includes calls to allocation_s and initializes the caching system
  !
  call initgolem95(5)
  !
  ! Definition of the S matrix
  !
  s_mat(1,1) = 0.0_ki
  s_mat(1,2) = 0.0_ki
  s_mat(1,3) = -3.0_ki
  s_mat(1,4) = -4.0_ki
  s_mat(1,5) = 0.0_ki
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = 0.0_ki
  s_mat(2,3) = 0.0_ki
  s_mat(2,4) = 6.0_ki
  s_mat(2,5) = 15.0_ki
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = 0.0_ki
  s_mat(3,4) = 0.0_ki
  s_mat(3,5) = 2.0_ki
  !
  s_mat(4,1) = s_mat(1,4)
  s_mat(4,2) = s_mat(2,4)
  s_mat(4,3) = s_mat(3,4)
  s_mat(4,4) = 0.0_ki
  s_mat(4,5) = 0.0_ki
  !
  s_mat(5,1) = s_mat(1,5)
  s_mat(5,2) = s_mat(2,5)
  s_mat(5,3) = s_mat(3,5)
  s_mat(5,4) = s_mat(4,5)
  s_mat(5,5) = 0.0_ki
  !
  ! This call fills the internal array s_mat_r.
  ! It also assigns the integers in s_mat_p, which encode the positions
  ! of complex mass entries and zero mass entries. It includes call to init_invs
  !
  call preparesmatrix()
  !
  write (*,*) 'Choose what the program should compute:'
  write (*,*) '0) form factor for five-point function, rank 0'
  write (*,*) '1) form factor for five-point function, rank 3 (z1*z2*z4)'
  write (*,*) '2) form factor for five-point function, rank 5 (z1*z2*z3*z4*z5)'
  write (*,*) '3) form factor for diagram with propagator 3 pinched, rank 0'
  write (*,*) '4) form factor for diagram with propagators 1 and 4 pinched, rank 0'
  write (*,*) '5) form factor for five-point function, rank 6 (z1^2*z2*z3*z4*z5)'
  read (*,*) choix
!
  if (choix == 0) then
    !
    write (*,*) 'calculating form factor for 5-point function rank 0'
    !
  else if   (choix == 1) then
    !
    write (*,*) 'calculating form factor A_124 for 5-point function rank 3'
    !
  else if   (choix == 2) then
    !
    write (*,*) 'calculating form factor A_12345 for 5-point function rank 5'
    !
  else if   (choix == 3) then
    !
    write (*,*) 'calculating form factor for a box stemming from '
    write (*,*) 'the pinch of propagator 3 of a 5-point funct., rank0'
    !
  else if   (choix == 4) then
    !
    write (*,*) 'calculating form factor for a triangle stemming from '
    write (*,*) 'the pinch of propagators 1 and 4 of a 5-point funct., rank 0'
    !
  else if (choix == 5) then
    !
    write (*,*) 'calculating form factor A_112345 for 5-point function rank 6'
  end if
  !
  write (*,*) 'The result has been written to the file test5point.txt'
  !
  call cpu_time(t1)
  !
  ! To change the value of mu^2 (in GeV) (set to 1. by default)
  ! uncomment this line
  !
  !mu2_scale_par = 12._ki
  !
  write (17,*) 'The kinematics is:'
  write (17,*) ''
  write (17,*) ' p1           p5    '
  write (17,*) '   \           /        '
  write (17,*) '    \   (5    /         '
  write (17,*) '     |---<---\ (4)     '
  write (17,*) '     |        \        '
  write (17,*) '     |         \____ p4'
  write (17,*) '  (1)|         /       '
  write (17,*) '     |        / (3)   '
  write (17,*) '     |--->---/         '
  write (17,*) '    /   (2)  \         '
  write (17,*) '   /          \        '
  write (17,*) ' p2            p3    '
  write (17,*) ''
  write (17,*) 'p1+p2+p3+p4+p5 = 0'
  write (17,*) ''
  write (17,*) ' S(1,3) = (p2+p3)^2=',s_mat(1,3)
  write (17,*) ' S(2,4) = (p3+p4)^2=',s_mat(2,4)
  write (17,*) ' S(2,5) = (p1+p2)^2=',s_mat(2,5)
  write (17,*) ' S(3,5) = (p4+p5)^2=',s_mat(3,5)
  write (17,*) ' S(1,4) = (p1+p5)^2=',s_mat(1,4)
  write (17,*) ' S(1,2) = p2^2=',s_mat(1,2)
  write (17,*) ' S(2,3) = p3^2=',s_mat(2,3)
  write (17,*) ' S(3,4) = p4^2=',s_mat(3,4)
  write (17,*) ' S(4,5) = p5^2=',s_mat(4,5)
  write (17,*) ' S(1,5) = p1^2=',s_mat(1,5)
  write (17,*) '(mu)^2 =',mu2_scale_par
  write (17,*) ''
  !
  if (choix == 0) then
    !
    ! form factor for five-point function, rank 0
    !
    res6 = a50(s_null)
    !
  else if (choix == 1) then
    !
    ! form factor for five-point function, rank 3
    !
    res6 = a53(1,2,4,s_null)
    !
  else if (choix == 2) then
    !
    ! form factor for five-point function, rank 5
    !
    res6 = a55(1,2,3,4,5,s_null)
    !
  else if (choix == 3) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagator 3 is pinched
    !
    !
    !   p1            p4
    !    \             /
    !     \   (5)    /
    !      |---<---|
    !      |          |
    ! (1) |          | (4)
    !      |          |
    !      |--->---|
    !     /   (2)     \
    !    /              \
    !  p2            p3+p4
    !
    write (17,*) 'Since the propagator 3 is pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) ' p1             p5'
    write (17,*) '   \           /  '
    write (17,*) '    \   (5)   /   '
    write (17,*) '     |---<---|    '
    write (17,*) '     |       |    '
    write (17,*) '  (1)|       |(4) '
    write (17,*) '     |       |    '
    write (17,*) '     |--->---|    '
    write (17,*) '    /   (2)   \   '
    write (17,*) '   /           \  '
    write (17,*) ' p2            p3+p4'
    write (17,*) ''
    !
    !
    res6 = a40( (/3/)  )
    !
  else if (choix == 4) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagators 1 and 4 are pinched
    !
    !           |
    !           | p1+p2
    !           |
    !           /\
    !          /  \
    !     (2) /    \ (5)
    !        /      \
    !       /--->----\
    !      /   (3)    \
    ! p3  /            \ p4+p5
    !
    write (17,*) 'Since the propagators 1 and 4 are pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) '            |             '
    write (17,*) '            | p1+p2       '
    write (17,*) '            |             '
    write (17,*) '            /\            '
    write (17,*) '           /  \           '
    write (17,*) '     (2)  /    \ (5)      '
    write (17,*) '         /      \         '
    write (17,*) '        /--->----\        '
    write (17,*) '       /   (3)    \       '
    write (17,*) '   p3 /            \ p4+p5'
    write (17,*) ''
    !
    !
    res6 = a30( (/1,4/) )
    !
   else if (choix == 5) then
    !
    ! form factor for five-point function, rank 6
    !
    res6 = a56(1,1,2,3,4,5,s_null)
    !
  end if
  !
  call cpu_time(t2)
  !!
  write (17,*) 'normalisation:'
  write (17,*) 'defining I_N^n= mu^(4-n) \int d^n k/(i*Pi^(n/2))*func(k,p_i)'
  write (17,*) '= r_Gam *(P2/eps^2+P1/eps+P0),'
  write (17,*) 'n = 4-2*eps,'
  write (17,*) 'r_Gam = Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2eps)'
  write (17,*) 'the program gives numbers for P2,P1,P0'
  write (17,*) ''
  write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
  write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki),aimag(res6%b)
  write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%c,ki),aimag(res6%c)
  write (17,*) ''
  write (17,*) 'CPU time=',t2-t1
  !
  ! routine to free the cache and allocated memory
  !
  call exitgolem95()
  !
  !
  close(17)
  close(19)
  !
end program main
