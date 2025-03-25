!
! This program computes form factors for the six-point functions and related
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
  use form_factor_6p ! export all
  use form_factor_higher_ranks ! export all
  use constante, only: s_null
  use spinor, only: scalar  ! to get the function scalar
  use parametre, only: mu2_scale_par
  !
  implicit none
  !
  type(form_factor) :: res6
  real(ki) :: t1,t2,tg
  real(ki), dimension(4) :: p1,p2,p3,p4,p5,p6
  real(ki), dimension(4) :: p12,p23,p34,p45,p56,p61,p123,p234,p345,p123456
  integer :: choix,i
  !
  ! Opening of the error file
  !
  open(unit=19,file='error_6point.txt',status='unknown')
  !
  ! Opening of the files containing the results
  !
  open(unit=17,file='test6point.txt',status='unknown')
  !
  !
  ! This is the entries of the S matrix
  ! defined as S(i,j) = (q_i - q_j)^2 where the q's are
  ! the momentum flowing in the propagators
  ! It is related to the cuts of the following diagram
  ! All the momenta are incoming : p1+p2+p3+p4+p5+p6 = 0
  !
  !         p6            p5
  !          \           /
  !           \   (5)   /
  !            /---<---\ (4)
  !       (6) /         \
  !   p1 ____/           \____ p4
  !          \           /
  !       (1) \         / (3)
  !            \--->---/
  !           /   (2)   \
  !          /           \
  !         p2            p3
  !
  ! S(1,3) = (p2+p3)^2
  ! S(1,4) = (p2+p3+p4)^2
  ! S(1,5) = (p1+p6)^2
  ! S(2,4) = (p3+p4)^2
  ! S(2,5) = (p3+p4+p5)^2
  ! S(2,6) = (p1+p2)^2
  ! S(3,5) = (p4+p5)^2
  ! S(3,6) = (p4+p5+p6)^2
  ! S(4,6) = (p5+p6)^2
  ! S(1,2) = p2^2
  ! S(2,3) = p3^2
  ! S(3,4) = p4^2
  ! S(4,5) = p5^2
  ! S(5,6) = p6^2
  ! S(1,6) = p1^2
  !
  ! we define a set of four momenta which obey the different constraints
  p1 = (/ 0.5_ki , 0.0_ki , 0.0_ki , 0.5_ki /)
  p2 = (/ 0.5_ki , 0.0_ki , 0.0_ki , -0.5_ki /)
 p3= (/  -0.185986759664737864811227397267961_ki,  4.452834165071030414484400000000003E-0002_ki,&
 &-4.312345594415791655285000000000000E-0002_ki, 1.112066130492704585375100000000000E-0002_ki/)
 p4= (/  -0.326011971996721328093116000000000_ki, -8.549197530383799426001400000000000E-0002_ki,&
 & -9.126117670433854522915600000000000E-0002_ki,  0.244992501253361860680258000000000_ki /)
 p5= (/  -0.233498949364477610783116773003482_ki, -0.116776431857305545980452000000000_ki,&
 & 0.191553925634582344894596000000000_ki,  -6.474656663462434458278200000000000E-0002_ki /)
 p6= (/  -0.254502318974063196312539829728557_ki,   0.157740065510433236095622000000000_ki,&
 & -5.716929298608588311259000000000000E-0002_ki, -0.191366595923664561951227000000000_ki /)
! p3 = -(/ 0.19178191094778038_ki , 0.12741179719516801_ki , 0.08262476614744381_ki , 0.11713105190921771_ki /)
!  p4 = -(/ 0.33662712284553753_ki , -0.06648281097623857_ki , -0.3189378514746887_ki , -0.08471424069583446_ki /)
!  p5 = -(/ 0.21604814388379073_ki , -0.20363139428835617_ki , 0.044157623555325_ki , 0.0571065672034082_ki /)
!  p6 = -(/ 0.2555428223228916_ki , 0.1427024080694266_ki , 0.19215546177191994_ki , -0.08952337841679145_ki /)
  !
  p12 = p1 + p2
  p23 = p2 + p3
  p234 = p2+p3+p4
  p61 = p1 + p6
  p34 = p3 + p4
  p345 = p34 + p5
  p45 = p4 + p5
  p123 = p12 + p3
  p56 = p5 + p6
  !
  ! Allocates memory to store the set of initial propagators, the S matrix,
  ! its inverse and the b coefficients.
  ! This call will allocate a derived type s_mat_p object.
  ! Includes calls to allocation_s and initializes the caching system
  !
  call initgolem95(6)
  !
  !
  ! Definition of the S matrix
  !
  s_mat(1,1) = 0.0_ki
  s_mat(1,2) = scalar(p2,p2)
  s_mat(1,3) = scalar(p23,p23)
  s_mat(1,4) = scalar(p234,p234)
  s_mat(1,5) = scalar(p61,p61)
  s_mat(1,6) = scalar(p1,p1)
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = 0.0_ki
  s_mat(2,3) = scalar(p3,p3)
  s_mat(2,4) = scalar(p34,p34)
  s_mat(2,5) = scalar(p345,p345)
  s_mat(2,6) = scalar(p12,p12)
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = 0.0_ki
  s_mat(3,4) = scalar(p4,p4)
  s_mat(3,5) = scalar(p45,p45)
  s_mat(3,6) = scalar(p123,p123)
  !
  s_mat(4,1) = s_mat(1,4)
  s_mat(4,2) = s_mat(2,4)
  s_mat(4,3) = s_mat(3,4)
  s_mat(4,4) = 0.0_ki
  s_mat(4,5) = scalar(p5,p5)
  s_mat(4,6) = scalar(p56,p56)
  !
  s_mat(5,1) = s_mat(1,5)
  s_mat(5,2) = s_mat(2,5)
  s_mat(5,3) = s_mat(3,5)
  s_mat(5,4) = s_mat(4,5)
  s_mat(5,5) = 0.0_ki
  s_mat(5,6) = scalar(p6,p6)
  !
  s_mat(6,1) = s_mat(1,6)
  s_mat(6,2) = s_mat(2,6)
  s_mat(6,3) = s_mat(3,6)
  s_mat(6,4) = s_mat(4,6)
  s_mat(6,5) = s_mat(5,6)
  s_mat(6,6) = 0.0_ki
  !
  ! This call fills the internal array s_mat_r.
  ! It also assigns the integers in s_mat_p, which encode the positions
  ! of complex mass entries and zero mass entries. It includes call to init_invs
  !
  call preparesmatrix()
  write (*,*) 'Choose what the program has to compute:'
  write (*,*) '0) form factor for six-point function, rank 0'
  write (*,*) '1) form factor for six-point function, rank 4 (z1^2*z2*z3) '
  write (*,*) '2) form factor for diagram where propagator 3 is pinched, rank 0'
  write (*,*) '3) form factor for diagram where props. 2,5 are pinched, rank 0'
  write (*,*) '4) form factor for diagram where props. 2,4,6 are pinched, rank 0'
  write (*,*) '5) form factor for six-point function, rank 7 (z1^2*z2^2*z3^2*z4) '
  read (*,*) choix
  write (*,*) 'The result has been written to the file test6point.txt'
  !
  !
  call cpu_time(t1)
  !
  !
  ! To change the value of mu^2 (in GeV) (set to 1. by default)
  ! uncomment this line
  !
  !mu2_scale_par = 12._ki
  !
  write (17,*) 'The kinematics is:'
  write (17,*) ''
  write (17,*) '          p6            p5      '
  write (17,*) '           \           /        '
  write (17,*) '            \   (5)   /         '
  write (17,*) '            /---<----\         '
  write (17,*) '       (6) /          \ (4)       '
  write (17,*) ' p1 ______/            \____ p4 '
  write (17,*) '          \            /        '
  write (17,*) '       (1) \          / (3)     '
  write (17,*) '            \--->----/          '
  write (17,*) '            /   (2)   \         '
  write (17,*) '           /           \        '
  write (17,*) '           p2            p3      '
  write (17,*) ''
  write (17,*) ' S(1,3) = (p2+p3)^2=',s_mat(1,3)
  write (17,*) ' S(1,4) = (p2+p3+p4)^2=',s_mat(1,4)
  write (17,*) ' S(1,5) = (p1+p6)^2=',s_mat(1,5)
  write (17,*) ' S(2,4) = (p3+p4)^2=',s_mat(2,4)
  write (17,*) ' S(2,5) = (p3+p4+p5)^2=',s_mat(2,5)
  write (17,*) ' S(2,6) = (p1+p2)^2=',s_mat(2,6)
  write (17,*) ' S(3,5) = (p4+p5)^2=',s_mat(3,5)
  write (17,*) ' S(3,6) = (p4+p5+p6)^2=',s_mat(3,6)
  write (17,*) ' S(4,6) = (p5+p6)^2=',s_mat(4,6)
  write (17,*) ' S(1,2) = p2^2=',s_mat(1,2)
  write (17,*) ' S(2,3) = p3^2=',s_mat(2,3)
  write (17,*) ' S(3,4) = p4^2=',s_mat(3,4)
  write (17,*) ' S(4,5) = p5^2=',s_mat(4,5)
  write (17,*) ' S(5,6) = p6^2=',s_mat(5,6)
  write (17,*) ' S(1,6) = p1^2=',s_mat(6,1)
  write (17,*) '(mu)^2 =',mu2_scale_par
  write (17,*) ''
  !
  if (choix == 0) then
    !
    ! form factor for six-point function, rank 0
    !
    res6 =  a60(s_null)
    !
  else if (choix == 1) then
    !
    ! form factor for six-point function, rank 4
    !
    res6 = a64(1,1,2,3,s_null)
    !
  else if (choix == 2) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagator 3 is pinched
    !
    !
    ! p1            p6
    !  \           /
    !   \   (6)   /
    !    |---<---\ (5)
    !    |        \
    !    |         \____ p5
    ! (1)|         /
    !    |        / (4)
    !    |--->---/
    !   /   (2)   \
    !  /           \
    ! p2            p3+p4
    !
    !
    write (17,*) 'Since the propagator 3 is pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) ' p1            p6     '
    write (17,*) '   \           /       '
    write (17,*) '    \   (6)   /        '
    write (17,*) '     |---<----\ (5)     '
    write (17,*) '     |         \        '
    write (17,*) '     |          \____ p5'
    write (17,*) ' (1) |          /       '
    write (17,*) '     |         / (4)    '
    write (17,*) '     |--->----/         '
    write (17,*) '    /   (2)   \        '
    write (17,*) '   /           \       '
    write (17,*) ' p2            p3+p4  '
    !
    !
    res6 = a50((/3/))
    !
  else if (choix == 3) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagators 2 and 5 are pinched
    !
    !
    ! p1            p5+p6
    !  \           /
    !   \   (6)   /
    !    |---<---|
    !    |       |
    ! (1)|       |(4)
    !    |       |
    !    |--->---|
    !   /   (3)   \
    !  /           \
    ! p2+p3          p4
    !
    write (17,*) 'Since the propagators 2 and 5 are pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) '   p1              p5+p6'
    write (17,*) '     \             /     '
    write (17,*) '      \     (6)   /      '
    write (17,*) '       |---<-----|       '
    write (17,*) '       |         |       '
    write (17,*) '  (1)  |         | (4)   '
    write (17,*) '       |         |       '
    write (17,*) '       |--->-----|       '
    write (17,*) '      /     (3)   \      '
    write (17,*) '     /             \     '
    write (17,*) '  p2+p3            p4   '
    write (17,*) ''
    !
    !
    res6 = a40((/2,5/))
    !
  else if (choix == 4) then
    !
    ! form factor for pinched diagram, rank 0
    ! the propagator 2, 4 and 6 are pinched
    !
    !           |
    !           | p1+p6
    !           |
    !           /\
    !          /  \
    !     (1) /    \ (5)
    !        /      \
    !       /--->----\
    !      /   (3)    \
    !     /            \
    ! p2+p3            p4+p5
    !
    write (17,*) 'Since the propagators 2, 4 and 6 are pinched'
    write (17,*) 'the reduced kinematics is:'
    write (17,*) ''
    write (17,*) '           |             '
    write (17,*) '           | p1+p6       '
    write (17,*) '           |             '
    write (17,*) '           /\            '
    write (17,*) '          /  \           '
    write (17,*) '     (1) /    \ (5)      '
    write (17,*) '        /      \         '
    write (17,*) '       /--->----\        '
    write (17,*) '      /   (3     \       '
    write (17,*) '     /            \      '
    write (17,*) ' p2+p3            p4+p5  '
    write (17,*) ''
    !
    !
    res6 = a30((/2,4,6/))
    !
  else if (choix == 5) then
    !
    ! form factor for six-point function, rank 7
    !
    res6 = a67(1,1,2,2,3,3,4,s_null)
  end if
  !
  call cpu_time(t2)
  !
  write (17,*) ''
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
  !
  write (17,*) 'CPU time=',t2-t1
  !
  !
  ! routine to free the cache and allocated memory
  !
  call exitgolem95()
  !
  close(17)
  close(19)
  !
end program main
