!
! This program computes the n, n+2 and n+4 dimensional four-point functions
! with or without Feynman parameters in the numerator,
! from the four-point form factors. This program can be
! used as a test to check this type of integrals.
! The normalisation is as follows:
! We define I_N^n= mu^(4-n) \int d^n k/(i*Pi^(n/2))*func(k,p_i)
! = r_Gam *(P2/eps^2+P1/eps+P0)
! n=4-2*eps
! r_Gam= Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2eps)
! the program gives numbers for P2,P1,P0
! the default value for the scale \mu is 1,
! it can be changed by the user by defining mu2_scale_par = xxx
!
program main
  !
  use precision_golem
  use matrice_s
  use form_factor_type
  use form_factor_4p
  use form_factor_3p
  use constante
  use generic_function_4p
  use function_4p1m
  use function_4p2m_opp
  use function_4p2m_adj
  use function_4p3m
  use function_4p4m
  use parametre
  use generic_function_np ! module with the generic four-point higher-rank functions
  use form_factor_higher_ranks ! corresponding form factors

  !use parametre, only: mu2_scale_par
  !"only" statements are recommended if you need to export only certain functions
  !
  implicit none
  !
  type(form_factor) :: res6,res6a
  complex(ki), dimension(3) :: verif1
  real(ki), dimension(6) :: verif2
  real(ki) :: ti1,ti2
  real(ki) :: p1sq,p2sq,p3sq,p4sq,m1sq,m2sq,m3sq,m4sq
  real(ki) :: s_var,t_var
  integer :: choix,choix_kinem
  real(ki) :: mu2,mu02,lmu2
  !
  write (*,*) 'Choose from the following kinematics:'
  write (*,*) '1) all external legs light-like, no internal masses'
  write (*,*) '2) one external leg not light-like, no internal masses'
  write (*,*) '3) two opposite external legs not light-like, no internal masses'
  write (*,*) '4) two adjacent external legs not light-like, no internal masses'
  write (*,*) '5) three external legs not light-like, no internal masses'
  write (*,*) '6) four external legs not light-like, no internal masses'
  write (*,*) '7) two external legs not light-like, internal masses, IR divergent'
  !(example QL Box 8)
  write (*,*) '8) three external legs not light-like, internal masses, IR divergent'
  !(example QL Box 13)
  read (*,*) choix_kinem
  !
  !
  ! Opening of the error files
  !
  !open(unit=19,file='error_4point.txt',status='unknown')
  !
  ! Opening of the files containing the results
  !
  open(unit=17,file='test4point.txt',status='unknown')
  !
  ! choose to calculate all or rational part only (default is set to 'tot')
  ! rat_or_tot_par='rat'
  !
  ! These are the entries of the S matrix
  ! They are related to the cuts of the following diagram
  ! S(1,3) = (p1+p4)^2 = t_var
  ! S(2,4) = (p1+p2)^2 = s_var
  ! S(1,2) = p2^2 = p2sq
  ! S(2,3) = p3^2 = p3sq
  ! S(3,4) = p4^2 = p4sq
  ! S(4,1) = p1^2 = p1sq
  !
  !  p1            p4
  !    \           /
  !     \   (4)   /
  !      |---<---|
  !      |       |
  ! (1)  |       | (3)
  !      |       |
  !      |--->---|
  !     /   (2)   \
  !    /           \
  !  p2             p3
  !
  !
  ! Allocates memory to store the set of initial propagators, the S matrix,
  ! its inverse and the b coefficients.
  ! This call will allocate a derived type s_mat_p object.
  ! Includes calls to allocation_s and initializes the caching system
  !
  call initgolem95(4)
  !
  m1sq = zero
  m2sq = zero
  m3sq = zero
  m4sq = zero
  !
  if (choix_kinem == 1) then
    !
    p1sq = 0.0_ki
    p2sq = 0.0_ki
    p3sq = 0.0_ki
    p4sq = 0.0_ki
    !
  else if (choix_kinem == 2) then
    !
    p1sq = 0.0_ki
    p2sq = 0.0_ki
    p3sq = 0.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 3) then
    !
    p1sq = 0.0_ki
    p2sq = 50.0_ki
    p3sq = 0.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 4) then
    !
    p1sq = 0.0_ki
    p2sq = 0.0_ki
    p3sq = 50.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 5) then
    !
    p1sq = 0.0_ki
    p2sq = 50.0_ki
    p3sq = 80.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 6) then
    !
    p1sq = 20.0_ki
    p2sq = 50.0_ki
    p3sq = 80.0_ki
    p4sq = 60.0_ki
    !
  else if (choix_kinem == 7) then
    !
    p1sq = 0.0_ki
    p2sq = 0.0_ki
    p3sq = -80.0_ki
    p4sq = -60.0_ki
    m1sq=0.0_ki
    m2sq=0.0_ki
    m3sq=20.0_ki
    m4sq=0.0_ki
    !
  else if (choix_kinem == 8) then
    !
    p1sq = 0.0_ki
    p2sq = 50.0_ki
    p3sq = -80.0_ki
    p4sq = -60.0_ki
    m1sq=0.0_ki
    m2sq=5.0_ki
    m3sq=10.0_ki
    m4sq=0.0_ki
    !
  end if
  !
  s_var = 200.0_ki
  t_var = -123.0_ki
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
  ! This call fills the internal array s_mat_r.
  ! It also assigns the integers in s_mat_p, which encode the positions
  ! of complex mass entries and zero mass entries. It includes call to init_invs
  !
  call preparesmatrix()
  !
  write (*,*) 'Choose what the program should compute:'
  write (*,*) '0) scalar four-point function in n dimensions'
  write (*,*) '1) four-point function in n dimensions with one Feynman parameter (z1)'
  write (*,*) '2) four-point function in n dimensions with two Feynman parameters (z1*z4)'
  write (*,*) '3) four-point function in n dimensions with three Feynman parameters (z1^2*z3)'
  write (*,*) '4) four-point fctn. in n dimensions with four Feynman parameters (z1*z2*z3*z4)'
  write (*,*) '5) scalar four-point function in n+2 dimensions'
  write (*,*) '6) four-point function in n+2 dimensions with two Feynman parameters (z1*z2)'
  write (*,*) '7) scalar four-point function in n+4 dimensions'
  write (*,*) '8) the mu dependence'
  write (*,*) '9) four-point fctn. in n dimensions with five Feynman parameters (z1*z2^2*z3^2)'
  write (*,*) '10) four-point fctn. in n+2 dimensions with three Feynman parameters (z1*z2^2)'
  write (*,*) '11) four-point fctn. in n+4 dimensions with one Feynman parameters (z1)'
  write (*,*) '12) scalar four-point fctn. in n+6 dimensions'

  read (*,*) choix
  !
  if (choix == 0) then
    !
    write (*,*) 'calculating n-dim scalar box with '
    !
  else if   (choix == 1) then
    !
    write (*,*) 'calculating n-dim rank one (z1) 4-point fctn. with'
    !
  else if   (choix == 2) then
    !
    write (*,*) 'calculating n-dim rank two (z1*z4) 4-point fctn. with'
    !
  else if   (choix == 3) then
    !
    write (*,*) 'calculating n-dim rank three (z1^2*z3) 4-point fctn. with'
    !
  else if   (choix == 4) then
    !
    write (*,*) 'calculating n-dim rank four (z1*z2*z3*z4) 4-point fctn. with'
    !
  else if   (choix == 5) then
    !
    write (*,*) 'calculating (n+2)-dim scalar 4-point fctn. with'
    !
  else if   (choix == 6) then
    !
    write (*,*) 'calculating (n+2)-dim rank two (z1*z2) 4-point fctn. with'
    !
  else if (choix == 7) then
    !
    write (*,*) 'calculating (n+4)-dim scalar 4-point function with'
    !
  else if (choix == 8) then
    !
    write (*,*) 'example of a test on the renormalisation scale mu dependence'
    !
  else if   (choix == 9) then
    !
    write (*,*) 'calculating n-dim rank five (z1*z2^2*z3^2) 4-point fctn. with'
    !
  else if   (choix == 10) then
    !
    write (*,*) 'calculating n+2-dim rank three  (z1*z2^2) 4-point fctn. with'
    !
  else if   (choix == 11) then
    !
    write (*,*) 'calculating n+4-dim rank one (z1) 4-point fctn. with'
    !

  end if
  !
  if (choix_kinem == 1) then
    !
    write (*,*) 'no off-shell leg'
    !
  else if (choix_kinem == 2) then
    !
    write (*,*) 'one off-shell leg'
    !
  else if (choix_kinem == 3) then
    !
    write (*,*) 'two opposite off-shell legs'
    !
  else if (choix_kinem == 4) then
    !
    write (*,*) 'two adjacent off-shell legs'
    !
  else if (choix_kinem == 5) then
    !
    write (*,*) 'three off-shell legs'
    !
  else if (choix_kinem == 6) then
    !
    write (*,*) 'four off-shell legs'
    !
  else if (choix_kinem == 7) then
    !
    write (*,*) 'two off-shell legs and one internal mass'
    !
  else if (choix_kinem == 8) then
    !
    write (*,*) 'three off-shell legs and two internal masses'
    !
  end if
  !
  write (*,*) 'The result has been written to the file test4point.txt'
  !
  !
  call cpu_time(ti1)
  !
  ! To change the value of mu^2 (in GeV) (set to 1. by default)
  ! uncomment this line
  !
  !mu2_scale_par = 12._ki
  mu02 = mu2_scale_par
  !
  if (choix == 0) then
    !
    ! Result for the scalar integral in n dimension
    !
    res6 =  a40(s_null)
    !
    ! the labels 1,2,..,4 correspond to Feynman parameter z1,z2,...,z4
  else if (choix == 1) then
    !
    ! Results for integrals in n dimension with one Feynman parameter
    ! at the numerator: z1
    !
    res6 =  - a41(1,s_null)
    !
  else if (choix == 2) then
    !
    ! Results for integrals in n dimension with two Feynman parameters
    ! at the numerator: z1*z4
    !
    res6 =  a42(1,4,s_null)
    !
  else if (choix == 3) then
    !
    ! Results for integrals in n dimension with three Feynman parameters
    ! at the numerator: z1^2*z3
    !
    res6 =  -a43(1,1,3,s_null)
    !
  else if (choix == 4) then
    !
    ! Results for integrals in n dimension with four Feynman parameters
    ! at the numerator: z1*z2*z3*z4
    !
    res6 =  a44(1,2,3,4,s_null)
    !
  else if (choix == 5) then
    !
    ! Results for integrals in n+2 dimension with no Feynman parameters
    ! in the numerator:
    ! There are two ways to compute a n+2 dimension four-point function,
    ! either use the generic function f4p_np2 (in module generic_function_4p),
    ! there are two arguments: the kinematics matrix and the set of the four unpinched
    ! propagators, in addition there are optional arguments for the label of the Feynman
    ! parameters
    ! or use the dedicated function for the corresponding kinematics, the label of the Feynman
    ! parameters is mandatory (if there is none put 0)
    ! In the case of the dedicated function, the kinematics is fixed, for example, for
    ! the one mass four-point function, the external legs are labelled as specified in
    ! picture at the beginning and the massive leg is assumed to be p4. If it is not the case,
    ! the user must do the crossing by hand. On the contrary, in the case of the generic function
    ! the program performs automatically the crossing, the user just has to provide the S matrix.
    ! Concerning the set of unpinched propagators, the S matrix defines the propagators because
    ! S_{ij} = (q_i - q_j)^2 (where the qi are the momenta flowing through the propagator i).
    ! It can happen that one wants to compute a four-point function related to a diagram having
    ! N propagators (just by pinching N-4 propagators), in this case, the argument s_mat is the
    ! complete S matrix (NxN) and set_ref is the set of unpinched propagators, IT MUST HAVE
    ! RANK 1 AND SHAPE 4.
    !
    ! In the following the integrals f3p_x have to be called with argument
    ! s_mat_p. This was defined with the call to preparesmatrix.
    !
    verif1 = czero
    verif1(3) =  f4p_np2(s_mat_p,b_ref,b_null)
    res6 = verif1
    !
    verif2 = zero
    if ( (choix_kinem == 1) .or. (choix_kinem == 2) ) then
      !
      verif2(3:6) =  f4p1m("n+2",s_var,t_var,p4sq,0,0,0,0)
      !
    else if (choix_kinem == 3) then
      !
      verif2(3:6) =  f4p2m_opp("n+2",s_var,t_var,p2sq,p4sq,0,0,0,0)
      !
    else if (choix_kinem == 4) then
      !
      verif2(3:6) =  f4p2m_adj("n+2",s_var,t_var,p3sq,p4sq,0,0,0,0)
      !
    else if (choix_kinem == 5) then
      !
      verif2(3:6) =  f4p3m("n+2",s_var,t_var,p2sq,p3sq,p4sq,0,0,0,0)
      !
    else if (choix_kinem == 6) then
      !
      verif2(3:6) =  f4p4m("n+2",s_var,t_var,p1sq,p2sq,p3sq,p4sq,0,0,0,0)
      !
    end if
    !
  else if (choix == 6) then
    !
    ! Results for integrals in n+2 dimensions with two Feynman parameters
    ! in the numerator: z1*z2
    ! The same as the preceeding case. Note that in the case of the generic function
    ! the labels of the Feynman parameters can be put in any order, instead in the dedicated
    ! functions, they have to be ordered
    !
    verif1 = czero
    verif1(3) =  f4p_np2(s_mat_p,b_ref,b_null,2,1)
    res6 = verif1
    !
    verif2 = zero
    !
    if ( (choix_kinem == 1) .or. (choix_kinem == 2) ) then
      !
      verif2(3:6) =  f4p1m("n+2",s_var,t_var,p4sq,0,0,1,2)
      !
    else if (choix_kinem == 3) then
      !
      verif2(3:6) =  f4p2m_opp("n+2",s_var,t_var,p2sq,p4sq,0,0,1,2)
      !
    else if (choix_kinem == 4) then
      !
      verif2(3:6) =  f4p2m_adj("n+2",s_var,t_var,p3sq,p4sq,0,0,1,2)
      !
    else if (choix_kinem == 5) then
      !
      verif2(3:6) =  f4p3m("n+2",s_var,t_var,p2sq,p3sq,p4sq,0,0,1,2)
      !
    else if (choix_kinem == 6) then
      !
      verif2(3:6) =  f4p4m("n+2",s_var,t_var,p1sq,p2sq,p3sq,p4sq,0,0,1,2)
      !
    end if
    !
  else if (choix == 7) then
    !
    ! Results for integrals in n+4 dimension with zero Feynman parameters
    ! at the numerator:
    !
    verif1 = czero
    verif1(2:3) =  f4p_np4(s_mat_p,b_ref,b_null)
    res6 = verif1
    !
    verif2 = zero
    !
    if ( (choix_kinem == 1) .or. (choix_kinem == 2) ) then
      !
      verif2(3:6) =  f4p1m("n+4",s_var,t_var,p4sq,0,0,0,0)
      !
    else if (choix_kinem == 3) then
      !
      verif2(3:6) =  f4p2m_opp("n+4",s_var,t_var,p2sq,p4sq,0,0,0,0)
      !
    else if (choix_kinem == 4) then
      !
      verif2(3:6) =  f4p2m_adj("n+4",s_var,t_var,p3sq,p4sq,0,0,0,0)
      !
    else if (choix_kinem == 5) then
      !
      verif2(3:6) =  f4p3m("n+4",s_var,t_var,p2sq,p3sq,p4sq,0,0,0,0)
      !
    else if (choix_kinem == 6) then
      !
      verif2(3:6) =  f4p4m("n+4",s_var,t_var,p1sq,p2sq,p3sq,p4sq,0,0,0,0)
      !
    end if
    !
  else if (choix == 8) then
    !
    ! By default, mu2_scale_par = 1._ki
    !
    res6 =  c44(s_null)+b33(2,(/3/))
    !
    ! note that we have to reset the cache in order that the new value
    ! of mu2_scale_par will be effective.  A call to preparesmatrix
    ! is sufficient.
    !
    call preparesmatrix()
    !
    ! we change the value of mu^2
    !
    mu2 = 34._ki
    lmu2 = log(mu2/mu2_scale_par)
    mu2_scale_par = mu2
    res6a =  c44(s_null) +b33(2,(/3/))
  else if (choix == 9) then
    !
    ! Results for integrals in n dimension with five Feynman parameter
    ! in the numerator: z1*z2^2*z3^2
    !
    res6 = -a45(1,2,2,3,3,s_null)
    !
  else if (choix == 10) then
    !
    ! Results for integrals in n+2 dimension with three Feynman parameter
    ! in the numerator: z1*z2^2
    !
    res6 = 2._ki*b45(1,2,2,s_null)
    !
  else if (choix == 11) then
    !
    ! Results for integrals in n+4 dimensions with one Feynman parameter
    ! in the numerator: z1
    !
    res6 =  -4._ki*c45(1,s_null)
    !
  else if (choix == 12) then
    !
    ! Results for scalar integral in n+6 dimensions
    !
    res6 =  -8._ki*d46(s_null)
    !
  !
  end if
  !
  call cpu_time(ti2)
  !
  ! The results are written to the file with unit 17
  !
  write (17,*) 'The kinematics is:'
  write (17,*) ''
  write (17,*) '  p1              p4'
  write (17,*) '    \                /  '
  write (17,*) '     \     (4)     /   '
  write (17,*) '      |----<----|    '
  write (17,*) '      |            |    '
  write (17,*) ' (1) |            |(3) '
  write (17,*) '      |            |    '
  write (17,*) '      |---->----|    '
  write (17,*) '     /     (2)     \   '
  write (17,*) '    /                \  '
  write (17,*) '  p2              p3'
  write (17,*) ''
  write (17,*) 'p1+p2+p3+p4 = 0'
  write (17,*) ''
  write (17,*) '(p1+p2)^2 =',s_var
  write (17,*) '(p2+p3)^2 =',t_var
  write (17,*) '(p1)^2 =',p1sq
  write (17,*) '(p2)^2 =',p2sq
  write (17,*) '(p3)^2 =',p3sq
  write (17,*) '(p4)^2 =',p4sq
  write (17,*) 'm1^2 =',m1sq
  write (17,*) 'm2^2 =',m2sq
  write (17,*) 'm3^2 =',m3sq
  write (17,*) 'm4^2 =',m4sq
  write (17,*) '(mu)^2 =',mu2_scale_par
  write (17,*) ''
  write (17,*) 'defining I_N^n= mu^(4-n) \int d^n k/(i*Pi^(n/2))*func(k,p_i)'
  write (17,*) '= r_Gam *(P2/eps^2+P1/eps+P0),'
  write (17,*) 'n = 4-2*eps,'
  write (17,*) 'r_Gam = Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2eps)'
  write (17,*) 'the program gives numbers for P2,P1,P0'
  write (17,*) ''
  write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
  write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki),aimag(res6%b)
  write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%c,ki),aimag(res6%c)
  !
  write (6,*) ''
  write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
  write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki),aimag(res6%b)
  write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%c,ki),aimag(res6%c)
  write (17,*) ''
  !
  if ( ((choix == 5) .or. (choix == 6) .or. (choix == 7)) .and. (choix_kinem < 7) ) then
    !
    write (17,*) 'Check with dedicated function:'
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif2(1),verif2(2)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif2(3),verif2(4)
    write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif2(5),verif2(6)
    !
    !
    write (6,*) 'Check with dedicated function:'
    write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif2(1),verif2(2)
    write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif2(3),verif2(4)
    write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif2(5),verif2(6)
    !
  end if
  if (choix == 8) then
    !
    write (17,*) 'The preceding result has been computed with mu^2=',mu02
    write (17,*) ' '
    write (17,*) 'Now setting by hand mu^2=',mu2
    write (17,*) 'and expanding (',mu2,'/',mu02,')^epsilon around epsilon=0'
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki)+ &
         &real(res6%a,ki)*lmu2, aimag(res6%b)+aimag(res6%a)*lmu2
    write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') &
    &real(res6%c,ki) + lmu2*real(res6%b,ki) + lmu2**2*real(res6%a,ki)/2._ki,&
    &aimag(res6%c) + lmu2*aimag(res6%b) + lmu2**2*aimag(res6%a)/2._ki
    write (17,*) ''
    write (17,*) 'check with direct calculation using the global variable mu2_scale_par=',mu2
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a%a,ki),aimag(res6a%a)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a%b,ki),aimag(res6a%b)
    write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a%c,ki),aimag(res6a%c)
    ! ***********************************
    write (6,*) ' '
    write (6,*) 'The preceding result has been computed with mu^2=',mu02
    write (6,*) 'Now setting by hand mu^2=',mu2
    write (6,*) 'and expanding (',mu2,'/',mu02,')^epsilon around epsilon=0'
    write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
    write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki)+ &
         & real(res6%a,ki)*lmu2, aimag(res6%b)+aimag(res6%a)*lmu2
    write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') &
    &real(res6%c,ki) + lmu2*real(res6%b,ki) + lmu2**2*real(res6%a,ki)/2._ki,&
    &aimag(res6%c) + lmu2*aimag(res6%b) + lmu2**2*aimag(res6%a)/2._ki
    !
    write (6,*) ''
    write (6,*) 'check with direct calculation using the global variable mu2_scale_par=',mu2
    write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a%a,ki),aimag(res6a%a)
    write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a%b,ki),aimag(res6a%b)
    write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6a%c,ki),aimag(res6a%c)
    !
    write (6,*) ''
  end if
  !
  write (17,*) 'CPU time=',ti2-ti1
  !
  ! routine to free the cache and allocated memory
  !
  call exitgolem95()
  !
  close(17)
  !close(19)
  !
end program main
