!
! this program computes three-point functions
! in n and n+2 dimensions.
! One can use the form factors
! or directly the dedicated three-point functions.
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
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type
  use form_factor_3p ! module containing the three-point form factors
  use constante
  use generic_function_3p ! module containing the generic three-point function
  use parametre
  use generic_function_np ! module with the generic three-point higher-rank functions
  use form_factor_higher_ranks
  !
  implicit none
  !
  type(form_factor) :: res6,res6a
  complex(ki), dimension(3) :: verif
  real(ki) :: mass_sq_1,mass_sq_2,mass_sq_3
  !real(ki) :: mass_int_sq_1,mass_int_sq_2,mass_int_sq_3,mu2,mu02,lmu2
  complex(ki) :: mass_int_sq_1,mass_int_sq_2,mass_int_sq_3
  real(ki) :: mu2,mu02,lmu2
  integer :: choix,choix_kinem
  integer, dimension(1) :: array
  !
  write (*,*) 'Choose from the following kinematics:'
  write (*,*) '1) one off-shell leg, no internal masses'
  write (*,*) '2) two off-shell legs, no internal masses'
  write (*,*) '3) three off-shell legs, no internal masses (finite)'
  write (*,*) '4) two off-shell legs, one internal mass (QL 3), '
  write (*,*) '5) one off-shell leg, one on-shell massive leg (one internal mass, QL 4)'
  write (*,*) '6) two on-shell massive legs (one internal mass, QL 5)'
  write (*,*) '7) one off-shell leg, two on-shell massive legs (two internal masses, QL 6)'
  read (*,*) choix_kinem
  !
  ! Opening of the error files
  !
  open(unit=19,file='error_3point.txt',status='unknown')
  !
  ! opening of the files containing the results
  !
  open(unit=17,file='test3point.txt',status='unknown')
  !
  ! These are the entries of the S matrix
  ! They are related to the cuts of the following diagram
  ! All the momenta are incoming : p1+p2+p3 = 0
  !
  !            |
  !            | p1
  !            |
  !           /\
  !          /  \
  !     (1) /    \ (3)
  !        /      \
  !       /--->----\
  !      /   (2)    \
  ! p2  /            \ p3
  !
  ! S(1,2) = p2^2
  ! S(2,3) = p3^2
  ! S(3,1) = p1^2
  !
  ! Allocates memory to store the set of initial propagators, the S matrix,
  ! its inverse and the b coefficients.
  ! This call will allocate a derived type s_mat_p object.
  ! Includes calls to allocation_s and initializes the caching system
  !
  call initgolem95(3)
  !
  !
  if (choix_kinem == 1) then
    !
    mass_int_sq_1 = 0.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = 0.0_ki
    mass_sq_2 = 0.0_ki
    mass_sq_3 = -2.0_ki
    !
  else if (choix_kinem == 2) then
    !
    mass_int_sq_1 = 0.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = 0.0_ki
    mass_sq_2 = 10.0_ki
    mass_sq_3 = -60.0_ki
    !
  else if (choix_kinem == 3) then
    !
    mass_int_sq_1 = 0.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = -50.0_ki
    mass_sq_2 = 10.0_ki
    mass_sq_3 = -60.0_ki
    !
  ! case p1^2, p2^2 /= m1^2, internal line 1 massive (QL3)
  else if (choix_kinem == 4) then
    !
    mass_int_sq_1 = 20.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = -123.0_ki
    mass_sq_2 = -60.0_ki
    mass_sq_3 = 0.0_ki
    !
  ! case p1^2 = m1^2, p2^2 /= m1^2, internal line 1 massive (QL4)
  else if (choix_kinem == 5) then
    !
    mass_int_sq_1 = 5.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = mass_int_sq_1
    mass_sq_2 = 7.0_ki
    mass_sq_3 = 0.0_ki
    !
  ! case p1^2 = m1^2, p2^2 = m1^2, internal line 1 massive (QL5)
  else if (choix_kinem == 6) then
    !
    mass_int_sq_1 = 5.0_ki
    mass_int_sq_2 = 0._ki
    mass_int_sq_3 = 0.0_ki
    mass_sq_1 = mass_int_sq_1
    mass_sq_2 = mass_int_sq_1
    mass_sq_3 = 0.0_ki
    !
  ! case p2^2 = m1^2; p1^2 /= m1^2,m3^2,nonzero; p3^2=m3^2, internal lines 1,2 massive (QL6)
  else if (choix_kinem == 7) then
    !
    mass_int_sq_1 = 9.0_ki
    mass_int_sq_2 = 0.0_ki
    mass_int_sq_3 = 3.0_ki
    mass_sq_1 = 25.0_ki
    mass_sq_2 = mass_int_sq_1
    mass_sq_3 = mass_int_sq_3
    !
  ! case p2^2=m1^2; p1^2 /= m1^2,nonzero; p3^2=m1^2, internal lines 1,2 massive (QL6 with internal masses equal)
  else if (choix_kinem == 8) then
    !
    mass_int_sq_1 = 5.0_ki
    mass_int_sq_2 = 0.0_ki
    mass_int_sq_3 = mass_int_sq_1
    mass_sq_1 = 25.0_ki
    mass_sq_2 = mass_int_sq_1
    mass_sq_3 = mass_int_sq_3
    !
  else if (choix_kinem == 9) then
    !
    mass_int_sq_1 = (7.6751408576965332_ki,0.0000000000000000_ki)
    mass_int_sq_2 = (7.6751408576965332_ki,-0.14250832796096802_ki)
    mass_int_sq_3 = (0.0000000000000000_ki,0.0000000000000000_ki)
    mass_sq_1 = 4.3248424530029297_ki
    mass_sq_2 = 0._ki
    mass_sq_2 = 1.e-9_ki
    mass_sq_3 = -1.7544794082641602_ki
    !
  end if
  !
  ! Definition of the S matrix
  !
  s_mat(1,1) = -2.0_ki*mass_int_sq_1
  s_mat(1,2) = mass_sq_2 - mass_int_sq_1 - mass_int_sq_2
  s_mat(1,3) = mass_sq_1 - mass_int_sq_1 - mass_int_sq_3
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = -2.0_ki*mass_int_sq_2
  s_mat(2,3) = mass_sq_3 - mass_int_sq_2 - mass_int_sq_3
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = -2.0_ki*mass_int_sq_3
  !
  ! This call fills the internal array s_mat_r.
  ! It also assigns the integers in s_mat_p, which encode the positions
  ! of complex mass entries and zero mass entries. It includes call to init_invs
  !
  call preparesmatrix()
  !
  !
  write (*,*) 'Choose what the program should compute:'
  write (*,*) '0) scalar three-point function in n dimensions'
  write (*,*) '1) three-point function in n dimensions with one Feynman parameter'
  write (*,*) '2) three-point function in n dimensions with two Feynman parameters'
  write (*,*) '3) three-point function in n dimensions with three Feynman parameters'
  write (*,*) '4) scalar three-point function in n+2 dimensions'
  write (*,*) '5) three-point function in n+2 dimensions with one Feynman parameter'
  write (*,*) '6) test of the mu independence'
  write (*,*) '7) three-point function in n dimensions with four Feynman parameters'
  write (*,*) '8) three-point function in n+2 dimensions with two Feynman parameters'
  write (*,*) '9) scalar three-point function in n+4 dimensions'
  read (*,*) choix
  !
  ! info for user
  if (choix == 0) then
   !
    write (*,*) 'calculating n-dim scalar 3-point fctn. with '
    !
  else if   (choix == 1) then
    !
    write (*,*) 'calculating n-dim rank one (z1) 3-point fctn. with'
    !
  else if   (choix == 2) then
    !
    write (*,*) 'calculating n-dim rank two (z1*z2) 3-point fctn. with'
    !
  else if   (choix == 3) then
    !
    write (*,*) 'calculating n-dim rank three (z1^2*z3) 3-point fctn. with'
    !
  else if   (choix == 4) then
    !
    write (*,*) 'calculating (n+2)-dim scalar 3-point fctn. with'
    !
 else if (choix == 5) then
    !
    write (*,*) 'calculating (n+2)-dim rank one (z2) 3-point fctn. with'
    !
  end if
  !
  if (choix_kinem == 1) then
    !
    write (*,*) 'one off-shell leg'
    !
  else if (choix_kinem == 2) then
    !
    write (*,*) 'two off-shell legs'
    !
  else if (choix_kinem == 3) then
    !
    write (*,*) 'three off-shell legs'
    !
  end if
  !
  write (*,*) 'the result has been written to the file test3point.txt'
  !
  ! start calculation
  !
  ! To change the value of mu^2 (in GeV) (set to 1. by default)
  ! uncomment this line
  ! mu2_scale_par = 12._ki
  !
  ! store original mu^2
  mu02 = mu2_scale_par
  !
  !
  ! In the following the integrals f3p_x have to be called with argument
  ! s_mat_p. This was defined with the call to preparesmatrix.
  !
  if (choix == 0) then
    !
    ! Result for the scalar integral in n dimension
    !
    ! both are working
    !
    res6 = a30(s_null)
    verif = f3p(s_mat_p,b_ref)
    !
    ! the labels 1,2,3 correspond to Feynman parameters z1,z2,z3
  else if (choix == 1) then
    !
    ! Results for integrals in n dimensions with one Feynman parameter
    ! in the numerator: z1
    !
    res6 = -a31(2,s_null)
    verif = f3p(s_mat_p,b_ref,2)
    !~ res6 = -a31(3,s_null)
    !~ verif = f3p(s_mat,b_ref,3)
    !
  else if (choix == 2) then
    !
    ! Results for integrals in n dimensions with two Feynman parameters
    ! in the numerator: z1*z2
    !
    res6 =  a32(1,2,s_null)
    verif = f3p(s_mat_p,b_ref,1,2)
    !~ res6 =  a32(2,2,s_null)
    !~ verif = f3p(s_mat,b_ref,2,2)
    !
  else if (choix == 3) then
    !
    ! Results for integrals in n dimensions with three Feynman parameters
    ! at the numerator: z1^2*z3
    !
    res6 =  -a33(2,2,2,s_null)
    verif = f3p(s_mat_p,b_ref,2,2,2)
    !
  else if (choix == 4) then
    !
    ! Results for integrals in n+2 dimensions with no Feynman parameters
    ! at the numerator:
    !
    res6 =  -2.0_ki*b32(s_null)
    verif = czero  ! complex zero. defined in module constante.
    verif(2:3) = f3p_np2(s_mat_p,b_ref)
    !
  else if (choix == 5) then
    !
    ! Results for integrals in n+2 dimensions with one Feynman parameters
    ! at the numerator: z2
    !
    res6 =  2.0_ki*b33(2,s_null)
    verif = czero
    verif(2:3) = f3p_np2(s_mat_p,b_ref,2)
    !
  else if (choix == 6) then
    !
    ! by default, mu2_scale_par = 1._ki
    ! take scalar triangle as example
    res6=A30(s_null)
    ! we have to reset the cache in order that the new value
    ! of mu2_scale_par will be effective. A call to preparesmatrix
    ! is sufficient.
    !
    call preparesmatrix()
    !
    ! we change the value of mu^2
    !
    mu2 = 34._ki
    lmu2 = log(mu2/mu2_scale_par)
    mu2_scale_par = mu2
    res6a =  A30(s_null)
  else if (choix == 7) then
    !
    ! Results for integrals in n dimensions with four Feynman parameters
    ! in the numerator: z1^2*z2*z3
    !
    res6 =  a34(1,1,2,3,s_null)
    res6a = fnp_generic(3,0,0,4,(/1,1,2,3/))
    verif = (/ res6a%a , res6a%b , res6a%c /)
  else if (choix == 8) then
    !
    ! Results for integrals in n+2 dimensions with two Feynman parameters
    ! in the numerator: z1*z2
    !
    res6 =  -2._ki*b34(1,2,s_null)
    res6a = fnp_generic(3,2,0,2,(/1,2/))
    verif = (/ res6a%a , res6a%b , res6a%c /)
  else if (choix == 9) then
    !
    ! Results for scalar integrals in n+4 dimensions
    !
    res6 =  4._ki*c34(s_null)
    res6a = fnp_generic(3,4,0,0, (/ integer :: /) )
    verif = (/ res6a%a , res6a%b , res6a%c /)

 !
  end if
  !
  write (17,*) 'The kinematics is:'
  write (17,*) ''
  write (17,*) '             |             '
  write (17,*) '             | p1        '
  write (17,*) '             |             '
  write (17,*) '            /\            '
  write (17,*) '           /  \           '
  write (17,*) '     (1) /     \ (3)     '
  write (17,*) '        /        \         '
  write (17,*) '       /--->----\       '
  write (17,*) '      /   (2)      \      '
  write (17,*) 'p2 /               \ p3'
  write (17,*) ''
  write (17,*) 'p1+p2+p3 = 0'
  write (17,*) ''
  write (17,*) '(p1)^2 =',mass_sq_1
  write (17,*) '(p2)^2 =',mass_sq_2
  write (17,*) '(p3)^2 =',mass_sq_3
  write (17,*) 'm1^2 =',mass_int_sq_1
  write (17,*) 'm2^2 =',mass_int_sq_2
  write (17,*) 'm3^2 =',mass_int_sq_3
  write (17,*) 'mu^2 =',mu2_scale_par
  write (17,*) ''
  write (17,*) 'defining I_N^n= mu^(4-n) \int d^n k/(i*Pi^(n/2))*func(k,p_i)'
  write (17,*) '= r_Gam *(P2/eps^2+P1/eps+P0),'
  write (17,*) 'n = 4-2*eps,'
  write (17,*) 'r_Gam = Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2eps)'
  write (17,*) 'the program gives numbers for P2,P1,P0'
  write (17,*) ''
    !
  write (17,*) 'result='
  write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
  write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki),aimag(res6%b)
  write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%c,ki),aimag(res6%c)
  write (17,*) ''
  write (6,*) 'result='
  write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
  write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki),aimag(res6%b)
  write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%c,ki),aimag(res6%c)
  !
  if ( choix < 6 .or. choix>=7) then
  !
  write (6,*) 'Check with dedicated function:'
  write (6,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif(1)
  write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif(2)
  write (6,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif(3)
  write (17,*) 'Check with dedicated function:'
  write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif(1)
  write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif(2)
  write (17,'("+ 1           * (",e16.10,1x,"+ I*",1x,e16.10,")")') verif(3)
  !
  else if ( choix ==6 ) then
  !
    write (17,*) 'The preceding result has been computed with mu^2=',mu02
    write (17,*) ' '
    write (17,*) 'Now setting by hand mu^2=',mu2
    write (17,*) 'and expanding (',mu2,'/',mu02,')^epsilon around epsilon=0'
    write (17,'("  1/epsilon^2 * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%a,ki),aimag(res6%a)
    write (17,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki)+ &
         &real(res6%a,ki)*lmu2,aimag(res6%b)+aimag(res6%a)*lmu2
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
    write (6,'("+ 1/epsilon   * (",e16.10,1x,"+ I*",1x,e16.10,")")') real(res6%b,ki) + &
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
    write (6,*) ' '
  !
  else
  !
  write (6,*) 'invalid choice, option number must be <= 9'
  !
  endif
  ! routine to free the cache and allocated memory
  !
  call exitgolem95()
  !
  close(17)
  close(19)
  !
end program main
!
!
