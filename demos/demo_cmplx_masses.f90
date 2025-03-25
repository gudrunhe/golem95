! 
! computes form factors with complex masses
!
!
program main
  !
  use precision_golem ! to get the type ki (for real and complex)
  use matrice_s
  use form_factor_type, only: form_factor, operator(-), operator(*), operator(+), assignment(=)
  use form_factor_3p ! module containing the three-point form factors (export all)
  use form_factor_1p
  use form_factor_2p
  use form_factor_6p
  use form_factor_5p
  use form_factor_4p
  use spinor, only: scalar
  use parametre 
  !
  implicit none
  !
  type(form_factor) :: res
  !
  complex(ki), dimension(6) :: mass_int_sq
  real(ki) :: t1, t2, s_in
  real(ki) :: a1, a2, a3, m1, m2, m3
  complex(ki) :: a1c, a2c, a3c, cz, m1c, m2c,m3c
  real(ki), dimension(4) :: p1,p2,p3,p4,p5,p6
  real(ki), dimension(4) :: p12,p23,p34,p45,p56,p61,p123,p234,p345
  !
  integer :: b_pin, i, i1, i2, i3, i4, i5 ,i6
  real(ki), dimension(:), allocatable :: set_s, set_m1, set_m2, set_m3, set_arg1, set_arg2, set_arg3
  complex(ki), dimension(:), allocatable :: set_m1c, set_m2c, set_m3c, set_arg1c, set_arg2c, set_arg3c
  logical :: six_test, at, cc
  !
  character(len=*), parameter :: formatff = '(a12,":   ","(",2ES22.15") , (",2ES22.15") , (",2ES22.15,")")'
  character(len=*), parameter :: formatff5 = '(a3,"-",5i0," :   ","(",2ES22.15") , (",2ES22.15") , (",2ES22.15,")")'
  character(len=*), parameter :: formatff6 = '(a3,"-",6i0," :   ","(",2ES22.15") , (",2ES22.15") , (",2ES22.15,")")'
  !
!  if_print_warn_par = .true.
  !
  call cpu_time(t1)
  !
!!! if only real masses are given, the program will either run 6pt-form factors (six_test) or mainly ff1 to ff5.
!!! results are written to test_ff6_r.txt or test_ff_r.txt or test_ff_r_at.txt
!!! at triggers 'all triangles' which calculates also triangles that cannot have complex entries.
!!!
!!! if masses have an imaginary part, the program writes to test_ff(6)_c.txt.
  !
  six_test = .true.
  cc = .true.
  at = .false.
  !
!!! Allocates memory to store the set of initial propagators, the S matrix, 
!!! its inverse and the b coefficients.
!!! This call will allocate a derived type s_mat_p object. It includes associated pointers
!!! to s_mat_r (case rmass and cmass) and s_mat_c (cmass only)
!!! includes call to allocation_s and initializes the caching system
  !
  call initgolem95(6)
  !
  p1 = (/ 0.5_ki , 0._ki , 0._ki , 0.5_ki /)
  p2 = (/ 0.5_ki , 0._ki , 0._ki , -0.5_ki /)
  p3 = -(/ 0.19178191094778038_ki , 0.12741179719516801_ki , 0.08262476614744381_ki , 0.11713105190921771_ki /)
  p4 = -(/ 0.33662712284553753_ki , -0.06648281097623857_ki , -0.3189378514746887_ki , -0.08471424069583446_ki /)
  p5 = -(/ 0.21604814388379073_ki , -0.20363139428835617_ki , 0.044157623555325_ki , 0.0571065672034082_ki /)
  p6 = -(/ 0.2555428223228916_ki , 0.1427024080694266_ki , 0.19215546177191994_ki , -0.08952337841679145_ki /)
  !
  !  p3 = -(/ 0.24721209116616516_ki, 0.11438327254920486_ki, 0.06268384608884053_ki, 0.10558040234674371_ki /)
  !  p4 = -(/ 0.07213339280895899_ki, 0.06871472648494759_ki, -0.0017424243541284828_ki, 0.021228339189685863_ki /)
  !  p5 = -(/ 0.1418614527530244_ki, 0.1103700408743443_ki, 0.0107542809493891_ki, 0.0884729975521105_ki /)
  !  p6 = -(/ 0.5387930632718515_ki, -0.29346803990849674_ki, -0.07169570268410115_ki, -0.21528173908854006_ki /)
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
  !
!!! Definition of the S matrix
!!! Fill out only s_mat or s_mat_c
  ! 
  !
  s_mat(1,1) = 0._ki
  s_mat(1,2) = scalar(p2,p2)
  !  s_mat(1,2) = scalar(p12,p12)
  s_mat(1,3) = scalar(p23,p23)
  s_mat(1,4) = scalar(p234,p234)
  s_mat(1,5) = scalar(p61,p61)
  s_mat(1,6) = scalar(p1,p1)
  !
  s_mat(2,1) = s_mat(1,2)
  s_mat(2,2) = 0._ki
  s_mat(2,3) = scalar(p3,p3)
  s_mat(2,4) = scalar(p34,p34)
  s_mat(2,5) = scalar(p345,p345)
  s_mat(2,6) = scalar(p12,p12)
  !
  s_mat(3,1) = s_mat(1,3)
  s_mat(3,2) = s_mat(2,3)
  s_mat(3,3) = 0._ki
  s_mat(3,4) = scalar(p4,p4)
  s_mat(3,5) = scalar(p45,p45)
  s_mat(3,6) = scalar(p123,p123)
  !
  s_mat(4,1) = s_mat(1,4)
  s_mat(4,2) = s_mat(2,4)
  s_mat(4,3) = s_mat(3,4)
  s_mat(4,4) = 0._ki
  s_mat(4,5) = scalar(p5,p5)
  s_mat(4,6) = scalar(p56,p56)
  !
  s_mat(5,1) = s_mat(1,5)
  s_mat(5,2) = s_mat(2,5)
  s_mat(5,3) = s_mat(3,5)
  s_mat(5,4) = s_mat(4,5)
  s_mat(5,5) = 0._ki
  s_mat(5,6) = scalar(p6,p6)
  !
  s_mat(6,1) = s_mat(1,6)
  s_mat(6,2) = s_mat(2,6)
  s_mat(6,3) = s_mat(3,6)
  s_mat(6,4) = s_mat(4,6)
  s_mat(6,5) = s_mat(5,6)
  s_mat(6,6) = 0._ki
  !
  !
  ! Definition of internal masses.
  !
  mass_int_sq(1) = cmplx(.037_ki,.0_ki,ki)
  mass_int_sq(2) = cmplx(.0_ki,.0_ki,ki)
  mass_int_sq(3) = cmplx(.015_ki,.0_ki,ki)
  mass_int_sq(4) = cmplx(.05_ki,.0_ki,ki)
  mass_int_sq(5) = cmplx(.0_ki,.0_ki,ki)
  mass_int_sq(6) = cmplx(.072_ki,.0_ki,ki)
  !
  if (cc) then
     !
     mass_int_sq(6) = cmplx(.072_ki,-.0000000001_ki,ki)
!     mass_int_sq(3) = cmplx(.015_ki,-.004101_ki,ki) 
!     mass_int_sq(1) = cmplx(.037_ki,-.0180000000001_ki,ki)
     !
  end if
  !
!  mass_int_sq(:) = cmplx(0._ki,0._ki,ki)
  !
  !  Finalize entries in s matrix.
  !
  do i=1,6 
     !
     s_mat(i,:) = s_mat(i,:) - mass_int_sq(:) - mass_int_sq(i)
     !
  end do
  !
!!! This call fills in s_mat_r.
!!! It also assigns the bit integers in s_mat_p, which describe the positions
!!! of complex mass entries and zero mass entries. It includes call to init_invs
  !
  call preparesmatrix()
  !
  call cpu_time(t2)
  !
  write (*,*) 'time for set-up: ', t2-t1
  !
  write (*,*) '*********************************'
  write (*,*) 'form factors'
  !
  if (six_test) then
     !
     if (cc) then
        !
        open(unit=17,file='test_ff6_c.txt',status='unknown')
        !
     else
        !
        open(unit=17,file='test_ff6_r.txt',status='unknown')
        !
     end if
     !
     write (*,*) '*********************************'
     write (*,*) ' '
     write (*,*) '6pt - function'
     write (17,*) '*********************************'
     write (17,*) ' '
     write (17,*) '6pt - function'
     !
     call cpu_time(t1)
     !
     i1=0
     i2=0
     i3=0
     i4=0
     i5=0
     i6=0
     i=0
     !
     b_pin = 0
     !
     res = a60(b_pin)
     !
     write (17,formatff6) 'a60',0,0,0,0,0,0, res
     !
     do while (i1 .lt. 6)
        !
        i1 = i1 + 1
        if (i1 .eq. i) i1 = i1 + 1
        write(*,*)'i1',i1   
        res = a61(i1,b_pin)
        write (17,formatff6) 'a61',0,0,0,0,0,i1, res

        i2 = i1 - 1

        do while (i2 .lt. 6)

           i2 = i2 + 1
           if (i2 .eq. i) i2 = i2 + 1
           !write(*,*)'i2',i2
           res = a62(i1,i2,b_pin)
           write (17,formatff6) 'a62',0,0,0,0,i1,i2, res

           i3 = i2 - 1

           do while (i3 .lt. 6)

              i3 = i3 + 1
              if (i3 .eq. i) i3 = i3 + 1
              !    write(*,*)'i3',i3     
              res = a63(i1,i2,i3,b_pin)
              write (17,formatff6) 'a63',0,0,0,i1,i2,i3, res

              i4 = i3 - 1

              do while (i4 .lt. 6)

                 i4 = i4 + 1
                 if (i4 .eq. i) i4 = i4 + 1
                 !            write(*,*)'i4',i4
                 res = a64(i1,i2,i3,i4,b_pin)
                 write (17,formatff6) 'a64',0,0,i1,i2,i3,i4, res
                 i5 = i4 - 1

                 do while (i5 .lt. 6)
                    i5 = i5 + 1
                    if (i5 .eq. i) i5 = i5 + 1
                    ! write(*,*)'i5',i5
                    res = a65(i1,i2,i3,i4,i5,b_pin)
                    write (17,formatff6) 'a65',0,i1,i2,i3,i4,i5, res
                    i6 = i5 - 1
                    do while (i6 .lt. 6)
                       i6 = i6 + 1
                       if (i6 .eq. i) i6 = i6 + 1
                       ! write(*,*)'i5',i5
                       res = a66(i1,i2,i3,i4,i5,i6,b_pin)
                       write (17,formatff6) 'a66',i1,i2,i3,i4,i5,i6, res

                    end do

                 end do

              end do

           end do

        end do

     end do
     !
     call cpu_time(t2)
     !
     write(*,*) "time-6pt-gfortran ", t2-t1
     write(17,*) "time-6pt-gfortran ", t2-t1
     !
  else
     !
     if (.not. cc) then
        !
        if (at) then
           !
           open(unit=17,file='test_ff_r_at.txt',status='unknown')
           !
        else
           !
           open(unit=17,file='test_ff_r.txt',status='unknown')
           !
        end if
        !
     else
        !
        open(unit=17,file='test_ff_c.txt',status='unknown')
        !
     end if
     !
     !
     write (*,*) '*********************************'
     write (*,*) ' '
     write (*,*) '6pt - function'
     write (17,*) '*********************************'
     write (17,*) ' '
     write (17,*) '6pt - function'
     !
     ! just a few 6pts
     !
     b_pin = 0
     !
     res = a60(b_pin)
     write (17,formatff) 'a60', res
     !
     res = a61(2,b_pin)
     write (17,formatff) 'a61-2', res
     !
     res = a62(2,3,b_pin)
     write (17,formatff) 'a62-23', res
     !
     res = a63(2,3,5,b_pin)
     write (17,formatff) 'a63-235', res
     !
     res = a64(1,1,2,3,b_pin)
     write (17,formatff) 'a64-1123 ', res
     !
     res = a65(1,3,4,5,6,b_pin)
     write (17,formatff) 'a65-13456', res
     !
     res = a66(1,2,3,4,5,6,b_pin)
     write (17,formatff) 'a66-123456', res
     !
     write (*,*) '*********************************'
     write (*,*) ' '
     write (*,*) '5pt - function'
     write (17,*) '*********************************'
     write (17,*) ' '
     write (17,*) '5pt - function'
     !
     !
     b_pin = 32
     !
     i = 5
     i1 = 0
     i2 = 0
     i3 = 0
     i4 = 0
     i5 = 0
     !
     res = a50(b_pin)
     write (17,formatff5) 'a50',0,0,0,0,0, res
     res = b52(b_pin)
     write (17,formatff5) 'b52',0,0,0,0,0, res
     res = c54(b_pin)
     write (17,formatff5) 'c54',0,0,0,0,0, res
     !
     do while (i1 .lt. 6)
        !
        i1 = i1 + 1
        if (i1 .eq. i) i1 = i1 + 1
        !write(*,*)'i1',i1   
        res = a51(i1,b_pin)
        write (17,formatff5) 'a51',0,0,0,0,i1, res
        res = b53(i1,b_pin)
        write (17,formatff5) 'b53',0,0,0,0,i1, res
        res = c55(i1,b_pin)
        write (17,formatff5) 'c55',0,0,0,0,i1, res

        i2 = i1 - 1

        do while (i2 .lt. 6)

           i2 = i2 + 1
           if (i2 .eq. i) i2 = i2 + 1
           !write(*,*)'i2',i2
           res = a52(i1,i2,b_pin)
           write (17,formatff5) 'a52',0,0,0,i1,i2, res
           res = b54(i1,i2,b_pin)
           write (17,formatff5) 'b54',0,0,0,i1,i2, res
           i3 = i2 - 1

           do while (i3 .lt. 6)

              i3 = i3 + 1
              if (i3 .eq. i) i3 = i3 + 1
              !    write(*,*)'i3',i3     
              res = a53(i1,i2,i3,b_pin)
              write (17,formatff5) 'a53',0,0,i1,i2,i3, res
              res = b55(i1,i2,i3,b_pin)
              write (17,formatff5) 'b55', 0,0,i1,i2,i3, res
              i4 = i3 - 1

              do while (i4 .lt. 6)

                 i4 = i4 + 1
                 if (i4 .eq. i) i4 = i4 + 1
                 !            write(*,*)'i4',i4
                 res = a54(i1,i2,i3,i4,b_pin)
                 write (17,formatff5) 'a54',0,i1,i2,i3,i4, res
                 i5 = i4 - 1

                 do while (i5 .lt. 6)
                    i5 = i5 + 1
                    if (i5 .eq. i) i5 = i5 + 1
                    ! write(*,*)'i5',i5
                    res = a55(i1,i2,i3,i4,i5,b_pin)
                    write (17,formatff5) 'a55',i1,i2,i3,i4,i5, res

                 end do

              end do

           end do

        end do

     end do
     !
     write (*,*) '*********************************'
     write (*,*) ' '
     write (*,*) '1pt - function'
     write (17,*) '*********************************'
     write (17,*) ' '
     write (17,*) '1pt - function'
     !
     res = a10(62)
     write (17,formatff) 'a10(62)',res
     !
     write (*,*) '*********************************'
     write (*,*) ' '
     write (*,*) '2pt - function'
     write (17,*) '*********************************'
     write (17,*) ' '
     write (17,*) '2pt - function'
     !
     b_pin = 120
     cz = cmplx(0._ki,0._ki,ki)
     m1c = mass_int_sq(6)
     m2c = mass_int_sq(3)
     !
     ! We change the values of the s matrix frequently now.
     ! Each time preparesmatrix should be called and
     ! the cache should be reset. 
     !
     s_mat(1,1) = -2._ki*m1c
     s_mat(2,2) = -2._ki*m2c
     !
     !call ffini
     ! write (*,*)'lt b0:', B0i(1,real(.234,ki_avh),real(m1,ki_avh),real(m2,ki_avh))
     !
     !
     allocate(set_s(8))
     set_s = (/ -2.45*m1, 0._ki, (sqrt(m1)-sqrt(m2))**2, m1, m2, 1.0003_ki*m1, &
          &     (sqrt(m1)+sqrt(m2))**2, 3.12_ki*(sqrt(m1)+sqrt(m2))**2  /)
     !
     do i=1,8
        !
        s_in = set_s(i)
        s_mat(1,2) = cmplx(s_in,0._ki,ki)-m1c-m2c
        s_mat(2,1) = s_mat(1,2)
        !
        ! We still use a submatrix of the six-dimensional s_matrix....
        !
        call preparesmatrix()
        !
        !res = a20(b_pin)
        !temp0 = cmplx(B0i(1,real(s_in,ki_avh),real(m1,ki_avh),real(m2,ki_avh)),kind=ki)
        !write (17,formatff) 'a20',res
        !write (*,*) i,' : ', (res%c-temp0)/abs(temp0)
        !
        res = a20(b_pin)
        write (17,formatff) 'a20',res
        res = a21(1,b_pin)
        write (17,formatff) 'a21-1',res
        res = a21(2,b_pin)
        write (17,formatff) 'a21-2',res
        res = a22(1,1,b_pin)
        write (17,formatff) 'a22-11',res
        res = a22(1,2,b_pin)
        write (17,formatff) 'a22-12',res
        res = a22(2,2,b_pin)
        write (17,formatff) 'a22',res
        
        res = b22(b_pin)
        write (17,formatff) 'b22',res
        
        write (17,*) '******'
        
     end do
     !
     deallocate(set_s)
     !
     m2 = 0._ki
     m2c = cz
     !
     s_mat(2,2) = cz
     !
     allocate(set_s(5))
     set_s = (/ -2.45_ki*m1, 0._ki,.672_ki*m1, m1, 1.2003_ki*m1 /)
     !
     do i=1,5
        !
        s_in = set_s(i)
        !
        s_mat(1,2) = cmplx(s_in,0._ki,ki)-m1c-m2c
        s_mat(2,1) = s_mat(1,2)
        !
        call preparesmatrix()
        !
        res = a20(b_pin)
        write (17,formatff) 'a20',res
        res = a21(1,b_pin)
        write (17,formatff) 'a21-1',res
        res = a21(2,b_pin)
        write (17,formatff) 'a21-2',res
        res = a22(1,1,b_pin)
        write (17,formatff) 'a22-11',res
        res = a22(1,2,b_pin)
        write (17,formatff) 'a22-12',res
        res = a22(2,2,b_pin)
        write (17,formatff) 'a22',res
        !
        res = b22(b_pin)
        write (17,formatff) 'b22',res
        !
        write (17,*) '******'
        !
     end do
     !
     deallocate(set_s)
     m1 = 0._ki
     m1c = cz
     !
     s_mat(1,1) = cz
     !
     allocate(set_s(2))
     set_s = (/ -.175_ki,.175_ki /)
     !
     do i=1,2
        !
        s_in = set_s(i)
        !
        s_mat(1,2) = cmplx(s_in,0._ki,ki)-m1c-m2c
        s_mat(2,1) = s_mat(1,2)
        !
        call preparesmatrix()
        !
        res = a20(b_pin)
        write (17,formatff) 'a20',res
        res = a21(1,b_pin)
        write (17,formatff) 'a21-1',res
        res = a21(2,b_pin)
        write (17,formatff) 'a21-2',res
        res = a22(1,1,b_pin)
        write (17,formatff) 'a22-11',res
        res = a22(1,2,b_pin)
        write (17,formatff) 'a22-12',res
        res = a22(2,2,b_pin)
        write (17,formatff) 'a22',res
        !
        res = b22(b_pin)
        write (17,formatff) 'b22',res
        !
        write (17,*) '******'
        !
     end do
     !
     deallocate(set_s)
     !
     write (*,*) '*********************************'
     write (*,*) ' '
     write (*,*) '3pt - function'
     write (17,*) '*********************************'
     write (17,*) ' '
     write (17,*) '3pt - function'
     !
     b_pin = 112
     !
     allocate(set_arg1(17))
     allocate(set_arg2(17))
     allocate(set_arg3(17))
     allocate(set_m1(17))
     allocate(set_m2(17))
     allocate(set_m3(17))
     allocate(set_arg1c(5))
     allocate(set_arg2c(5))
     allocate(set_arg3c(5))
     allocate(set_m1c(5))
     allocate(set_m2c(5))
     allocate(set_m3c(5))
     !
     m1c = cmplx(real(mass_int_sq(1),ki),0._ki,ki)
     m2c = cmplx(real(mass_int_sq(3),ki),0._ki,ki)
     m3c = mass_int_sq(6)
     !
     m1 = real(m1c,ki)
     m2 = real(m2c,ki)
     m3 = real(m3c,ki)
     !
     a1c = cmplx(.238751_ki,0._ki,ki)! - m1c - m2c
     a2c = cmplx(1.23062_ki,0._ki,ki)! - m2c - m3c
     a3c = cmplx(.1930573_ki,0._ki,ki)! - m1c - m3c
     !
     a1 = real(a1c,ki)
     a2 = real(a2c,ki)
     a3 = real(a3c,ki)
     !
     !
     if (at) then

        set_m1 = (/ 0._ki,0._ki,0._ki,0._ki,m1,0._ki,m1,0._ki,0._ki,0._ki,0._ki,m1,m1,0._ki,0._ki,0._ki,m1 /)
        set_m2 = (/ m2,0._ki,m2,0._ki,0._ki,m2,m2,0._ki,m2,0._ki,m2,m2,m2,0._ki,m2,m2,m2 /)
        set_m3 = (/ 0._ki,0._ki,0._ki,m3,m3,m3,m3,0._ki,0._ki,m3,m3,0._ki,m3,0._ki,0._ki,m3,m3 /)
        set_arg1 = (/ 0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,a1,a1,a1,a1 /)
        set_arg2 = (/ 0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,0._ki,a2,a2,a2,a2,a2,a2,a2,a2,a2,a2 /)
        set_arg3 = (/ 0._ki,a3,a3,a3,a3,a3,a3,a3,a3,a3,a3,a3,a3,a3,a3,a3,a3 /)
        !
        do i=1,17
           !
           s_mat(1,1) = -2._ki*set_m1(i)
           s_mat(2,2) = -2._ki*set_m2(i)
           s_mat(3,3) = -2._ki*set_m3(i)
           !              
           s_mat(1,3) = set_arg3(i)
           s_mat(2,3) = set_arg2(i)
           s_mat(1,2) = set_arg1(i)
           s_mat(3,1) = s_mat(1,3)
           s_mat(2,1) = s_mat(1,2)
           s_mat(3,2) = s_mat(2,3)
           !
           call preparesmatrix()
           !
           res = a30(b_pin)
           write (17,formatff) 'a30', res
           res = a31(1,b_pin)
           write (17,formatff) 'a31-1', res
           res = a31(2,b_pin)
           write (17,formatff) 'a31-2', res
           res = a31(3,b_pin)
           write (17,formatff) 'a31-3', res
           res = a32(1,1,b_pin)
           write (17,formatff) 'a32-11', res
           res = a32(1,2,b_pin)
           write (17,formatff) 'a32-12', res
           res = a32(1,3,b_pin)
           write (17,formatff) 'a32-13', res
           res = a32(2,2,b_pin)
           write (17,formatff) 'a32-22', res
           res = a32(2,3,b_pin)
           write (17,formatff) 'a32-23', res
           res = a32(3,3,b_pin)
           write (17,formatff) 'a32-33', res
           res = a33(1,1,1,b_pin)
           write (17,formatff) 'a33-111', res
           res = a33(1,1,2,b_pin)
           write (17,formatff) 'a33-112', res
           res = a33(1,1,3,b_pin)
           write (17,formatff) 'a33-113', res
           res = a33(1,2,2,b_pin)
           write (17,formatff) 'a33-122', res
           res = a33(1,2,3,b_pin)
           write (17,formatff) 'a33-123', res
           res = a33(1,3,3,b_pin)
           write (17,formatff) 'a33-133', res
           res = a33(2,2,2,b_pin)
           write (17,formatff) 'a33-222', res
           res = a33(2,2,3,b_pin)
           write (17,formatff) 'a33-223', res
           res = a33(2,3,3,b_pin)
           write (17,formatff) 'a33-233', res
           res = a33(3,3,3,b_pin)
           write (17,formatff) 'a33-333', res
           res = b32(b_pin)
           write (17,formatff) 'b32', res
           res = b33(1,b_pin)
           write (17,formatff) 'b33-1', res
           res = b33(2,b_pin)
           write (17,formatff) 'b33-2', res
           res = b33(3,b_pin)
           write (17,formatff) 'b33-3', res
           !
        end do
        !
     else
        !
        cz = cmplx(0._ki,0._ki,ki)
        set_m1c = (/ cz,cz,m1c,cz,m1c /)
        set_m2c = (/ cz,m2c,m2c,m2c,m2c /)
        set_m3c = m3c
        set_arg1c = (/ cz,cz,cz,a1c,a1c /)
        set_arg2c = a2c
        set_arg3c = a3c

        do i=1,5
           !
           s_mat(1,1) = -2._ki*set_m1c(i)
           s_mat(2,2) = -2._ki*set_m2c(i)
           s_mat(3,3) = -2._ki*set_m3c(i)

           s_mat(1,3) = set_arg3c(i)-set_m1c(i)-set_m3c(i)
           s_mat(2,3) = set_arg2c(i)-set_m2c(i)-set_m3c(i)
           if (set_arg1c(i) == cz) then
              s_mat(1,2) = cz
           else
              s_mat(1,2) = set_arg1c(i)-set_m1c(i)-set_m2c(i)
           end if
           s_mat(3,1) = s_mat(1,3)
           s_mat(2,1) = s_mat(1,2)
           s_mat(3,2) = s_mat(2,3)
           !
           !
           call preparesmatrix()
           !
           res = a30(b_pin)
           write (17,formatff) 'a30', res
           res = a31(1,b_pin)
           write (17,formatff) 'a31-1', res
           res = a31(2,b_pin)
           write (17,formatff) 'a31-2', res
           res = a31(3,b_pin)
           write (17,formatff) 'a31-3', res
           res = a32(1,1,b_pin)
           write (17,formatff) 'a32-11', res
           res = a32(1,2,b_pin)
           write (17,formatff) 'a32-12', res
           res = a32(1,3,b_pin)
           write (17,formatff) 'a32-13', res
           res = a32(2,2,b_pin)
           write (17,formatff) 'a32-22', res
           res = a32(2,3,b_pin)
           write (17,formatff) 'a32-23', res
           res = a32(3,3,b_pin)
           write (17,formatff) 'a32-33', res
           res = a33(1,1,1,b_pin)
           write (17,formatff) 'a33-111', res
           res = a33(1,1,2,b_pin)
           write (17,formatff) 'a33-112', res
           res = a33(1,1,3,b_pin)
           write (17,formatff) 'a33-113', res
           res = a33(1,2,2,b_pin)
           write (17,formatff) 'a33-122', res
           res = a33(1,2,3,b_pin)
           write (17,formatff) 'a33-123', res
           res = a33(1,3,3,b_pin)
           write (17,formatff) 'a33-133', res
           res = a33(2,2,2,b_pin)
           write (17,formatff) 'a33-222', res
           res = a33(2,2,3,b_pin)
           write (17,formatff) 'a33-223', res
           res = a33(2,3,3,b_pin)
           write (17,formatff) 'a33-233', res
           res = a33(3,3,3,b_pin)
           write (17,formatff) 'a33-333', res
           res = b32(b_pin)
           write (17,formatff) 'b32', res
           res = b33(1,b_pin)
           write (17,formatff) 'b33-1', res
           res = b33(2,b_pin)
           write (17,formatff) 'b33-2', res
           res = b33(3,b_pin)
           write (17,formatff) 'b33-3', res

        end do
        
     end if
     
     deallocate(set_arg1)
     deallocate(set_arg2)
     deallocate(set_arg3)
     deallocate(set_m1)
     deallocate(set_m2)
     deallocate(set_m3)
     deallocate(set_arg1c)
     deallocate(set_arg2c)
     deallocate(set_arg3c)
     deallocate(set_m1c)
     deallocate(set_m2c)
     deallocate(set_m3c)
     !
     call exitgolem95()
     !
     !
     write (*,*) '*********************************'
     write (*,*) ' '
     write (*,*) '4pt - function'
     write (17,*) '*********************************'
     write (17,*) ' '
     write (17,*) '4pt - function'
     !
     b_pin = 0
     !
     ! we set up a new s matrix
     !
     call initgolem95(4)
     !
!!! We scan over QCDLoop Boxes 8,10,12,13 and one finite box
     !
     allocate(set_arg1(6))
     allocate(set_arg2(6))
     allocate(set_m1(6))
     !
     m1c = mass_int_sq(1)
     m3c = mass_int_sq(6)
     !
     m1 = real(m1c,ki)
     m3 = real(m3c,ki)
     !
     a1c = cmplx(.238751_ki,0._ki,ki)! - m1c - m2c
     a2c = cmplx(1.23062_ki,0._ki,ki)! - m2c - m3c
     !
     a1 = real(a1c,ki)
     a2 = real(a2c,ki)
     !
     s_mat(1,1) = cz
     s_mat(1,3) = cmplx(.1329563_ki,0._ki,ki)
     s_mat(1,4) = cmplx(.573522_ki,0._ki,ki) - m3c
     s_mat(2,2) = cz
     s_mat(2,4) = cmplx(.1493572_ki,0._ki,ki) - m3c
     s_mat(3,4) = cmplx(1.96387_ki,0._ki,ki) - m3c
     s_mat(4,4) = -2._ki*m3c
     s_mat(3,1) = s_mat(1,3)
     s_mat(4,1) = s_mat(1,4)
     s_mat(4,2) = s_mat(2,4)
     s_mat(4,3) = s_mat(3,4)
     !
     set_arg1 = (/ 0._ki,a1,0._ki,a1,a1,a1 /)
     set_arg2 = (/ 0._ki,0._ki,0._ki,0._ki,a2,a2 /)
     set_m1 = (/ 0._ki,0._ki,m1,m1,m1,0._ki /)
     !
     do i=1,6 !!!QCDLoop Boxes + finite 
        !
        s_mat(3,3) = cmplx(-2._ki*set_m1(i),0._ki,ki)
        s_mat(2,3) = cmplx(set_arg1(i),0._ki,ki)
        s_mat(1,2) = cmplx(set_arg2(i),0._ki,ki)
        s_mat(3,2) = s_mat(2,3)
        s_mat(2,1) = s_mat(1,2)
        !
        if (i == 6 ) then
           s_mat(4,4) = 0._ki
           s_mat = cmplx(real(s_mat,ki),0._ki,ki)
        end if
        !        
        call preparesmatrix()
        !
        write (17,*) 'i: ', i
        !
        res = a40(b_pin)
        write (17,formatff) 'a40', res
        res = a41(1,b_pin)
        write (17,formatff) 'a41-1', res
        res = a41(2,b_pin)
        write (17,formatff) 'a41-2', res
        res = a41(3,b_pin)
        write (17,formatff) 'a41-3', res
        res = a41(4,b_pin)
        write (17,formatff) 'a41-4', res
        res = a42(1,1,b_pin)
        write (17,formatff) 'a42-11', res
        res = a42(1,2,b_pin)
        write (17,formatff) 'a42-12', res
        res = a42(1,3,b_pin)
        write (17,formatff) 'a42-13', res
        res = a42(1,4,b_pin)
        write (17,formatff) 'a42-14', res
        res = a42(2,2,b_pin)
        write (17,formatff) 'a42-22', res
        res = a42(2,3,b_pin)
        write (17,formatff) 'a42-23', res
        res = a42(2,4,b_pin)
        write (17,formatff) 'a42-24', res
        res = a42(3,3,b_pin)
        write (17,formatff) 'a42-33', res
        res = a42(3,4,b_pin)
        write (17,formatff) 'a42-34', res
        res = a42(4,4,b_pin)
        write (17,formatff) 'a42-44', res
        res = a43(1,1,1,b_pin)
        write (17,formatff) 'a43-111', res
        res = a43(1,1,2,b_pin)
        write (17,formatff) 'a43-112', res
        res = a43(1,1,3,b_pin)
        write (17,formatff) 'a43-113', res
        res = a43(1,1,4,b_pin)
        write (17,formatff) 'a43-114', res
        res = a43(1,2,2,b_pin)
        write (17,formatff) 'a43-122', res
        res = a43(1,2,3,b_pin)
        write (17,formatff) 'a43-123', res
        res = a43(1,2,4,b_pin)
        write (17,formatff) 'a43-124', res
        res = a43(1,3,3,b_pin)
        write (17,formatff) 'a43-133', res
        res = a43(1,3,4,b_pin)
        write (17,formatff) 'a43-134', res
        res = a43(1,4,4,b_pin)
        write (17,formatff) 'a43-144', res
        res = a43(2,2,2,b_pin)
        write (17,formatff) 'a43-222', res
        res = a43(2,2,3,b_pin)
        write (17,formatff) 'a43-223', res
        res = a43(2,2,4,b_pin)
        write (17,formatff) 'a43-224', res
        res = a43(2,3,3,b_pin)
        write (17,formatff) 'a43-233', res
        res = a43(2,3,4,b_pin)
        write (17,formatff) 'a43-234', res
        res = a43(2,4,4,b_pin)
        write (17,formatff) 'a43-244', res
        res = a43(3,3,3,b_pin)
        write (17,formatff) 'a43-333', res
        res = a43(3,3,4,b_pin)
        write (17,formatff) 'a43-334', res
        res = a43(3,4,4,b_pin)
        write (17,formatff) 'a43-344', res
        res = a43(4,4,4,b_pin)
        write (17,formatff) 'a43-444', res

        res = a44(1,1,1,1,b_pin)
        write (17,formatff) 'a44-1111', res
        res = a44(1,1,1,2,b_pin)
        write (17,formatff) 'a44-1112', res
        res = a44(1,1,1,3,b_pin)
        write (17,formatff) 'a44-1113', res
        res = a44(1,1,1,4,b_pin)
        write (17,formatff) 'a44-1114', res
        res = a44(1,1,2,2,b_pin)
        write (17,formatff) 'a44-1122', res
        res = a44(1,1,2,3,b_pin)
        write (17,formatff) 'a44-1123', res
        res = a44(1,1,2,4,b_pin)
        write (17,formatff) 'a44-1124', res
        res = a44(1,1,3,3,b_pin)
        write (17,formatff) 'a44-1133', res
        res = a44(1,1,3,4,b_pin)
        write (17,formatff) 'a44-1134', res
        res = a44(1,1,4,4,b_pin)
        write (17,formatff) 'a44-1144', res
        res = a44(1,2,2,2,b_pin)
        write (17,formatff) 'a44-1222', res
        res = a44(1,2,2,3,b_pin)
        write (17,formatff) 'a44-1223', res
        res = a44(1,2,2,4,b_pin)
        write (17,formatff) 'a44-1224', res
        res = a44(1,2,3,3,b_pin)
        write (17,formatff) 'a44-1233', res
        res = a44(1,2,3,4,b_pin)
        write (17,formatff) 'a44-1234', res
        res = a44(1,2,4,4,b_pin)
        write (17,formatff) 'a44-1244', res
        res = a44(1,3,3,3,b_pin)
        write (17,formatff) 'a44-1333', res
        res = a44(1,3,3,4,b_pin)
        write (17,formatff) 'a44-1334', res
        res = a44(1,3,4,4,b_pin)
        write (17,formatff) 'a44-1344', res
        res = a44(1,4,4,4,b_pin)
        write (17,formatff) 'a44-1444', res

        res = a44(2,2,2,2,b_pin)
        write (17,formatff) 'a44-2222', res
        res = a44(2,2,2,3,b_pin)
        write (17,formatff) 'a44-2223', res
        res = a44(2,2,2,4,b_pin)
        write (17,formatff) 'a44-2224', res
        res = a44(2,2,3,3,b_pin)
        write (17,formatff) 'a44-2233', res
        res = a44(2,2,3,4,b_pin)
        write (17,formatff) 'a44-2234', res
        res = a44(2,2,4,4,b_pin)
        write (17,formatff) 'a44-2244', res
        res = a44(2,3,3,3,b_pin)
        write (17,formatff) 'a44-2333', res
        res = a44(2,3,3,4,b_pin)
        write (17,formatff) 'a44-2334', res
        res = a44(2,3,4,4,b_pin)
        write (17,formatff) 'a44-2344', res
        res = a44(2,4,4,4,b_pin)
        write (17,formatff) 'a44-2444', res

        res = a44(3,3,3,3,b_pin)
        write (17,formatff) 'a44-3333', res
        res = a44(3,3,3,4,b_pin)
        write (17,formatff) 'a44-3334', res
        res = a44(3,3,4,4,b_pin)
        write (17,formatff) 'a44-3344', res
        res = a44(3,4,4,4,b_pin)
        write (17,formatff) 'a44-3444', res

        res = a44(4,4,4,4,b_pin)
        write (17,formatff) 'a44-4444', res

        res = b42(b_pin)
        write (17,formatff) 'b42', res
        res = b43(1,b_pin)
        write (17,formatff) 'b43-1', res
        res = b43(2,b_pin)
        write (17,formatff) 'b43-2', res
        res = b43(3,b_pin)
        write (17,formatff) 'b43-3', res
        res = b43(4,b_pin)
        write (17,formatff) 'b43-4', res
        res = b44(1,1,b_pin)
        write (17,formatff) 'b44-11', res
        res = b44(1,2,b_pin)
        write (17,formatff) 'b44-12', res
        res = b44(1,3,b_pin)
        write (17,formatff) 'b44-13', res
        res = b44(1,4,b_pin)
        write (17,formatff) 'b44-14', res
        res = b44(2,2,b_pin)
        write (17,formatff) 'b44-22', res
        res = b44(2,3,b_pin)
        write (17,formatff) 'b44-23', res
        res = b44(2,4,b_pin)
        write (17,formatff) 'b44-24', res
        res = b44(3,3,b_pin)
        write (17,formatff) 'b44-33', res
        res = b44(3,4,b_pin)
        write (17,formatff) 'b44-34', res

        res = c44(b_pin)
        write (17,formatff) 'c44', res

     end do
     !     
     close(17)
     !
  end if
  !
  call exitgolem95()
  !
end program main
