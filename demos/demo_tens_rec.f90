! this program illustrates the use of the function "evaluate"
! for a given numerator, "evaluate" combines the coefficients 
! C_{i,alpha} in eq.(15) with the corresponding tensor integrals
! 
! In the example below, the "amplitude" is given by a six-point integral
! where the numerator consists of three propagators which are also present 
! in the denominator.
! The scalar three-point function resulting from the direct cancellation of 
! these propagators is compared to the result obtained by expanding
! the numerator into contracted loop and external momenta, leading to 
! rank 6 hexagons if no propagators are cancelled, 
! rank 4 pentagons if one propagator is cancelled (3 possibilities to cancel 1 prop.),
! rank 2 boxes if two propagators are cancelled (3 possibilities to cancel 2 props.).
!
! the numerator is defined in the module reduction_module in demo_tens_mod.f90
!
program reduction_test
   use precision_golem
   use tens_rec
   use tens_comb
   use form_factor_type, only: form_factor, assignment(=)
   use reduction_module
   use array, only: packb, punion, pminus
   use matrice_s, only: b_ref
   use parametre
   implicit none

   integer :: prop1, prop2, prop3, prop4, prop5, prop6
   integer, dimension(6) :: prop_set
   integer :: b_prop_set_3, b_prop_set_4, b_prop_set_5, b_prop_set_6
   integer :: i, case4, case5

   type(form_factor) :: res3
   type(form_factor), dimension(3) :: res4
   type(form_factor), dimension(3) :: res5
   type(form_factor) :: res6

   if_print_warn_par = .true.

   call init_kinematics()

   do prop1 = 1, 4
      prop_set(1) = prop1
      do prop2 = prop1 + 1, 5
         prop_set(2) = prop2
         do prop3 = prop2 + 1, 6
            prop_set(3) = prop3
            pinched(:) = .false.
            b_prop_set_3 = packb(prop_set(1:3))
!
            res3 = evaluate(numerator, rvecs, pminus(b_ref, b_prop_set_3), 0)

            case4 = 1
            case5 = 1
            do prop4 = 1, 6
               if((prop1.eq.prop4).or.(prop2.eq.prop4).or.(prop3.eq.prop4)) &
               & cycle

               prop_set(4) = prop4
               b_prop_set_4 = packb(prop_set(1:4))
               pinched(prop4) = .true.
!	       
               res4(case4) = evaluate(numerator, rvecs, &
               & pminus(b_ref, b_prop_set_4), 2)
               case4 = case4 + 1

               do prop5 = 1, 6
                  if((prop1.eq.prop5).or.(prop2.eq.prop5).or.&
                  &  (prop3.eq.prop5).or.(prop4.ge.prop5)) cycle

                  prop_set(5) = prop5
                  b_prop_set_5 = packb(prop_set(1:5))
                  pinched(prop5) = .true.
!		  
                  res5(case5) = evaluate(numerator, rvecs, &
                  & pminus(b_ref, b_prop_set_5), 4)
                  case5 = case5 + 1
                  do prop6 = 1, 6
                     if((prop1.eq.prop6).or.(prop2.eq.prop6).or.&
                     &  (prop3.eq.prop6).or.(prop4.ge.prop6).or.&
                     &  (prop5.ge.prop6)) cycle!

                     prop_set(6) = prop6
                     b_prop_set_6 = packb(prop_set(1:6))
                     pinched(prop6) = .true.
!		     
                     res6 = evaluate(numerator, rvecs, &
                     & pminus(b_ref, b_prop_set_6), 6)
                     pinched(prop6) = .false.
                  end do
                  pinched(prop5) = .false.
               end do
               pinched(prop4) = .false.
            end do

            write(*,'(A11,1x,I1,I1,I1)'), "Propagators", prop1, prop2, prop3 
            write(*,'(A11)'), "Double Pole"
            write(*,'(A10,F24.16,1x,F24.16)'), "tri", res3%a
            do i = 1,3
               write(*,'(A9,I1,F24.16,1x,F24.16)'), "box", i, res4(i)%a
            end do
            do i = 1,3
               write(*,'(A9,I1,F24.16,1x,F24.16)'), "pentagon", i, res5(i)%a
            end do
            write(*,'(A10,F24.16,1x,F24.16)'), "hexagon", res6%a

            write(*,'(A11)'), "Single Pole"
            write(*,'(A10,F24.16,1x,F24.16)'), "tri", res3%b
            do i = 1,3
               write(*,'(A9,I1,F24.16,1x,F24.16)'), "box", i, res4(i)%b
            end do
            do i = 1,3
               write(*,'(A9,I1,F24.16,1x,F24.16)'), "pentagon", i, res5(i)%b
            end do
            write(*,'(A10,F24.16,1x,F24.16)'), "hexagon", res6%b

            write(*,'(A11)'), "Finite Part"
            write(*,'(A10,F24.16,1x,F24.16)'), "tri", res3%c
            do i = 1,3
               write(*,'(A9,I1,F24.16,1x,F24.16)'), "box", i, res4(i)%c
            end do
            do i = 1,3
               write(*,'(A9,I1,F24.16,1x,F24.16)'), "pentagon", i, res5(i)%c
            end do
            write(*,'(A10,F24.16,1x,F24.16)'), "hexagon", res6%c
         end do
      end do
   end do
   call done_kinematics()
end program reduction_test
