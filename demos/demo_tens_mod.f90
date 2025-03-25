module     reduction_module
   use precision_golem, only: ki
   use matrice_s

   implicit none

   private

   real(ki), dimension(6), parameter :: s6   = &
            & (/0.0_ki, .0_ki, 6._ki, 0._ki, 0.0_ki, 0.0_ki/)
   ! loop propagator masses
   complex(ki), dimension(6), parameter :: msq6 = (/ &
            & (0.0_ki, 0.0_ki), &
            & (.96_ki, -0.0234678_ki), &
            & (100.230_ki, 0.0_ki), &
!            & (8.0230698238746_ki, 0.0_ki), &
            & (8.0234343_ki, -20.03452_ki), &
            & (1.54987724_ki, -0.030000001_ki), &
            & (0.0_ki, 0.0_ki) &
   & /)
!   complex(ki), dimension(6), parameter :: msq6 = (/ &
!            & (0.0_ki, 0.0_ki), &
!            & (.0_ki, 0.0_ki), &
!            & (.0_ki, 0.0_ki), &
!            & (.0_ki, 0.0_ki), &
!            & (0.0_ki, 0.0_ki), &
!            & (0.0_ki, 0.0_ki) &
!   & /)

   real(ki), dimension(6, 0:3), public :: rvecs
   real(ki), dimension(6, 0:3) :: vecs
   complex(ki), dimension(6,6) :: inv_c, aai
   real(ki) :: error6
   logical, dimension(1:6), public :: pinched

   public :: init_kinematics, done_kinematics, numerator,sdot

   interface sdot
      module procedure sdot_rr
      module procedure sdot_rc
      module procedure sdot_cr
      module procedure sdot_cc
   end interface

contains

  subroutine     init_kinematics()
    implicit none
    
    integer :: i,j,k
    
    real(ki), parameter :: phi = 0.75_ki * 3.141596_ki
    real(ki), parameter :: chi = 0.25_ki * 3.141596_ki
    real(ki) :: E12, E1, E2, E3, E4, E5, E6, x
    
    ! E12 + E1 + E2 = 0
    ! E1^2 - x^2 = s1
    ! E2^2 - x^2 = s2
    ! E1^2 + E2^2 + 2*E1*E2 = E12^2
    ! s1 + x^2 + s2 + x^2 + 2*E1*E2 = E12^2
    ! 2*x^2 = E12^2 - s1 - s2 - 2*E1*(-E1-E12)
    ! 2*x^2 = E12^2 - s1 - s2 + 2*E1^2 + 2*E1*E12
    ! 0 = s1 + 0.5 * (E12^2 - s1 - s2) + E1*E12
    ! E1 = - 1/(2*E12) * (E12^2 + s1 - s2)
    
    
    E3 = -sqrt(1.0_ki+s6(3))
    E4 = -sqrt(1.0_ki+s6(4))
    E5 = -sqrt(1.0_ki+s6(5))
    
    vecs(3,:)  = (/E3,  0.0_ki,  0.8_ki,  0.6_ki/)
    vecs(4,:)  = (/E4,  0.6_ki,  0.0_ki, -0.8_ki/)
    vecs(5,:)  = (/E5,  sin(chi)*sin(phi), cos(chi)*sin(phi), -cos(phi)/)
    vecs(6,1:3)  = -vecs(5,1:3)-vecs(4,1:3)-vecs(3,1:3)
    
    
    x  = vecs(6,1)**2 + vecs(6,2)**2 + vecs(6,3)**2
    E6 = -sqrt(x+s6(6))
    vecs(6,0) = E6
    E12 = +E3+E4+E5+E6
    
    E1 = 0.5_ki * ((s6(2) - s6(1))/E12 - E12)
    E2 = - E12 - E1
    x = sqrt(0.5_ki) * sqrt(E12*E12 - s6(1) - s6(2) - 2.0_ki * E1 * E2)
    
    vecs(1,:)  = (/E1, 0.0_ki,    0.0_ki,         x/)
    vecs(2,:)  = (/E2, 0.0_ki,    0.0_ki,        -x/)
    
    !print*, "k1.k1", sdot(vecs(1,:), vecs(1,:))
    !print*, "k2.k2", sdot(vecs(2,:), vecs(2,:))
    !print*, "k3.k3", sdot(vecs(3,:), vecs(3,:))
    !print*, "k4.k4", sdot(vecs(4,:), vecs(4,:))
    !print*, "k5.k5", sdot(vecs(5,:), vecs(5,:))
    !print*, "k6.k6", sdot(vecs(6,:), vecs(6,:))
    
    rvecs(1,:) = vecs(1,:)
    do i = 2, 6
       rvecs(i,:) = rvecs(i-1,:) + vecs(i,:)
    end do
    
    call initgolem95(6)
    
    do i = 1, 6
       
       do j = 1, 6
          s_mat(i,j) = &!
               & sdot(rvecs(j,:)-rvecs(i,:),rvecs(j,:)-rvecs(i,:)) &
               & - msq6(i) - msq6(j)
          
          if (abs(s_mat(i,j)) .lt. 1.0E-14) then
             s_mat(i,j) = 0.0_ki
          end if
       end do
    end do
    !
    call preparesmatrix
    !
  end subroutine init_kinematics
  
  subroutine     done_kinematics()
    implicit none
    call exitgolem95()
  end subroutine done_kinematics
  
   function     numerator(Q, mu2) result(ans)
      implicit none
      real(ki), dimension(0:3), intent(in) :: Q
      real(ki), intent(in) :: mu2
      complex(ki) :: ans
      integer :: i

      ans = 1.0_ki

      do i = 1, 6
         if (pinched(i)) ans = ans * &
            & (sdot(Q(:)+rvecs(i,:),Q(:)+rvecs(i,:)) - msq6(i) - mu2)
      end do
   end function numerator
    
   pure function sdot_rr(v,w) result(r)
      implicit none
      real(ki), dimension(0:3), intent(in) :: v
      real(ki), dimension(0:3), intent(in) :: w
      real(ki) :: r

      r = v(0)*w(0) - v(1)*w(1) - v(2)*w(2) - v(3)*w(3)
   end  function sdot_rr

   pure function sdot_cc(v,w) result(r)
      implicit none
      complex(ki), dimension(0:3), intent(in) :: v
      complex(ki), dimension(0:3), intent(in) :: w
      complex(ki) :: r

      r = v(0)*w(0) - v(1)*w(1) - v(2)*w(2) - v(3)*w(3)
   end  function sdot_cc

   pure function sdot_rc(v,w) result(r)
      implicit none
      real(ki), dimension(0:3), intent(in) :: v
      complex(ki), dimension(0:3), intent(in) :: w
      complex(ki) :: r

      r = v(0)*w(0) - v(1)*w(1) - v(2)*w(2) - v(3)*w(3)
   end  function sdot_rc

   pure function sdot_cr(v,w) result(r)
      implicit none
      complex(ki), dimension(0:3), intent(in) :: v
      real(ki), dimension(0:3), intent(in) :: w
      complex(ki) :: r

      r = v(0)*w(0) - v(1)*w(1) - v(2)*w(2) - v(3)*w(3)
   end  function sdot_cr


   subroutine rnd_c0(arr)
      implicit none
      complex(ki), intent(out) :: arr
      real(ki) :: re, im

      call random_number(re)
      call random_number(im)
      arr = re + (0,1) * im
   end subroutine rnd_c0

   subroutine rnd_c_arr(arr, sz)
      implicit none
      complex(ki), dimension(:,:), intent(out) :: arr
      integer, dimension(2), intent(in) :: sz
      real(ki),dimension(sz(1),sz(2)) :: re, im

      call random_number(re)
      call random_number(im)
      arr = re + (0,1) * im
   end subroutine rnd_c_arr

end module reduction_module
