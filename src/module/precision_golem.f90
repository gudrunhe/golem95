! 
!****h* src/module/precision_golem
! NAME
!
!  Module precision_golem
!
! USAGE
!
!  use precision_golem
!
! DESCRIPTION
!
!  This module defines the parameter ki which gives the representation 
!  of the real and complex numbers in golem
!
! OUTPUT
!
!  The integer parameter ki
!  The integer parameter ki_avh, which is the real kind used in avh_olo
!  The integer parameter ki_lt, which is the real kind used in LoopTools
!
! USES
!
!  No uses
!
!*****
module precision_golem
    use avh_olo_dp_kinds, only: ki_avh => kindr2
    implicit none
    private
    !
    integer, parameter, public :: ki=kind(1.0d0)
!AC!integer, parameter, public :: ki_lt=kind(1.0d0)

    public :: ki_avh
end module precision_golem

