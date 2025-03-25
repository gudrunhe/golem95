! this program demonstrates how to call Golem form factors 
! using the same syntax as for the call of 
! LoopTools integrals, to facilitate comparisons 
! In addition, the last argument iep has been added 
! to indicate the order of the Laurent expansion in epsilon,
! in the spirit of QCDLoop (iep = -2,-1,0) 
!
! For the divergent integrals, a global prefactor 
! Gamma(1-eps)^2/Gamma(1-2*eps)*Gamma(1+eps)*(4*pi)^eps 
! has been extracted
!
!
      program main
      implicit none
      real*8 s12,s23,s1,s2,s3,s4,m1,m2,m3,m4,mu2
      real*8 s34,s45,s51,s5,m5
      real*8 s56,s61,s234,s345,s123,s6,m6
      integer iep
      complex*16 result
      complex*16 gb0i,gc0i,gd0i,ge0i,gf0i
      integer nlegs
      character*6 labels,fname
      parameter (labels = '123456')
c
      nlegs = 4
c      
      iep = 0
c
      mu2 = 1.d0
c
      open(unit=17,file='result_demoLT.txt',status='unknown')
c
      if (nlegs .eq. 6) then
            s12  = 100.d0
            s23  = -2.d0
            s34 = 6.d0
            s45 = 7.d0
            s56 = 4.d0
            s61 = -9.d0
            s123 = -8.d0
            s234 = 11.d0
            s345 = 13.d0
            s1 =  0.d0
            s2 =  0.d0
            s3 =  0.d0
            s4 =  0.d0
            s5 = 0.d0
            s6 = 0.d0
            m1 =  0.d0
            m2 =  0.d0
            m3 =  0.d0
            m4 =  0.d0
            m5 =  0.d0
            m6 =  0.d0
            result = gf0i('ff0',s1,s2,s3,s4,s5,s6,s12,s23,s34,s45,s56,
     #                  s61,s123,s234,s345,m1,m2,m3,m4,m5,m6,mu2,iep)
            write(*,*) 'F0i :',result
            write(17,*) 'F0i =',result
      elseif (nlegs .eq. 5) then
            s12  = 100.d0
            s23  = -2.d0
            s34 = 6.d0
            s45 = 7.d0
            s51 = -4.d0
            s1 =  0.d0
            s2 =  0.d0
            s3 =  0.d0
            s4 =  0.d0
            s5 = 0.d0
            m1 =  0.d0
            m2 =  0.d0
            m3 =  0.d0
            m4 =  0.d0
            m5 =  0.d0
            result = ge0i('ee0',s1,s2,s3,s4,s5,s12,s23,s34,s45,s51,
     #                  m1,m2,m3,m4,m5,mu2,iep)
            write(*,*) 'E0i :',result
            write(17,*) 'E0i =',result
      elseif (nlegs .eq. 4) then
            s12  = 100.d0
            s23  = -2.d0
            s1 =  3.d0
            s2 =  4.d0
            s3 =  5.d0
            s4 =  6.d0
            m1 =  0.d0
            m2 =  0.d0
            m3 =  0.d0
            m4 =  0.d0
            result = gd0i('dd0',s1,s2,s3,s4,s12,s23,m1,m2,m3,m4,mu2,iep)
            write(*,*) 'D0i :',result
            write(17,*) 'D0i =',result
      elseif (nlegs .eq. 3) then
            s1 =  3.d0
            s2 =  4.d0
            s3 =  5.d0
            m1 =  0.d0
            m2 =  0.d0
            m3 =  0.d0
            result = gc0i('cc0',s1,s2,s3,m1,m2,m3,mu2,iep)
            write(*,*) 'C0i :',result
            write(17,*) 'C0i =',result
      elseif (nlegs .eq. 2) then
            s1 =  3.d0
            m1 =  0.d0
            m2 =  0.d0
            result = gb0i('bb0',s1,m1,m2,mu2,iep)
            write(*,*) 'B0i :',result
            write(17,*) 'B0i =',result
      endif
c
      fname=labels(nlegs:nlegs)//'point'
c      
      if (iep.eq.0) then
      write(17,*)'    = golem result for finite part of ', fname
      elseif (iep.eq.-1) then
      write(17,*)'    = golem result for 1/eps part of ' , fname
      elseif (iep.eq.-2) then
      write(17,*)'    = golem result for 1/eps^2 part of ' , fname
      endif
      write(6,*) 'The result has been written to result_demoLT.txt'
c
      end
