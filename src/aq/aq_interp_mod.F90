! (C) Copyright 2021-2022 CERFACS.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module floats

   use aq_constants_mod
   !
   ! Floating point precision parameters
   !

   integer, parameter :: low = aq_single
   integer, parameter :: dbl = aq_real
   integer, parameter :: ext = 2*aq_real

   integer, parameter :: std = dbl

   real(std), parameter :: Zero   = 0.0_std
   real(std), parameter :: Half   = 0.5_std
   real(std), parameter :: One    = 1.0_std
   real(std), parameter :: Two    = 2.0_std
   real(std), parameter :: Three  = 3.0_std

   real(std), parameter :: Small = tiny(Zero)

end module floats

module MathLib

   implicit none

   public :: LINTERP1, INTERPOL, Invert, MatSolve, linterp_h

contains
   !
   function linear_int1(y0,x0,x1,nodata)
      !
      use floats
      implicit none
      !
      ! x0 is an array of pressure levels and the linear interpolation
      ! is performed in pressure
      !
      real(std), intent(in), dimension(:) :: y0,x0,x1
      real(std), intent(in)               :: nodata
      real(std), dimension(size(x1))      :: linear_int1
      !
      ! local
      !
      integer                   :: i,j,npt0,npt1
      real(std)                      :: z1,z2
      real(std)                      :: z0_j,z0_jp1
      real(std)                      :: y0_j,y0_jp1
      real(std)                      :: dz1,dz2
      real(std), dimension(size(x1)) :: y1
      !
      if(size(x0).ne.size(y0)) then
         print *,'incompatibility of dimensions in linear_int1'
         print *,size(x0),size(y0)
         call abort()
      endif
      !
      npt0 = size(x0)
      npt1 = size(x1)
      !
      y1 = nodata
      do i=2,npt1
         do j=1,npt0-1
            z0_j   = x0(j)
            y0_j   = y0(j)
            z0_jp1 = x0(j+1)
            y0_jp1 = y0(j+1)
            dz1 = x1(i) - z0_j
            dz2 = x1(i) - z0_jp1
            if((dz1*dz2).le.0.) then
               dz1 = abs(dz1) / abs(z0_jp1-z0_j)
               dz2 = abs(dz2) / abs(z0_jp1-z0_j)
               y1(i) = dz2 * y0_j + dz1 * y0_jp1
               exit
            endif
         end do
      end do
      !
      y1(1) = y0(npt0)
      !
      linear_int1 = y1
      !
      return
   end function linear_int1

   !
   function linear_int(y0,x0,x1,nodata)
      !
      use floats
      implicit none
      !
      ! x0 is an array of pressure levels and the linear interpolation
      ! is performed in pressure
      !
      real(std), intent(in), dimension(:) :: y0,x0,x1
      real(std), intent(in)               :: nodata
      real(std), dimension(size(x1))      :: linear_int
      !
      ! local
      !
      integer                   :: i,j,npt0,npt1
      real(std)                      :: z1,z2
      real(std)                      :: z0_j,z0_jp1
      real(std)                      :: y0_j,y0_jp1
      real(std)                      :: dz1,dz2
      real(std), dimension(size(x1)) :: y1
      !
      integer :: wk_npt0
      real(std), dimension(:), pointer :: wk_x0,wk_y0
      !
      if(size(x0).ne.size(y0)) then
         print *,'incompatibility of dimensions in linear_int'
         print *,size(x0),size(y0)
         call abort()
      endif
      !
      npt0 = size(x0)
      npt1 = size(x1)
      !
      wk_npt0 = 0
      do i=1,npt0
         if(y0(i).ne.nodata) then
            wk_npt0 = wk_npt0 + 1
         endif
      end do
      !
      allocate(wk_x0(wk_npt0))
      allocate(wk_y0(wk_npt0))
      wk_npt0 = 0
      do i=1,npt0
         if(y0(i).ne.nodata) then
            wk_npt0 = wk_npt0 + 1
            wk_x0(wk_npt0) = x0(i)
            wk_y0(wk_npt0) = y0(i)
         endif
      end do
      !
      y1 = nodata
      do i=1,npt1
         do j=1,wk_npt0
            z0_j = wk_x0(j)
            y0_j = wk_y0(j)
            if(j.le.(wk_npt0-1)) then
               z0_jp1 = wk_x0(j+1)
               y0_jp1 = wk_y0(j+1)
            else
               z0_jp1 = 0.
               y0_jp1 = 0.
            endif
            dz1 = x1(i) - z0_j
            dz2 = x1(i) - z0_jp1
            if((dz1*dz2).le.0.) then
               dz1 = abs(dz1) / abs(z0_jp1-z0_j)
               dz2 = abs(dz2) / abs(z0_jp1-z0_j)
               y1(i) = dz2 * y0_j + dz1 * y0_jp1
               exit
            endif
         end do
      end do
      !
      linear_int = y1
      !
      return
   end function linear_int
   !
   function linterp_h(y0,x0,x1,nodata)
      !
      use floats
      implicit none
      !
      ! x0 is an array of pressure levels and the linear interpolation
      ! is performed in pressure
      !
      real(std), intent(in), dimension(:) :: y0,x0,x1
      real(std), intent(in)               :: nodata
      real(std), dimension(size(x1))      :: linterp_h
      !
      ! local
      !
      integer                        :: i,j,npt0,npt1
      real(std)                      :: z1,z2,z0,dz1,dz2
      real(std), dimension(size(x1)) :: y1
      !
      if(size(x0).ne.size(y0)) then
         print *,'incompatibility of dimensions in linterp_h'
         print *,size(x0),size(y0)
         call abort()
      endif
      !
      npt0 = size(x0)
      npt1 = size(x1)
      !
      ! print *,'npt0 ',npt0,x0
      ! print *,'npt1 ',npt1,x1
      !
      y1 = nodata
      do j=1,npt0
         z0 = x0(j)
         do i=1,npt1-1
            dz1 = x1(i  ) - z0
            dz2 = x1(i+1) - z0
            if((dz1*dz2).le.0.) then
               dz1 = abs(dz1) / abs(x1(i)-x1(i+1))
               dz2 = abs(dz2) / abs(x1(i)-x1(i+1))
               y1(i  ) = y1(i  ) + dz2 * y0(j)
               y1(i+1) = y1(i+1) + dz1 * y0(j)
               exit
            endif
         end do
      end do
      !
      linterp_h = y1
      !
      return
   end function linterp_h
   !
   function LINTERP1( y, x, x1, NoData )
      !+
      ! PURPOSE:
      !      1-D Linear interpolation and extrapolation.
      !
      ! Input:
      ! y:       The input vector
      ! x:       The absicissae values for y
      ! x1:      The absicissae values for the result (y1)
      ! NoData:  Missing data value for y
      !
      ! Output:
      ! Returns vector of interpolated values. No extrapolation
      ! is performed if NoData is set - NoData value is assigned to
      ! the interpolated data outside of the x range.
      !-
      use floats
      implicit none
      !
      real(std), intent(IN) :: y(:), x(:), x1(:)
      real(std), optional, intent(IN) :: NoData
      real(std), dimension(size(x1)) :: LINTERP1
      !
      integer :: i, k, N
      integer, dimension(size(x)) :: Counter
      real(std), dimension(size(x1)) :: y1
      real(std), pointer, dimension(:) :: xx, yy
      real(std) :: xMin, xMax
      !
      if (size(x) .ne. size(y)) then
         print *, 'LInterp1% Incompatible x & y'
         stop
      endif

      if (present(NoData)) then

         Counter = 0

         where (y .ne. NoData)
            Counter = 1
         ENDWHERE

         N = sum( Counter )
         allocate( xx(N) )
         allocate( yy(N) )

         if (N .ge. 2) then
            k = 1
            do  i=1, size(x)
               if ( y(i) .ne. NoData ) then
                  xx(k) = x(i)
                  yy(k) = y(i)
                  k = k+1
               endif
            enddo
         else
            print *, 'LInterp1% Only one valid point.'
            y1 = NoData
            LINTERP1 = y1
            deallocate( xx )
            deallocate( yy )
            return
         endif

         y1=INTERPOL(yy,xx,x1)

         xMax = maxval(xx)
         xMin = minval(xx)

         where (x1 < xMin)
            y1 = NoData
         ENDWHERE

         where (x1 > xMax)
            y1 = NoData
         ENDWHERE
         deallocate( xx )
         deallocate( yy )
      else
         y1=INTERPOL(y,x,x1)
      endif

      LINTERP1 = y1

      !
      return
   end function LINTERP1

   function INTERPOL( V, X, U )
      !+
      ! NAME:
      ! INTERPOL
      !
      ! PURPOSE:
      ! Linearly interpolate and extrapolate vectors with a regular or irregular grid.
      !
      ! CATEGORY:
      ! E1 - Interpolation
      !
      ! CALLING SEQUENCE:
      ! Result = INTERPOL(V, X, U)
      !
      ! INPUTS:
      ! V:      The input vector can be any type except string.
      !
      ! X:      The absicissae values for V.  This vecotr must have same # of
      !         elements as V.  The values MUST be monotonically ascending
      !         or descending.
      !
      ! U:      The absicissae values for the result.  The result will have
      !         the same number of elements as U.  U does not need to be
      !         monotonic.
      !
      !
      ! OUTPUTS:
      ! INTERPOL returns a floating-point vector of N points determined
      ! by linearly interpolating the input vector.
      !
      ! PROCEDURE:
      ! Result(i) = V(x) + (x - FIX(x)) * (V(x+1) - V(x))
      !
      ! where        x = U(i).
      !         m = number of points of input vector.
      !
      ! MODIFICATION HISTORY:
      !                 IDL routine translated to F90 by Boris Khattatov, 9/25/96
      !-
      !
      use floats
      implicit none
      !
      real(std), intent(IN) :: V(:), X(:), U(:)
      real(std), dimension( size(U) ) :: INTERPOL
      !
      real(std), dimension( size(U) ) :: r
      real(std) :: s1, d
      integer :: M, N, ix, i
      !
      M = size(V)       !# of input pnts

      if (size(X) .ne. M) then
         print *, 'INTERPOL %  V and X must have same # of elements'
         stop
      endif

      N= size(U)        !# of output points

      r = V(1)  !floating, dbl or cmplx result

      if (X(2) - X(1) .ge. 0) then
         s1 = 1.0
      else
         s1 = -1.0
      endif !Incr or Decr X

      ix = 1                    !current point
      do i=1,N                     !point loop
         d = s1 * (U(i)-X(ix))  !difference
         if (abs(d) .lt. tiny(d)) then
            r(i)=V(ix)
         else   !at point
            if (d > 0) then
               do while ( (s1*(U(i)-X(ix+1)) > 0) .and. (ix < M-1) )
                  ix=ix+1
               enddo
            else
               do while ( (s1*(U(i)-X(ix)) < 0) .and.  (ix > 1) )
                  ix=ix-1
               enddo
            endif
            r(i) = V(ix) + (U(i)-X(ix))*(V(ix+1)-V(ix))/(X(ix+1)-X(ix))
         endif
      enddo

      INTERPOL = r
      !
      return
   end function INTERPOL

   function INTERPOL_2( V, X, U, nodata )
      !+
      ! NAME:
      ! INTERPOL_2
      !
      ! PURPOSE:
      ! Linearly interpolate  vectors with a regular or irregular grid.
      !       NO EXTRAPOLATION (AK)!
      !
      ! CATEGORY:
      ! E1 - Interpolation
      !
      ! CALLING SEQUENCE:
      ! Result = INTERPOL_2(V, X, U)
      !
      ! INPUTS:
      ! V:      The input vector can be any type except string.
      !
      ! X:      The absicissae values for V.  This vecotr must have same # of
      !         elements as V.  The values MUST be monotonically ascending
      !         or descending.
      !
      ! U:      The absicissae values for the result.  The result will have
      !         the same number of elements as U.  U does not need to be
      !         monotonic.
      !
      !
      ! OUTPUTS:
      ! INTERPOL_2 returns a floating-point vector of N points determined
      ! by linearly interpolating the input vector.
      !
      ! PROCEDURE:
      ! Result(i) = V(x) + (x - FIX(x)) * (V(x+1) - V(x))
      !
      ! where        x = U(i).
      !         m = number of points of input vector.
      !
      ! MODIFICATION HISTORY:
      !                 IDL routine translated to F90 by Boris Khattatov, 9/25/96
      !-
      !
      use floats
      implicit none
      !
      real(std), intent(IN) :: V(:), X(:), U(:)
      real(std), dimension( size(U) ) :: INTERPOL_2
      real(std), intent(in)               :: nodata
      !
      real(std), dimension( size(U) ) :: r
      real(std) :: s1, d
      integer :: M, N, ix, i
      !
      M = size(V)       !# of input pnts

      if (size(X) .ne. M) then
         print *, 'INTERPOL_2 %  V and X must have same # of elements'
         stop
      endif

      N= size(U)        !# of output points

      r = V(1)  !floating, dbl or cmplx result

      if (X(2) - X(1) .ge. 0) then
         s1 = 1.0
      else
         s1 = -1.0
      endif !Incr or Decr X

      r=nodata

      ix = 1                    !current point
      do i=1,N                     !point loop
         d = s1 * (U(i)-X(ix))  !difference
         if (abs(d) .lt. tiny(d)) then
            r(i)=V(ix)
         else   !at point
            if (d > 0) then
               do while ( (s1*(U(i)-X(ix+1)) > 0) .and. (ix < M-1) )
                  ix=ix+1
               enddo
            else
               do while ( (s1*(U(i)-X(ix)) < 0) .and.  (ix > 1) )
                  ix=ix-1
               enddo
            endif
            if ( (U(i)-X(ix))*(U(i)-X(ix+1)) .le. 0. ) then
               r(i) = V(ix) + (U(i)-X(ix))*(V(ix+1)-V(ix))/(X(ix+1)-X(ix))
            endif

         endif
      enddo

      INTERPOL_2 = r
      !
      return
   end function INTERPOL_2

   function Invert( N, A )
      !+
      !
      !-
      use floats
      implicit none
      !
      integer, intent(IN) :: N
      real(std), intent(IN) :: A(N,N)
      real(std), dimension(N,N) :: Invert
      !
      real(std), dimension(N,N) :: AInv, B, BInv
      real(std) :: d
      integer :: i,j, Indx(N), Flag
      !

      !  print *, 'In Invert.'

      do i=1,N
         do j=1,N
            BInv(i,j) = 0.
            B(i,j) = 0.
         enddo
         BInv(i,i) = 1.
      enddo

      do i=1,N
         do j=1,N
            B(i,j)=A(i,j)
         enddo
      enddo

      !      print *, ' CALL LUDCMP( B, N, N, Indx, d, Flag )'
      call LUDCMP( B, N, N, Indx, d, Flag )
      !      print *, 'Done'
      do i=1,N
         call LUBKSB( B, N, N, Indx, BInv(1,i) )
      enddo

      do i=1,N
         do j=1,N
            AInv(i,j)=BInv(i,j)
         enddo
      enddo

      Invert = AInv

      return
   end function Invert

   subroutine MatSolve( A, b, x)
      !+
      ! Solves A*x=b
      !-
      use floats
      implicit none
      !
      ! INPUTS:
      real(std), intent(IN) :: A(:,:), b(:)
      !
      ! OUTPUTS:
      real(std), intent(OUT) :: x(size(b))
      !
      ! INTERNALS:
      integer:: N
      real(std) :: d
      integer :: Indx(size(b)), Flag
      real(std) :: AA(size(b),size(b)), bb(size(b))
      !
      real(std) :: p(size(b))
      !
      N = size(b)
      !
      AA = A
      bb = b
      !
      call LUDCMP( AA, N, N, Indx, d, Flag )
      if(Flag.ne.0) call abort()
      call LUBKSB( AA, N, N, Indx, bb )
      x = bb
      !
      ! use Cholesky decomposition instead of LU-decomposition
      !
      ! call choldc(aa,n,n,p)
      ! call cholsl(aa,n,n,p,bb,x)
      !
      ! use conjugate gradient method
      !
      ! call conjugate_gradient(n,aa,b,x)
      !
      return
   end subroutine MatSolve

   subroutine LUBKSB(a,n,np,indx,b)
      !+
      !
      !-
      use floats
      implicit none
      !
      integer n,np,indx(n)
      real(std) a(np,np),b(n)
      !
      integer i,ii,j,ll
      real(std) sum_
      !
      ii=0
      do i=1,n
         ll=indx(i)
         sum_=b(ll)
         b(ll)=b(i)
         if (ii  .ne.  0) then
            do  j=ii,i-1
               sum_=sum_-a(i,j)*b(j)
            enddo
         else if (sum_  .ne.  0.) then
            ii=i
         endif
         b(i)=sum_
      enddo
      do i=n,1,-1
         sum_=b(i)
         do  j=i+1,n
            sum_=sum_-a(i,j)*b(j)
         enddo
         b(i)=sum_/a(i,i)
      enddo
      !
      return
   end subroutine LUBKSB



   subroutine LUDCMP(a,n,np,indx,d, Flag)
      !+
      !
      !-
      use floats
      implicit none
      !
      integer :: n,np,indx(n)
      real(std) :: d,a(np,np)
      integer, intent(OUT) :: Flag
      !
      integer i,imax,j,k
      real(std) aamax,dum,sum_,vv(np)
      !
      Flag = 0

      d=1.
      do i=1,n
         aamax=0.
         do  j=1,n
            if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j))
         enddo
         if (aamax .eq. 0.) then
            print *,  'Singular matrix in ludcmp ',i
            print *, a(i,:)
            Flag = -1
            return
         endif
         vv(i)=1./aamax
      enddo
      do  j=1,n
         do  i=1,j-1
            sum_=a(i,j)
            do  k=1,i-1
               sum_=sum_-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum_
         enddo
         aamax=0.
         do i=j,n
            sum_=a(i,j)
            do k=1,j-1
               sum_=sum_-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum_
            dum=vv(i)*abs(sum_)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo
         if (j  .ne.  imax) then
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            enddo
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if (a(j,j) .eq. 0.) a(j,j)=1.e-30
         if (j  .ne.  n) then
            dum=1./a(j,j)
            do  i=j+1,n
               a(i,j)=a(i,j)*dum
            enddo
         endif
      enddo
      !
      return
   end subroutine LUDCMP
   !
   ! the following subroutines for the Cholesky decompostion
   ! come from the Numerical Recipes (chapter 2.9)
   !
   subroutine choldc(a,n,np,p)
      !
      use floats
      implicit none
      !
      integer n,np
      real(std) a(np,np),p(n)
      integer i,j,k
      real(std) wk_sum
      !
      do i=1,n
         do j=i,n
            wk_sum=a(i,j)
            do k=i-1,1,-1
               wk_sum=wk_sum-a(i,k)*a(j,k)
            enddo
            if(i.eq.j)then
               if(wk_sum.le.0.) then
                  print *,'choldc failed ',wk_sum
                  do k=i,1,-1
                     print *,i,k,a(i,k)
                  end do
                  write(110) n
                  write(110) a
                  call abort()
               endif
               p(i)=sqrt(wk_sum)
            else
               a(j,i)=wk_sum/p(i)
            endif
         enddo
      enddo
      return
   end subroutine choldc
   !
   subroutine cholsl(a,n,np,p,b,x)
      !
      use floats
      implicit none
      !
      integer n,np
      real(std) a(np,np),b(n),p(n),x(n)
      integer i,k
      real(std) wk_sum
      do i=1,n
         wk_sum=b(i)
         do k=i-1,1,-1
            wk_sum=wk_sum-a(i,k)*x(k)
         enddo
         x(i)=wk_sum/p(i)
      enddo
      do i=n,1,-1
         wk_sum=x(i)
         do k=i+1,n
            wk_sum=wk_sum-a(k,i)*x(k)
         enddo
         x(i)=wk_sum/p(i)
      enddo
      return
   end subroutine cholsl
   !
   subroutine conjugate_gradient(n,a,b,x)
      !
      ! based on algorithm 10.2.7 from Golub and van Loon
      !
      use floats
      implicit none
      !
      ! arguments
      !
      integer n
      real(std) a(n,n),b(n),x(n)
      !
      ! local
      !
      integer i,j,niter,niter_max
      real(std) r(n),p(n),w(n),sum
      real(std) alpha,beta,epsilon,norm_b
      real(std), dimension(:), pointer :: rho
      !
      ! parameters
      !
      niter_max = 100
      epsilon   = 1.e-8
      allocate(rho(0:niter_max))
      !
      ! set x to 0, r to b
      !
      do i=1,n
         x(i) = 0.
         r(i) = b(i)
      end do
      !
      ! set counter to 0
      !
      niter = 0
      !
      ! iterate
      !
100   continue
      !
      ! compute 2-norm of r
      !
      rho(niter) = 0.d0
      do i=1,n
         rho(niter) = rho(niter) + (r(i)*r(i))
      end do
      !
      ! check convergence
      !
      ! ! rho(0) = ||b||^2
      !
      if(sqrt(rho(niter)).lt.(epsilon*sqrt(rho(0)))) then
         print    *,'Conjugate-gradient has converged after ',niter &
            ,' iterations'
         return
      endif
      !
      niter = niter + 1
      !
      ! test number of iterations
      !
      if(niter.gt.niter_max) then
         print *,'reached maximum number of iterations in CG, stopping'
         call abort()
      endif
      !
      if(niter.eq.1) then
         do i=1,n
            p(i) = r(i)
         end do
      else
         beta = rho(niter-1)/rho(niter-2)
         do i=1,n
            p(i) = r(i) + beta * p(i)
         end do
      endif
      !
      do i=1,n
         w(i) = 0.
         do j=1,n
            w(i) = w(i) + a(i,j) * p(j)
         end do
      end do
      sum = 0.
      do i=1,n
         sum = sum + p(i) * w(i)
      end do
      alpha = rho(niter-1)/sum
      do i=1,n
         x(i) = x(i) + alpha * p(i)
         r(i) = r(i) - alpha * w(i)
      end do
      !
      goto 100
      !
   end subroutine conjugate_gradient

   subroutine get_coef_interpol( X, U , weights, positions )
      ! NAME:
      !       get_coef_interpol
      !
      ! PURPOSE:
      !       this is a version of the interpolation routine INTERPOL contained
      !       in the original assimilation package from NCAR. Instead of returning
      !       the interpolated array, the routine returns the arrays with the
      !       obtained weights and positions.
      !
      !- revision 1.1 Jan 11, 2006: the routine now handles cases with 0 in the U array.
      ! this corresponds to cases when only a part of the U vector is used (the number of
      ! observed pressure levels is smaller than the dimension of the array.
      !
      use floats
      implicit none
      !
      ! INPUT VARIABLES:
      real(std), intent(IN) ::  X(:), U(:)
      ! OUTPUT VARIABLES:
      real(std), intent(out) :: weights(2*size(U))
      integer,intent(out) :: positions(2*size(U))
      !
      real(std) :: s1, d, factor,si1,si2
      integer :: M, N, ix, i, i1, i2
      !
      M = size(X)    !# of input pnts

      N= size(U)     !# of output points

      if (abs(X(2)-X(1)) <= Small) then
         write(*,*) 'ERROR in get_coef_interpol (in mathlib.f90). '
         write(*,*) 'Two input levels have the same value, cant decide if the order'
         write(*,*) 'is increasing or decreasing monotonically. EXITING...'
      else
         s1 =sign(1._std,X(2)-X(1))
      endif !Incr or Decr X

      ix = 1                 !current point
      do i=1,N                          !point loop
         i1=(i-1)*2+1
         i2=(i-1)*2+2
         if (abs(U(i)) > Small) then  !if the pressure is at 0, no interpolation at this level
            d = s1 * (U(i)-X(ix))       !difference
            if (abs(d) < Small) then
               positions(i1)=ix
               positions(i2)=ix
               weights(i1)=1.
               weights(i2)=0.
            else
               if (d > 0) then
                  do while ( (s1*(U(i)-X(ix+1)) > 0) .and. (ix < M-1) )
                     ix=ix+1
                  enddo
               else
                  do while ( (s1*(U(i)-X(ix)) < 0) .and.  (ix > 1) )
                     ix=ix-1
                  enddo
               endif
               positions(i1)=ix
               positions(i2)=ix+1
               ! no extrapolation :
               if (abs(X(ix)-U(i)) <= Small ) then
                  weights(i1)=1.
                  weights(i2)=0.
               elseif (abs(X(ix+1)-U(i)) <= Small) then
                  weights(i1)=0.
                  weights(i2)=1.
               else
                  si1 = sign(1.0_std,(X(ix+1)-U(i)))
                  si2 = sign(1.0_std,(X(ix)-U(i)))
                  !               si1 = (X(ix+1)-U(i))/ABS((X(ix+1)-U(i)))
                  !               si2 = (X(ix)-U(i))/ABS(X(ix)-U(i))
                  if ( (si1*si2) <= 0. ) then
                     if (abs(X(ix+1)-X(ix)) < Small) then
                        weights(i1)=1.
                        weights(i2)=0.
                     else
                        factor=1./(X(ix+1)-X(ix))
                        weights(i1)=(X(ix+1)-U(i))*factor
                        weights(i2)=(U(i)-X(ix))*factor
                     endif
                  else   !si1 and si2 have the same sign
                     if (si1.gt.0) then  !obs above the top model level
                        if (X(ix)>X(ix+1)) then
                           weights(i1)=0.
                           weights(i2)=1.
                        else
                           weights(i1)=1.
                           weights(i2)=0.
                        end if
                     else                !obs below the bottom model level
                        if (X(ix)>X(ix+1)) then
                           weights(i1)=1.
                           weights(i2)=0.
                        else
                           weights(i1)=0.
                           weights(i2)=1.
                        end if
                     end if

                  end if
               end if
            endif
         else
            positions(i1)=1  !arbitrary positions
            positions(i2)=2
            weights(i1)=0.
            weights(i2)=0.
         end if
      enddo
      !
      return
   end subroutine get_coef_interpol


   !
end module MathLib
module INTERP_MATRIX_STRUCTURE_MOD

   use floats

   implicit none

   integer, parameter :: ip_2d = 0
   integer, parameter :: ip_interp = 1
   integer, parameter :: ip_ground = 2
   integer, parameter :: ip_identity = 3

   type CSR_FORMAT
      integer*4                       :: NPRF
      integer*4                       :: ncol
      integer*4                       :: nrow
      integer*4                       :: nmax
      real(std),dimension(:,:),allocatable  :: H
      integer*4,dimension(:,:),allocatable  :: i
      integer*4,dimension(:,:),allocatable  :: j
   end type CSR_FORMAT

   type MODEL_DATA
      integer*4  :: nlon       !numb. of model longitudes
      integer*4  :: nlat       !numb. of model latitudes
      integer*4  :: nlev     !numb. of model veritcal levels

      real(std),dimension(2)                     :: date
      real(std), pointer                         :: lat(:,:)
      real(std), pointer                         :: lon(:,:)
      real(std), pointer                         :: traceur(:,:,:,:)
      real(std), pointer                         :: pres(:,:,:,:)
   end type MODEL_DATA

end module INTERP_MATRIX_STRUCTURE_MOD
module matrix_manipulations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !this module contains routines that facilitate manipulating CSR matrices.
   ! calls the csr libraries.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! VERSION: 1.0
   ! WRITTEN BY: NOVELTIS (A.KLONECKI)
   ! CREATION DATE: 26/09/2005
   ! REVISION DATE:

   private atmux

contains

   subroutine multiply_matrices_csr( iw             & !work array
      ,Ha_csr                  & !matrix A
      ,Hb_csr                  & !matrix B
      ,Hc_csr                  & !matrix C=A*B
      ,nmaxc                   &
      )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! routine that multiply two matrices in CSR format, it uses the routine amub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use floats
      use interp_matrix_structure_mod


      implicit none
      !
      !INPUT:
      integer, dimension(:)                  :: iw
      type(CSR_FORMAT), intent(IN)           :: Ha_csr
      type(CSR_FORMAT), intent(IN)           :: Hb_csr
      integer, intent(IN)                    :: nmaxc
      !OUTPUT:
      type(CSR_FORMAT), intent(OUT)          :: Hc_csr
      !LOCAL:
      integer*4        :: job,ierr,np,NPRF,ALLOC_ERR


      job = 1
      ! the problem with the iw variable is that it is dimensioned by the number of
      ! columns of the full matrix which is large: 1522800 for 90*180*47*2
      ! To avoid intializing iw to zero in amub (as coded in SPARSKIT, and which takes time),
      ! the initialization is done only once in this routine. The variable iw is reinitialized to
      ! zero in amub2 uniquely for the indicies which were modified.

      if (allocated(Hc_csr%H)) then
         call deallocate_operator(Hc_csr)
      end if

      NPRF=Ha_csr%NPRF

      Hc_csr%NPRF  = NPRF
      Hc_csr%nrow  = Ha_csr%nrow
      Hc_csr%ncol  = Hb_csr%ncol
      Hc_csr%nmax  = nmaxc
      !
      allocate(Hc_csr%H(nmaxc,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"multiply_matrices_csr (")
      allocate(Hc_csr%j(nmaxc,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"multiply_matrices_csr (")
      allocate(Hc_csr%i(Hc_csr%nrow+1,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"multiply_matrices_csr (")
      !
      !      Hc_csr%H=0.
      !      Hc_csr%j=0
      !      Hc_csr%i=0
      do np=1,NPRF

         !amub2 is a SPARSKIT routine modified by NOVELTIS to run faster. For details
         !see the routine amub2 in this file.
         call amub2(Ha_csr%nrow,Hb_csr%ncol,job,              &
            Ha_csr%H(:,np),Ha_csr%j(:,np),Ha_csr%i(:,np), &    !matrix A
            Hb_csr%H(:,np),Hb_csr%j(:,np),Hb_csr%i(:,np), &    !matrix B
            Hc_csr%H(:,np),Hc_csr%j(:,np),Hc_csr%i(:,np), &    !out:A*B
            nmaxc,iw,ierr)
         call check_error_sparse(ierr,"amub")
      end do

      return
   end subroutine multiply_matrices_csr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine multiply_matrix_csr_vector(           &
                                !  input
      H_csr                   & !matrix in format csr
      ,vector                  & !vector
      ,nrow                    & !dim 1
      ,NPRF                    & !dim 2
                                !  output
      ,vector_out              & !the resulting vector
      )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! routine that multiplies a matrix in CSR format by a vector, uses the routine amux
      ! from SPARSKIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use floats
      use interp_matrix_structure_mod


      implicit none
      !
      !INPUT:
      type(CSR_FORMAT)                       :: H_csr
      real(std),dimension(:)                 :: vector
      integer                                :: nrow
      integer                                :: NPRF
      !OUTPUT:
      real(std),dimension(nrow,NPRF)         :: vector_out
      !LOCAL:
      integer*4        :: np,ALLOC_ERR,ib_i,ib_k

      vector_out=0.
      do np=1,NPRF
         do ib_i = 1,H_csr%nrow
            do ib_k = H_csr%i(ib_i,np), H_csr%i(ib_i+1,np)-1
               vector_out(ib_i,np) = vector_out(ib_i,np) &
                  + H_csr%H(ib_k,np)*vector(H_csr%j(ib_k,np))
            end do
         end do
      end do

      return
   end subroutine multiply_matrix_csr_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine multiply_matrixT_csr_vector(           &
                                !  input
      H_csr                   & !matrix in format csr
      ,vector                  & !vector
      ,nrow                    & !dim 1
      ,NPRF                    & !dim 2
                                !  output
      ,vector_out              & !the resulting vector
      )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! routine that multiplies the transpose of a matrix in CSR format by a vector,
      ! uses the routine atmux
      ! from SPARSKIT
      ! written by Andrea Piacentini
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use floats
      use interp_matrix_structure_mod


      implicit none
      !
      !INPUT:
      type(CSR_FORMAT)                       :: H_csr
      real(std),dimension(:)                 :: vector
      integer                                :: nrow
      integer                                :: NPRF
      !OUTPUT:
      real(std),dimension(nrow,NPRF)         :: vector_out
      !LOCAL:
      integer*4        :: np,ALLOC_ERR

      do np=1,NPRF
         call atmux(H_csr%nrow,vector,vector_out(:,np),           &
            H_csr%H(:,np),H_csr%j(:,np),H_csr%i(:,np))
      end do

      return
   end subroutine multiply_matrixT_csr_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine addmult_matrixT_csr_vector(           &
                                !  input
      H_csr                   & !matrix in format csr
      ,vector                  & !vector
      ,ncol                    & !dim 1
      ,NPRF                    & !dim 2
                                !  output
      ,vector_out              & !the resulting vector
      ,NoData                  & !the bad value used for masking
      )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! routine that multiplies the transpose of a matrix in CSR format by a vector,
      ! uses the routine atmux
      ! from SPARSKIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use floats
      use interp_matrix_structure_mod


      implicit none
      !
      !INPUT:
      type(CSR_FORMAT)                       :: H_csr
      real(std),dimension(ncol,NPRF)         :: vector
      integer                                :: ncol
      integer                                :: NPRF
      real(std), optional                    :: NoData
      !OUTPUT:
      real(std),dimension(:),intent(inout)   :: vector_out
      !LOCAL:
      integer*4        :: np,ALLOC_ERR,ib_i,ib_k

      if (present(NoData)) then
         do np=1,NPRF
            do ib_i = 1,H_csr%nrow
               do ib_k = H_csr%i(ib_i,np), H_csr%i(ib_i+1,np)-1
                  if (vector(ib_i,np) /= NoData) &
                     vector_out(H_csr%j(ib_k,np)) = vector_out(H_csr%j(ib_k,np)) &
                     + vector(ib_i,np)*H_csr%H(ib_k,np)
               end do
            end do
         end do
      else
         do np=1,NPRF
            do ib_i = 1,H_csr%nrow
               do ib_k = H_csr%i(ib_i,np), H_csr%i(ib_i+1,np)-1
                  vector_out(H_csr%j(ib_k,np)) = vector_out(H_csr%j(ib_k,np)) &
                     + vector(ib_i,np)*H_csr%H(ib_k,np)
               end do
            end do
         end do
      end if
      return
   end subroutine addmult_matrixT_csr_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine copy_matrix_csr(Ha_csr    &
      ,Hb_csr)

      use floats
      use interp_matrix_structure_mod

      implicit none
      !
      !INPUT:
      type(CSR_FORMAT)                       :: Ha_csr
      ! OUTPUT
      type(CSR_FORMAT)                       :: Hb_csr
      ! LOCAL
      integer*4  :: i,np

      if (allocated(Hb_csr%H)) then
         call deallocate_operator(Hb_csr)
      end if
      !
      Hb_csr%nrow=Ha_csr%nrow
      Hb_csr%ncol=Ha_csr%ncol
      Hb_csr%nmax=Ha_csr%nmax
      Hb_csr%NPRF=Ha_csr%NPRF
      allocate(Hb_csr%H(Hb_csr%nmax,Hb_csr%NPRF))
      allocate(Hb_csr%j(Hb_csr%nmax,Hb_csr%NPRF))
      allocate(Hb_csr%i(Hb_csr%nrow+1,Hb_csr%NPRF))
      !
      do np=1,Ha_csr%NPRF
         do i=1,Hb_csr%nmax
            Hb_csr%H(i,np)=Ha_csr%H(i,np)
            Hb_csr%j(i,np)=Ha_csr%j(i,np)
         end do
         do i=1,Hb_csr%nrow+1
            Hb_csr%i(i,np)=Ha_csr%i(i,np)
         end do
      end do
      return
   end subroutine copy_matrix_csr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SPARSKIT routine modified slightly to run faster. The intialization
   ! of the array iw is very time consuming for arrays with large number
   ! of columns. Instead of intializing the array iw each time the routine
   ! is called, the inialization is done once in the calling routine. The
   ! indicies used by amub2 are reset to zero at the end of the routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine amub2 (nrow,ncol,job,a,ja,ia,b,jb,ib,  &
      c,jc,ic,nzmax,iw,ierr)
      integer nrow, ncol, job, nzmax, ierr
      real*8 a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
      real*8 scal
      integer ii, jj, jcol, jpos, k, ka, kb, len
      logical values
      values = (job .ne. 0)
      len = 0
      ic(1) = 1
      ierr = 0

      !      do j=1, ncol
      !         iw(j) = 0
      !      enddo
      !c
      do ii=1, nrow
         !c     row i
         do ka=ia(ii), ia(ii+1)-1
            if (values) scal = a(ka)
            jj   = ja(ka)
            do kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
               if (jpos .eq. 0) then
                  len = len+1
                  if (len .gt. nzmax) then
                     ierr = ii
                     return
                  endif
                  jc(len) = jcol
                  iw(jcol)= len
                  if (values) c(len)  = scal*b(kb)
               else
                  if (values) c(jpos) = c(jpos) + scal*b(kb)
               endif
            end do
         end do
         do k=ic(ii), len
            iw(jc(k)) = 0
         end do
         ic(ii+1) = len+1
      end do
      return
   end subroutine amub2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine check_error_sparse(ierr,str_routine)
      implicit none

      integer,intent(in)        ::ierr
      character(*),intent(in)   ::str_routine

      if (ierr /= 0) then
         write(*,*) 'Pb in the sparse matrix routine ',str_routine(1:len_trim(str_routine))
         write(*,*) 'error number ierr=',ierr,' exiting...'
         write(*,*) 'check if the number of non-zero elements in the result matrix not too small'
         call exit(-1)
      end if
      return
   end subroutine check_error_sparse
!!!!!!!!!!!!!!!!

   subroutine check_allocate(ierr,str_routine)
      implicit none

      integer,intent(in)        ::ierr
      character(*),intent(in)   ::str_routine

      if (ierr /= 0) then
         write(*,*) 'Pb while allocating memory in routine ',str_routine(1:len_trim(str_routine))
         write(*,*) 'error ierr= ',ierr
         write(*,*) 'Check if not out of memory, exiting...'
         call exit(-1)
      end if
      return
   end subroutine check_allocate



   subroutine deallocate_operator(H_csr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THIS ROUTINE DEALLOCATES THE ARRAYS IN THE CSR_FORMAT structure
      use interp_matrix_structure_mod

      implicit none
      !
      !INPUT:
      type(CSR_FORMAT)                       :: H_csr

      if(allocated(H_csr%H)) deallocate(H_csr%H)
      if(allocated(H_csr%j)) deallocate(H_csr%j)
      if(allocated(H_csr%i)) deallocate(H_csr%i)

   end subroutine deallocate_operator


   subroutine atmux ( n, x, y, a, ja, ia )

      !*****************************************************************************80
      !
      !! ATMUX computes A' * x for a CSR matrix A.
      !
      !  Discussion:
      !
      !    This routine multiplies the transpose of a matrix by a vector when the
      !    original matrix is stored in compressed sparse row storage. Can also be
      !    viewed as the product of a matrix by a vector when the original
      !    matrix is stored in the compressed sparse column format.
      !
      !  Modified:
      !
      !    07 January 2004
      !
      !  Author:
      !
      !    Youcef Saad
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
      !
      !    Input, real X(*), an array whose length is equal to the
      !    column dimension of A.
      !
      !    Output, real Y(N), the product A' * X.
      !
      !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
      !    Compressed Sparse Row format.
      !
      implicit none

      integer ( kind = 4 ) n

      real ( kind = 8 ) a(*)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ia(*)
      integer ( kind = 4 ) ja(*)
      integer ( kind = 4 ) k
      real ( kind = 8 ) x(*)
      real ( kind = 8 ) y(n)

      y(1:n) = 0.0D+00

      do i = 1, n
         do k = ia(i), ia(i+1)-1
            y(ja(k)) = y(ja(k)) + x(i) * a(k)
         end do
      end do

      return
   end subroutine atmux


end module matrix_manipulations


module transformations_matrix
   !-------------------------------------------
   !this module contains the subroutines used in the transformation of the MOCAGE
   ! data from model space to observational space. The routines return the
   ! interpolation matrices corresponding to the interpolation matrices. The matrices
   ! are represented in the Compressed Sparse Row (CSR) format to avoid creating
   ! and manipulating the matrices with full model dimensions.
   !
   !-------------------------------------------
   ! CSR FORMAT:
   ! IN THE CSR FORMAT EACH MATRIX IS REPRESENTED BY THREE ARRAYS (see the document
   ! contained in the CSR documentation):
   ! -THE REAL ARRAY CONTAINING THE NON-ZERO ELEMENTS (SIZE NNZ)
   ! -INTEGER ARRAY (JZ) CONTAINING THE COLUMN INDEX FOR EACH NON-ZERO ELEMENT
   !  (SIZE NNZ)
   ! -INTEGER ARRAY (IZ) CONTAINING THE POINTERS (POSITIONS) TO THE ELEMENTS IN
   !  THE REAL ARRAY. THE DIMENSIONS ARE NLEV+1, WITH THE FIRST NLEV ELEMENTS POINTING
   !  TO THE LOCATION OF THE FIRST ELEMENT IN EACH OF THE ROWS OF THE FULL MATRIX. THE
   !  ELEMENT N+1 CONTAINS : IA(1)+NNZ
   !-------------------------------------------

   !
   ! VERSION: 1.0
   ! WRITTEN BY: NOVELTIS (A.KLONECKI)
   ! CREATION DATE: 26/09/2005
   ! REVISION DATE:

contains

   !
   !******************************************************************
   ! THE ROUTINE transf_horiz_interp_matrix IS RESPONSIBLE FOR PREPARING THE MATRIX
   ! FOR THE HORIZONTAL INTERPOLATION OF THE MODEL DATA TOWARDS THE HORIZONTAL GRID
   ! OF THE OBSERVATIONS. THE INPUT DATA IS GIVEN FOR ALL MODEL LEVELS AND FOR ONE
   ! OR TWO MOCAGE TIME PERIODS.THE OUTPUT MATRIX IS GIVEN IN THE CRS FORMAT
   !******************************************************************
   !
   subroutine transf_horiz_interp_matrix(NPRF,nlon,nlat,nlev,ntime,        &
      latmod,lonmod,lat_o,lon_o,        &
      periodic_lon,                     &
      discard_obs,                      &
      H_csr)
      use floats
      use interp_matrix_structure_mod
      use matrix_manipulations

      implicit none

      ! INPUT VARIABLES:

      integer  ,intent(in)  :: NPRF       !numb. of obser. profiles present in
      !the current MOCAGE interval
      integer*4,intent(in)  :: nlon,nlat  !numb. of model longitudes and latitudes
      integer*4,intent(in)  :: nlev       !numb. of model vertical levels
      integer*4,intent(in)  :: ntime      !numb. of model outputs (1 if no temporal
      !interpolation, 2 if temp. interpolation
      !between two model outputs.
      ! Model lats and longs
      real(std), dimension(nlat),intent(in)         :: latmod
      real(std), dimension(nlon),intent(in)         :: lonmod
      ! IASI lats and longs
      real(std),dimension(NPRF),intent(in)                 :: lat_o
      real(std),dimension(NPRF),intent(in)                 :: lon_o
      logical                  ,intent(in)                 :: periodic_lon
      ! OUTPUT VARIABLES:
      !  Obs that are not in the model domain
      logical,dimension(NPRF),intent(out)             :: discard_obs
      !  The horizontal interp matrix in the CSR representation (structure with 3 arrays):
      type(CSR_FORMAT)                                :: H_csr

      ! LOCAL VARIABLES:
      !  variables used for horizontal interpolations:
      !  the 4 model pixels encircling the obs (NPRF,4) :
      integer,   pointer         :: nhxy(:,:)
      !  weights for the 4 pixels, used in horiz. interpol. (NPRF,4) :
      real(std),    pointer         :: hxy(:,:)
      !
      !  variables used for matrix representation:
      integer*4 nmax
      integer*4 i1,i2,j1,j2
      real(std)    wj1,wj2,wi1,wi2
      !  indexing variables:
      integer*4 :: i,n,k,it
      integer*4 :: icounter
      integer*4 :: ind1,ind2,ind3,ind4,indk
      !
      integer*4 :: ALLOC_ERR
      integer*4,parameter  :: NUMB_NON_ZERO=4  !number of non-zero elements
      !
      !-------------------------------------------
      !-------------------------------------------
      ! the code stars here:
      !-------------------------------------------
      !-------------------------------------------
      !
      allocate(hxy (NPRF,4),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_horiz_interp_matrix")
      allocate(nhxy(NPRF,4),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_horiz_interp_matrix")

      ! calculate the weights:
      call define_hxy(nlon,nlat,NPRF        &
         ,latmod,lonmod           &
         ,lat_o                   &
         ,lon_o                   &
         ,periodic_lon            &
         ,hxy,nhxy,discard_obs)
      !      print *,'define_hxy called'
      !
      ! define the dimensions of the CRS arrays
      H_csr%NPRF  = NPRF
      H_csr%nrow  = nlev*ntime
      H_csr%ncol  = nlat*nlon*nlev*ntime
      nmax        = NUMB_NON_ZERO*H_csr%nrow  !there are 4 non 0 elements in each row
      H_csr%nmax  = nmax

      allocate(H_csr%H(nmax,NPRF),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_horiz_interp_matrix")
      allocate(H_csr%j(nmax,NPRF),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_horiz_interp_matrix")
      allocate(H_csr%i(H_csr%nrow+1,NPRF),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_horiz_interp_matrix")


      do n=1,NPRF
         if (discard_obs(n)) then
            write(*,*) 'WARNING, found observations that are not in the model domain!'
            write(*,*) 'observation number ',n,'lat,lon',lat_o(n),lon_o(n)
            H_csr%H(:,n)=0.
            H_csr%i(:,n)=0
            H_csr%j(:,n)=0
         else
            ! construct the matrix for the horizontal interpolation in the
            ! CSR format
            !
            ! positions of the 4 model pixels:
            j1  = nhxy(n,1)   !lat1
            j2  = nhxy(n,2)   !lat2
            i1  = nhxy(n,3)   !lon1
            i2  = nhxy(n,4)   !lon2

            if(i1.lt.0.or.i2.lt.0.or.j1.lt.0.or.j2.lt.0) cycle

            ! weights for the 4 model pixels (in the j direction, and in i direction):
            wj1 = hxy (n,1)
            wj2 = hxy (n,2)
            wi1 = hxy (n,3)
            wi2 = hxy (n,4)
            !
            ind1=(j1-1)*nlon+i1
            ind2=(j1-1)*nlon+i2
            ind3=(j2-1)*nlon+i1
            ind4=(j2-1)*nlon+i2

            ! construct the operator matrix: (for explanation see the
            ! corresponding document)
            ! the dimensions of the full matrix are :
            ! Rows=ntime*nlev ,  Columns=nlat*nlon*nlev*ntime
            !
            icounter=1
            do it=1,ntime
               do k=1,nlev
                  ! assign the real array with non-zero weights
                  H_csr%H(icounter,  n)=wi1*wj1
                  H_csr%H(icounter+1,n)=wj1*wi2
                  H_csr%H(icounter+2,n)=wj2*wi1
                  H_csr%H(icounter+3,n)=wj2*wi2
                  !assign the integer array with the column positions for non-zero weights
                  indk=(it-1)*nlat*nlon*nlev+(k-1)*nlat*nlon
                  H_csr%j(icounter  ,n)=indk+ind1
                  H_csr%j(icounter+1,n)=indk+ind2
                  H_csr%j(icounter+2,n)=indk+ind3
                  H_csr%j(icounter+3,n)=indk+ind4
                  !
                  icounter=icounter+NUMB_NON_ZERO
                  !
               end do
               !assign the integer array with the positions of the first non-zero element
               ! in each row of the full matrix
               do i=1,nlev
                  H_csr%i((it-1)*nlev+i,n)=(it-1)*nlev*NUMB_NON_ZERO+  &
                     (i-1)*NUMB_NON_ZERO+1
               end do
            end do
            H_csr%i(ntime*nlev+1,n)=H_csr%i(1,n)+H_csr%nmax

            ! matrix construction finished
            !
         end if !if discard_obs
      end do  !n=1,NPRF
      !
      deallocate(hxy)
      deallocate(nhxy)

      return
   end subroutine transf_horiz_interp_matrix

   !******************************************************************
   !******************************************************************
   ! THE ROUTINE transf_vert_interp_matrix IS RESPONSIBLE FOR THE PREPARING THE
   ! MATRIX CORRESPONDING TO THE VERTICAL INTERPOLATION OF THE
   ! MODEL DATA TOWARDS THE VERTICAL (pressure) GRID OF THE OBSERVATIONS.
   ! THE INPUT DATA IS GIVEN FOR ONE or TWO MOCAGE TIME PERIODS.
   !******************************************************************
   !******************************************************************
   !******************************************************************
   subroutine transf_vert_interp_matrix(pres_inth,                    &
      pres_o,                       &
      discard_obs,                  &
      NPRF,nlev,mlev,ntime,         &
      H_csr)

      use floats
      use MathLib
      use interp_matrix_structure_mod
      use matrix_manipulations

      implicit none

      ! INPUT VARIABLES
      integer*4,intent(in)       :: nlev  !model levels
      integer*4,intent(in)       :: mlev  !obs levels
      integer*4,intent(in)       :: NPRF
      integer*4,intent(in)       :: ntime

      real(std), dimension(ntime*nlev,NPRF),intent(in)     :: pres_inth   !model pressures
      real(std), dimension(mlev,NPRF)      ,intent(in)     :: pres_o      !obs pressures
      logical,dimension(NPRF)              ,intent(in)     :: discard_obs !obs to discard
      !
      ! OUTPUT VARIABLES
      !  The vertical interpolation matrix in the CSR representation (structure with 3 arrays):
      type(CSR_FORMAT)    :: H_csr  !REAL ARRAY

      !  information on the size of the arrays in the CSR format:

      !
      ! LOCAL VARIABLES
      ! VARIABLES FOR THE CONSTRUCTION OF THE MATRIX
      integer                          :: nmax
      !
      real(std), dimension(nlev)       :: unit                !unit vector
      integer*4                        :: ind1,ind2,ind3,ind4 !indices
      real(std),dimension(mlev)        :: temp                !temporary array for storing
      !the contribution of a model level
      !to all observed levels

      integer*4                  :: k
      integer*4                  :: ip
      integer*4                  :: it
      integer*4                  :: i_nlev,i_mlev
      integer*4                  :: ierr
      integer*4 :: ALLOC_ERR

      real(std), pointer    ,dimension(:)   :: weights
      integer,pointer    ,dimension(:)   :: positions

      integer*4,parameter  :: NUMB_NON_ZERO=2  !number of non-zero elements in each row
      ! construct the matrix:
      !
      ! define the dimensions of the CRS arrays
      H_csr%NPRF  = NPRF
      H_csr%nrow  = mlev*ntime
      H_csr%ncol  = nlev*ntime
      nmax        = NUMB_NON_ZERO*H_csr%nrow      !there are 2 non 0 elements in each row
      H_csr%nmax=nmax
      !
      allocate(H_csr%H(nmax,NPRF),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_vert_interp_matrix")
      allocate(H_csr%j(nmax,NPRF),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_vert_interp_matrix")
      allocate(H_csr%i(H_csr%nrow+1,NPRF),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_vert_interp_matrix")
      !
      allocate(weights(nmax),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_vert_interp_matrix")
      allocate(positions(nmax),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_vert_interp_matrix")
      !
      do ip=1,NPRF
         if (.not.discard_obs(ip)) then
            do it=1,ntime
               ind1=(it-1)*mlev*NUMB_NON_ZERO
               ind2=it*mlev*NUMB_NON_ZERO
               ind3=(it-1)*nlev
               ind4=it*nlev
               ! get the weights and the indexes of the model levels used for interpolating
               ! to pressure levels:
               call get_coef_interpol(pres_inth(ind3+1:ind4,ip),     &
                  pres_o(:,ip),                                   &
                  weights(ind1+1:ind2),positions(ind1+1:ind2))
            end do
            !construct the matrix:
            ! fill the upper half of the matrix:
            do k=1,H_csr%nrow
               H_csr%H(k,ip)=weights(k)
               H_csr%j(k,ip)=positions(k)
            end do
            ! fill the lower half of the matrix:
            ! (if ntime=2, shift the columns by nlev, see the document)
            do k=H_csr%nrow+1,NUMB_NON_ZERO*H_csr%nrow
               H_csr%H(k,ip)=weights(k)
               H_csr%j(k,ip)=nlev*(ntime-1)+positions(k)
            end do
            do k=1,H_csr%nrow+1
               H_csr%i(k,ip)=(k-1)*NUMB_NON_ZERO+1
            end do
         end if
      end do
      deallocate(weights)
      deallocate(positions)

      return
   end subroutine transf_vert_interp_matrix

   !******************************************************************
   !******************************************************************
   ! THE ROUTINE transf_temp_interp_matrix IS RESPONSIBLE FOR THE PREPARING
   ! THE MATRIX FOR TEMPORAL INTERPOLATION OF THE MODEL DATA TO THE TIME OF
   ! OBS USING THE TWO MOCAGE INTERVALS at T and T+DT.
   !******************************************************************
   !******************************************************************
   !******************************************************************
   subroutine transf_temp_interp_matrix(mlev,NPRF,ntime,date_in,date_out, &
      H_csr  )
      !
      use floats
      use interp_matrix_structure_mod
      use matrix_manipulations
      !
      implicit none

      ! input variables:
      integer,intent(in)                       :: mlev   !obs levels
      integer,intent(in)                       :: ntime  !2 model times
      integer,intent(in)                       :: NPRF   !# of profiles
      !
      real(std), dimension(ntime),intent(in)   ::date_in !the two dates of the model
      real(std), dimension(NPRF),intent(in)    ::date_out !the dates of the obs

      ! output variable
      !  The vertical interp matrix in the CSR representation (structure with 3 ARRAYS):
      type(CSR_FORMAT),intent(out)                     :: H_csr  !REAL ARRAY


      ! local variables:
      real(std)                     :: fact
      real(std)                     :: w1,w2   !weights
      integer                       :: i,n,k
      integer                       :: nmax
      integer*4 :: ALLOC_ERR

      integer*4,parameter  :: NUMB_NON_ZERO=2  !number of non-zero elements


      !-------------------------------------------
      ! the code stars here:
      !-------------------------------------------
      ! construct the sparse matrix in the CRS format:

      ! define the dimensions of the CRS arrays
      H_csr%NPRF  = NPRF
      H_csr%nrow  = mlev
      H_csr%ncol  = mlev*ntime
      nmax        = NUMB_NON_ZERO*H_csr%nrow    !there are 2 non 0 elements in each row
      H_csr%nmax=nmax

      allocate(H_csr%H(nmax,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_temp_interp_matrix")
      allocate(H_csr%j(nmax,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_temp_interp_matrix")
      allocate(H_csr%i(H_csr%nrow+1,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_temp_interp_matrix")

      if (ntime.ne.2) then
         write(*,*) 'ERROR, the time interpolation routine was called for ntime which is not 2'
         call exit(-1)
      end if
      fact=1./(date_in(2)-date_in(1))
      do i=1,NPRF
         w1=(date_in(2)-date_out(i))*fact
         w2=(date_out(i)-date_in(1))*fact
         ! create matrix:
         !   first create the real array with non zero values: (see the corresponding document)
         do k=1,mlev
            H_csr%H((k-1)*2+1,i)=w1
            H_csr%H((k-1)*2+2,i)=w2
         end do
         !   next the integer array with the column position of the NZ elements of the matrix
         do k=1,mlev
            H_csr%j((k-1)*2+1,i)=k
            H_csr%j((k-1)*2+2,i)=k+mlev
         end do
         !   and finally the integer array with the pointers to the first element in
         !   each row of the full matrix.
         do k=1,mlev+1
            H_csr%i(k,i)=(k-1)*NUMB_NON_ZERO+1
         end do
      end do
      !
      return
   end subroutine transf_temp_interp_matrix

   ! DEVELOPPED AT NCAR
   !
   subroutine define_hxy(nlon,nlat,nobs &
      ,lat_model,lon_model,lat_obs,lon_obs,periodic_domain &
      ,weight,position,discard_obs)
      !
      use floats
      implicit none
      !
      integer                        , intent(in   ) :: nlon,nlat,nobs
      real(std), dimension(nlat)     , intent(in   ) :: lat_model
      real(std), dimension(nlon)     , intent(in   ) :: lon_model
      real(std), dimension(nobs)     , intent(in   ) :: lat_obs,lon_obs
      logical                        , intent(in   ) :: periodic_domain
      logical, dimension(nobs)  , intent(inout)      :: discard_obs
      integer, dimension(nobs,4), intent(inout)      :: position
      real(std), dimension(nobs,4), intent(inout)    :: weight
      !
      ! local
      !
      integer      :: i,j,k,n,i1,j1,i2,j2,ii,jj,ip
      real(std)    :: dx1,dx2,dy1,dy2
      logical      :: test_lat,test_lon,found_cell
      logical      :: ll_not_mono
      !
      ! check id domain is periodic in the longitudinal directions
      !
      !  call check_periodic_domain(nlon,nlat,lon_model,periodic_domain)
      !  print *,'periodic_domain = ',periodic_domain
      !  periodic_domain=.true.
      !
      ! reset output arrays
      !
      discard_obs = .true.
      position    = -99
      weight      = 0.
      !
      ! loop over longitude checking monotonicity
      !
      LON: do ip=1,nlon-1
         ll_not_mono = (lon_model(ip).gt.lon_model(ip+1))
         if (ll_not_mono) exit LON
      end do LON
      !
      ! loop over observations
      !
!$omp parallel default(none), &
!$omp shared(nobs,nlon,nlat,ip), &
!$omp shared(discard_obs,position,weight), &
!$omp shared(periodic_domain,ll_not_mono), &
!$omp shared(lon_obs,lon_model,lat_obs,lat_model), &
!$omp private(k,found_cell,test_lon,test_lat), &
!$omp private(i1,i2,j1,j2,i,j,ii,jj,dx1,dx2,dy1,dy2)
!$omp do
      OBS : do k=1,nobs
         !
         ! reset flag and cell coordinates
         !
         found_cell = .false.
         i1         = -99
         i2         = -99
         j1         = -99
         j2         = -99
         !aak modified to speed up, assumption of a regular grid domain...
         LOOP_I : do i=1,nlon-1

            test_lon = .false.
            test_lon = (lon_obs(k).ge.lon_model(i  )).and. &
               (lon_obs(k).le.lon_model(i+1))
            if (test_lon) exit LOOP_I
         end do LOOP_I
         ii=i
         LOOP_J : do j=1,nlat-1
            test_lat = .false.
            test_lat = (lat_obs(k).ge.lat_model(j  )).and. &
               (lat_obs(k).le.lat_model(j+1))
            if (test_lat) exit LOOP_J
         end do LOOP_J
         jj=j
         !
         ! easy case : found a cell
         !
         if(test_lat.and.test_lon) then
            i1 = ii
            i2 = ii + 1
            j1 = jj
            j2 = jj + 1
            found_cell = .true.
         endif
         !
         ! latitude found, but not longitude => check if periodic
         !
         if(test_lat.and.(.not.test_lon)) then
            j1 = j
            j2 = j + 1
            if(periodic_domain) then
               if(ll_not_mono) then
                  if ((lon_obs(k).ge.lon_model(ip).and.&
                     lon_obs(k).lt.lon_model(ip+1)+360.).or. &
                     (lon_obs(k).ge.lon_model(ip)-360..and.&
                     lon_obs(k).lt.lon_model(ip+1))) then
                     i1 = ip
                     i2 = ip+1
                     found_cell = .true.
                  elseif (lon_obs(k).ge.lon_model(nlon).and.&
                     lon_obs(k).lt.lon_model(1)) then
                     i1 = nlon
                     i2 =    1
                     found_cell = .true.
                  endif
               else
                  test_lon = (lon_obs(k).lt.lon_model(1)) .or. &
                     (lon_obs(k).ge.lon_model(nlon))
                  if(test_lon) then
                     i1 = nlon
                     i2 =    1
                     found_cell = .true.
                  endif
               endif
            endif
         endif

         !
         ! if no cell was found, that means that the observation
         ! is outside the model domain, so discard it
         !
         if(.not.found_cell) then
            discard_obs(k) = .true.
            cycle OBS
         endif
         !
         ! if you are here, that means that a cell was found.
         ! One can therefore calculate the weights for each corner
         ! and make sure that the observation is not discarded.
         !
         discard_obs(k) = .false.
         position(k,1) = j1
         position(k,2) = j2
         position(k,3) = i1
         position(k,4) = i2
         if(ll_not_mono) then
            if(i1.eq.ip) then
               dx1 = abs(lon_obs(k) - lon_model(i1))
               dx2 = abs(lon_obs(k) - lon_model(i2) - 360. )
            else
               dx1 = abs(lon_obs(k) - lon_model(i1))
               dx2 = abs(lon_obs(k) - lon_model(i2))
            endif
         else
            if(i1.lt.i2) then
               dx1 = abs(lon_obs(k) - lon_model(i1))
               dx2 = abs(lon_obs(k) - lon_model(i2))
            else
               dx1 = abs(lon_obs(k) - lon_model(i1))
               dx2 = abs(lon_obs(k) - lon_model(i2) - 360. )
            endif
         endif
         dy1 = abs(lat_obs(k) - lat_model(j1))
         dy2 = abs(lat_obs(k) - lat_model(j2))
         if (dy1.lt.1e-4) then
            weight(k,1) = 1.
            weight(k,2) = 0.
         elseif(dy2.lt.1e-4) then
            weight(k,1) = 0.
            weight(k,2) = 1.
         else
            weight(k,1) = dy2/(dy1+dy2)
            weight(k,2) = dy1/(dy1+dy2)
         end if
         if (dx1.lt.1e-4) then
            weight(k,3) = 1.
            weight(k,4) = 0.
         elseif(dx2.lt.1e-4) then
            weight(k,3) = 0.
            weight(k,4) = 1.
         else
            weight(k,3) = dx2/(dx1+dx2)
            weight(k,4) = dx1/(dx1+dx2)
         end if
         !
         !   if(mod(k,25).eq.0) then
         !     print *,k,lon_obs(k),lon_model(j1),lon_model(i2),dx1,dx2
         !     print *,k,lat_obs(k),lat_model(j1),lat_model(i1),dy1,dy2
         !   endif
         !
      end do OBS
!$omp end do
!$omp end parallel
      !
      return

   end subroutine define_hxy

end module transformations_matrix









module space_time_operator_mod

   ! Modules that only groups former single subroutines for the
   ! construction of noveltis spatiotemporal interpolator

   private dnscsr

contains
   !
   !******************************************************************
   !******************************************************************
   ! VERSION: 1.1
   ! WRITTEN BY: NOVELTIS (A.KLONECKI)
   ! THE PROGRAM IS BASED ON THE ALGORITHMS PRESENT IN THE ASSIMILATION CODE
   ! WRITTEN BY JF Lamarque
   !
   ! CREATION DATE: 10/08/2005
   ! REVISION DATE: 12/01/2006
   !    -added flags : flag_averaging_kernel,flag_temp_interp
   !    -handles non-constant observed pressure levels
   !    -instead of estimating the number of nonzero elements for the products
   !     of multiplication of two matrices (nmax), these values are computed
   !     based on analytic expressions.

   !******************************************************************
   subroutine observ_operator  (                 &
      ! INPUT MODEL DATA:
      nlat            & !dimensions:
      ,nlon            &
      ,nlev_m          &
      ,ntime           &
      ,latmod          & !data:
      ,lonmod          &
      ,pres_m          &
      ,date_MOCAGE     &
      ! INPUT OBS DATA:
      ,NPRF            & !dimensions:
      ,n_lvls          &
      ,n_lvls_out      &
      ,lat             & !data:
      ,lon             &
      ,periodic_lon    &
      ,date            &
      ,pres_o          &
      ! EXECUTION OPTIONS:
      ,vert_interp_kind  &
      ,flag_scale_input  &
      ,flag_averaging_kernel &
      ,flag_temp_interp      &
      ! OUTPUT DATA:
      ,Hres_csr        & !Observational Op. in CSR format
      ! INPUT COMPLEMENTARY DATA:
      ,scale_fact      &
      ,avg_kernel      &
      ) !
      !
      use floats
      use interp_matrix_structure_mod
      use transformations_matrix
      use matrix_manipulations
      implicit none
      !-------------------------------------------------------------------
      ! DECLARATIONS:
      !-------------------------------------------------------------------
      ! INPUT DATA:
      !numb. of model longitudes and latitudes
      integer*4,intent(in)                      :: nlon,nlat
      !numb. of model veritcal levels
      integer*4,intent(in)                      :: nlev_m
      !numb. of temporary MOCAGE files: if 1 no temporal interpolation
      !                                 if 2 - temporal interpolation :
      integer*4,intent(in)                      :: ntime

      ! model lats and lons:
      real(std), dimension(nlat),intent(in):: latmod
      real(std), dimension(nlon),intent(in):: lonmod
      ! model fields:
      real(std), dimension(:),intent(in) :: pres_m
      real(std), dimension(ntime),intent(in)    :: date_MOCAGE

      !numb. of obser. profiles in MOCAGE interval being treated:
      integer,intent(in)                        :: NPRF
      !maximum numb. of obs vertical levels
      integer*4,intent(in)                      :: n_lvls
      !this variable dimensions the AvKernel(n_lvls_out,n_lvls),=1 in case
      ! of integrated column data
      integer*4,intent(in)                      :: n_lvls_out

      real(std),dimension(NPRF),intent(in)      :: lat      !IASI lat
      real(std),dimension(NPRF),intent(in)      :: lon      !IASI lon
      logical,intent(in)                        :: periodic_lon
      real(std),dimension(NPRF),intent(in)      :: date     !IASI dates

      real(std),dimension(:,:),intent(in):: pres_o

      integer,intent(in)                      :: vert_interp_kind
      logical,intent(in)                      :: flag_temp_interp
      logical,intent(in)                      :: flag_averaging_kernel
      logical,intent(in)                      :: flag_scale_input

      !-------------------------------------------------------------------
      ! OUTPUT DATA:
      type(CSR_FORMAT),intent(out)              :: Hres_csr

      real(std), dimension(:,:,:),intent(in), optional :: avg_kernel
      real(std), dimension(:),intent(in), optional :: scale_fact

      !-------------------------------------------------------------------
      ! LOCAL VARIABLES:
      integer, pointer, dimension(:)        :: iw  ! work array for the routine amub
      ! THE FOUR MATRICES CORRESPONDING TO THE 4 OPERATORS, MATRICES IN CSR FORMAT:
      type(CSR_FORMAT)                          :: Hh_csr  !horizontal operator
      type(CSR_FORMAT)                          :: Hv_csr  !vertical operator
      type(CSR_FORMAT)                          :: Ht_csr  !temporal operator
      type(CSR_FORMAT)                          :: Ha_csr  !smoothing operator
      ! TEMPORARY MATRICES FOR STORING THE RESULTS OF MULTIPLICATIONS OF MATRICES:
      type(CSR_FORMAT)                          :: Hres1_csr  !result of Hv*Hh
      type(CSR_FORMAT)                          :: Hres2_csr  !result of Ht*(Hv*Hh)

      logical, dimension(NPRF)      :: discard_obs

      ! MODEL PRESSURE INTEGRATED TO THE LOCATION OF OBSERVATION
      real(std),pointer,dimension(:,:)      :: pres_inth

      !# of model levels used in vertical interpolation, if the averaging kernel is used,
      ! this number is used to store the  number of non-zero elements in the final
      !  matrix
      integer*4                                 :: np   !indexing variable
      ! number of non-zero elements in the CSR format:
      integer*4                                 :: nmax
      integer*4                                 :: ierr !error index
      ! actual number of levels for interpolation (depend on interp. kind)
      integer                                   :: nlev

      integer*4                                 :: ALLOC_ERR

      ! loop counters
      integer*4 :: ip, ib, k

      allocate(Hres1_csr%H(1,1))
      allocate(Hres1_csr%j(1,1))
      allocate(Hres1_csr%i(1,1))
      allocate(Hres2_csr%H(1,1))
      allocate(Hres2_csr%j(1,1))
      allocate(Hres2_csr%i(1,1))

      !******************************************************************
      ! END OF DECLARATIONS

      !-----------------------------------------------------------------
      ! 1) PREPARE THE MATRIX TO INTERPOLATE MODEL CO AND PRES FROM MODEL GRID TO
      !    OBS LOCATIONS
      if (vert_interp_kind == ip_2d) then
         nlev = 1
      else
         nlev = nlev_m
      end if

      call transf_horiz_interp_matrix(NPRF,nlon,nlat,nlev,ntime,         &
         latmod,lonmod,                        &
         lat,lon,  &
         periodic_lon, &
         discard_obs,             &
         Hh_csr)

      !-----------------------------------------------------------------
      ! 2) PREPARE THE MATRIX TO INTERPOLATE TO OBS VERTICAL LEVELS (coarse):
      !
      select case(vert_interp_kind)
      case (ip_interp)

         ! need the interpolated pressure levels for the vertical interpolation:
         allocate(pres_inth(Hh_csr%nrow,NPRF),STAT = ALLOC_ERR)
         call check_allocate(ALLOC_ERR,"multiply_matrix_csr_vector")
         call multiply_matrix_csr_vector(Hh_csr,pres_m,          &
            Hh_csr%nrow,NPRF,pres_inth)

         ! pres_inth now contains the model pressures on the obs horizontal grid
         call transf_vert_interp_matrix(pres_inth,pres_o,    &
            discard_obs,     &
            NPRF,nlev_m,n_lvls,ntime,       &
            Hv_csr)

         deallocate(pres_inth)
      case(ip_ground)

         Hv_csr%NPRF = NPRF
         Hv_csr%nrow = ntime
         Hv_csr%ncol = nlev_m*ntime
         nmax        = Hv_csr%nrow      !there is just 1 non 0 el in each row
         Hv_csr%nmax = nmax
         !
         allocate(Hv_csr%H(nmax,NPRF),STAT=ALLOC_ERR)
         call check_allocate(ALLOC_ERR,"transf_vert_interp_matrix")
         Hv_csr%H(:,:) = 0.
         allocate(Hv_csr%j(nmax,NPRF),STAT=ALLOC_ERR)
         call check_allocate(ALLOC_ERR,"transf_vert_interp_matrix")
         Hv_csr%j(:,:) = 0
         allocate(Hv_csr%i(Hv_csr%nrow+1,NPRF),STAT=ALLOC_ERR)
         call check_allocate(ALLOC_ERR,"transf_vert_interp_matrix")
         Hv_csr%i(:,:) = 0

         do ip=1,NPRF
            if (.not.discard_obs(ip)) then
               do k = 1,Hv_csr%nrow
                  Hv_csr%H(k,ip) = 1.
                  Hv_csr%j(k,ip) = k*nlev_m
               end do
               do k = 1,Hv_csr%nrow+1
                  Hv_csr%i(k,ip) = k
               end do
            end if
         end do
      end select
      !-----------------------------------------------------------------
      ! 3) PREPARE THE MATRIX TO INTERPOLATE WITH RESPECT TO TIME

      if (flag_temp_interp) then   !if the temporal interpolation is
         ! necessary:

         call transf_temp_interp_matrix(n_lvls,                 &
            NPRF,ntime,                              &
            date_MOCAGE(1:2),                        &
            date,                     &
            Ht_csr)
      end if
      !-------------------------------------------
      ! 4) the Averaging matrix is already calculated, convert to the CSR format
      if (flag_averaging_kernel) then
         Ha_csr%NPRF = NPRF
         Ha_csr%ncol = n_lvls
         Ha_csr%nrow = n_lvls_out
         nmax=Ha_csr%ncol*Ha_csr%nrow
         Ha_csr%nmax = nmax
         allocate(Ha_csr%H(nmax,NPRF))    !for now the same for all profiles
         allocate(Ha_csr%j(nmax,NPRF))
         allocate(Ha_csr%i(Ha_csr%nrow+1,NPRF))
         do np=1,NPRF
            call dnscsr(Ha_csr%nrow,Ha_csr%ncol,nmax,  &
               avg_kernel(:,:,np),Ha_csr%nrow,Ha_csr%H(:,np),  &
               Ha_csr%j(:,np),Ha_csr%i(:,np),ierr)
         end do
      end if

      !---
      !-------------------------------------------
      !-------------------------------------------
      ! FINISHED CONSTRUCTING THE MATRICES, NOW MULTIPLY THEM TO OBTAIN THE FINAL
      ! OBSERVATIONAL OPERATOR H=Ha*Ht*Hv*Hh
      ! FOR THE INDIVIDUAL MATRICES IT WAS STRAIGHTFORWARD TO FIND THE NUMBER OF
      ! NON ZERO ELEMENTS. FOR THE PRODUCTS OF TWO MATRICES, THE NUMBER OF NON ZERO
      ! ELEMENTS CAN BE ALSO ESTIMATED ANALYTICALLY. THE FORMULAS ARE GIVEN
      ! BELOW.
      !-------------------------------------------
      !-------------------------------------------
      ! the multplication of Hv*Hh gives 8 non zero elements in each row.
      select case(vert_interp_kind)
      case (ip_interp, ip_ground)
         nmax = 8*Hv_csr%nrow

         allocate(iw(Hh_csr%ncol))           !work array needed for routine amub
         iw=0

         call multiply_matrices_csr(iw &
            ,Hv_csr    & !matrix A
            ,Hh_csr    & !matrix B
            ,Hres1_csr & !out:A*B
            ,nmax)

         !
         call deallocate_operator(Hv_csr)
         deallocate(iw)
      case (ip_2d, ip_identity)
         call copy_matrix_csr(Hh_csr,Hres1_csr)
      end select

      call deallocate_operator(Hh_csr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !next calculate: H_t*H_res where H_res is H_v*H_h:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flag_temp_interp) then
         !the multiplication of Ht*Hres gives 16 non zero values in each row
         nmax = 16 * Ht_csr%nrow

         allocate(iw(Hres1_csr%ncol))    !work array needed for routine amub
         iw=0
         call multiply_matrices_csr(iw     &
            ,Ht_csr     &    !matrix A
            ,Hres1_csr  &    !matrix B
            ,Hres2_csr  &    !out:A*B
            ,nmax)
         call deallocate_operator(Ht_csr)
         deallocate(iw)
      else
         call copy_matrix_csr(Hres1_csr,Hres2_csr)
      end if  !if flag_temp_interp

      call deallocate_operator(Hres1_csr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( flag_averaging_kernel ) then
         ! calculate H_a*H_res where H_res is : H_t*H_v*H_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nmax = get_number_col(Hres2_csr,NPRF,discard_obs)*Ha_csr%nrow
         allocate(iw(Hres2_csr%ncol))    !work array needed for routine amub
         iw=0
         call multiply_matrices_csr(iw      &
            ,Ha_csr        &    !matrix A
            ,Hres2_csr     &    !matrix B
            ,Hres_csr      &    !out:A*B
            ,nmax)
         call deallocate_operator(Ha_csr)
         deallocate(iw)
      else
         call copy_matrix_csr(Hres2_csr,Hres_csr)
      end if

      call deallocate_operator(Hres2_csr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( flag_scale_input ) then
         ! and finaly scale the input by multiplication of the weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         do ip = 1, NPRF
            do ib = 1,Hres_csr%nrow
               do k = Hres_csr%i(ib,ip), Hres_csr%i(ib+1,ip)-1
                  Hres_csr%H(k,ip) = &
                     & Hres_csr%H(k,ip) * scale_fact(Hres_csr%j(k,ip))
               end do
            end do
         end do
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! finished calculating the H matrix,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
   end subroutine observ_operator

   subroutine observ_operator_tc  (                 &
      ! INPUT MODEL DATA:
      nlat            & !dimensions:
      ,nlon            &
      ,nlev_m          &
      ,ntime           &
      ,latmod          & !data:
      ,lonmod          &
      ,pres_m          &
      ,pres_m_int      &
      ,date_MOCAGE     &
      ! INPUT OBS DATA:
      ,NPRF            & !dimensions:
      ,n_lvls_out      & !output dim (after avg kern)
      ,lat             & !data:
      ,lon             &
      ,periodic_lon    &
      ,date            &
      ,pres_bounds     &
      ! EXECUTION OPTIONS:
      ,flag_averaging_kernel &
      ,flag_temp_interp      &
      ! OUTPUT DATA:
      ,Hres_csr        & !Observational Op. in CSR format
      ,avg_kernel      & !Optional averaging kernel
      ) !
      !
      use floats
      use interp_matrix_structure_mod
      use transformations_matrix
      use matrix_manipulations
      implicit none
      !-------------------------------------------
      ! DECLARATIONS:
      !-------------------------------------------
      ! INPUT DATA:
      !numb. of model longitudes and latitudes
      integer*4,intent(in)                      :: nlon,nlat
      !numb. of model veritcal levels
      integer*4,intent(in)                      :: nlev_m
      !numb. of temporary MOCAGE files: if 1 no temporal interpolation
      !                                 if 2 - temporal interpolation :
      integer*4,intent(in)                      :: ntime

      ! model lats and lons:
      real(std), dimension(nlat),intent(in):: latmod
      real(std), dimension(nlon),intent(in):: lonmod
      ! model fields:
      real(std), dimension(nlon*nlat*nlev_m*ntime),intent(in) :: pres_m
      real(std), dimension(nlon*nlat*nlev_m*ntime),intent(in) :: pres_m_int
      real(std), dimension(ntime),intent(in)    :: date_MOCAGE

      !numb. of obser. profiles in MOCAGE interval being treated:
      integer,intent(in)                        :: NPRF
      !this variable dimensions the AvKernel(n_lvls_out,n_lvls),=1 in case
      ! of integrated column data
      integer*4,intent(in)                      :: n_lvls_out

      real(std),dimension(NPRF),intent(in)      :: lat      !IASI lat
      real(std),dimension(NPRF),intent(in)      :: lon      !IASI lon
      logical,intent(in)                        :: periodic_lon
      real(std),dimension(NPRF),intent(in)      :: date     !IASI dates
      real(std),dimension(:,:,:),intent(in) :: pres_bounds !col. bounds (2,n_lvls,NPRF)

      logical,intent(in)                      :: flag_averaging_kernel
      logical,intent(in)                      :: flag_temp_interp

      real(std), dimension(:,:,:),intent(in), optional :: avg_kernel ! (n_lvls_out,n_lvls,NPRF)

      !-------------------------------------------
      ! OUTPUT DATA:
      type(CSR_FORMAT),intent(out)              :: Hres_csr

      !-------------------------------------------
      ! LOCAL VARIABLES:
      integer, pointer, dimension(:)        :: iw  ! work array for the routine amub
      ! THE FOUR MATRICES CORRESPONDING TO THE 4 OPERATORS, MATRICES IN CSR FORMAT:
      type(CSR_FORMAT)                          :: Hh_csr  !horizontal operator
      type(CSR_FORMAT)                          :: Ht_csr  !temporal operator
      type(CSR_FORMAT)                          :: Hc_csr  !column operator
      type(CSR_FORMAT)                          :: Ha_csr  !smoothing operator
      ! TEMPORARY MATRICES FOR STORING THE RESULTS OF MULTIPLICATIONS OF MATRICES:
      type(CSR_FORMAT)                          :: Hres2_csr  !result of Ht*Hh
      type(CSR_FORMAT)                          :: Hres3_csr  !result of Hc*Ht*Hh

      logical, dimension(NPRF)     :: discard_obs

      ! MODEL PRESSURE INTEGRATED TO THE LOCATION OF OBSERVATION
      real(std),pointer,dimension(:,:)      :: pres_inth
      real(std),pointer,dimension(:,:)      :: pres_inth_int

      real(std), dimension(:,:,:),allocatable  :: col_kernel

      !# of model levels used in vertical interpolation,
      ! if the averaging kernel is used,
      ! this number is used to store the  number of non-zero elements in the final
      !  matrix
      integer*4                                 :: n_lvls
      integer*4                                 :: lvls_out
      integer*4                                 :: np   !indexing variable
      ! number of non-zero elements in the CSR format:
      integer*4                                 :: nmax
      integer*4                                 :: ierr !error index

      integer*4                                 :: ALLOC_ERR


      !******************************************************************
      ! END OF DECLARATIONS

      if (flag_averaging_kernel) then
         n_lvls = size(avg_kernel,dim=2)
         if (n_lvls_out .ne. size(avg_kernel,DIM=1)) write(0,*) ' WRONG AVG KERN 1ST DIM'
         if (NPRF .ne. size(avg_kernel,DIM=3)) write(0,*) ' WRONG AVG_KERN 3RD DIM'
      else
         n_lvls = n_lvls_out
      end if

      if (2 .ne. size(pres_bounds,DIM=1)) write(0,*) ' WRONG PRES BOUNDS 1ST DIM'
      if (n_lvls .ne. size(pres_bounds,DIM=2)) write(0,*) ' WRONG PRES BOUNDS 2ND DIM'
      if (NPRF .ne. size(pres_bounds,DIM=3)) write(0,*) ' WRONG PRES BOUNDS 3RD DIM'

      !-------------------------------------------
      ! 1) PREPARE THE MATRIX TO INTERPOLATE MODEL CO AND PRES FROM MODEL GRID TO
      !    OBS LOCATIONS
      call transf_horiz_interp_matrix(NPRF,nlon,nlat,nlev_m,ntime,      &
         latmod,lonmod,                      &
         lat,lon,                            &
         periodic_lon,                       &
         discard_obs,                        &
         Hh_csr)
      !
      !-------------------------------------------
      ! 2) NO NEED OF VERTICAL INTERPOLATION:
      !
      !-------------------------------------------
      ! 3) PREPARE THE MATRIX TO INTERPOLATE WITH RESPECT TO TIME

      if (flag_temp_interp) then   !if the temporal interpolation is
         ! necessary:

         call transf_temp_interp_matrix(nlev_m,                &
            NPRF,ntime,            &
            date_MOCAGE(1:2),      &
            date,                  &
            Ht_csr)
      end if
      !
      !-------------------------------------------
      !-------------------------------------------
      ! FINISHED CONSTRUCTING THE MATRICES, NOW MULTIPLY THEM TO OBTAIN THE FINAL
      ! OBSERVATIONAL OPERATOR H=Ha*Ht*Hh
      ! FOR THE INDIVIDUAL MATRICES IT WAS STRAIGHTFORWARD TO FIND THE NUMBER OF
      ! NON ZERO ELEMENTS. FOR THE PRODUCTS OF TWO MATRICES, THE NUMBER OF NON ZERO
      ! ELEMENTS CAN BE ALSO ESTIMATED ANALYTICALLY. THE FORMULAS ARE GIVEN
      ! BELOW.
      !-------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculate: H_t*H_h:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flag_temp_interp) then
         !the multiplication of Ht*Hres gives 16 non zero values in each row
         nmax = 16 * Ht_csr%nrow

         allocate(iw(Hh_csr%ncol))    !work array needed for routine amub
         iw=0
         call multiply_matrices_csr(iw,         &
            Ht_csr,     &    !matrix A
            Hh_csr,     &    !matrix B
            Hres2_csr,  &    !out:A*B
            nmax)
         call deallocate_operator(Ht_csr)
         deallocate(iw)
      else
         call copy_matrix_csr(Hh_csr,Hres2_csr)
      end if!if flag_temp_interp
      !
      call deallocate_operator(Hh_csr)
      !
      !-------------------------------------------
      ! need the interpolated pressure levels for the total column:
      allocate(pres_inth(Hres2_csr%nrow,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"multiply_matrix_csr_vector")
      call multiply_matrix_csr_vector(Hres2_csr,pres_m,          &
         Hres2_csr%nrow,NPRF,pres_inth)
      !
      allocate(pres_inth_int(Hres2_csr%nrow,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"multiply_matrix_csr_vector")
      call multiply_matrix_csr_vector(Hres2_csr,pres_m_int,          &
         Hres2_csr%nrow,NPRF,pres_inth_int)
      !
      ! pres_inth now contains the model pressures on the obs horizontal grid
      ! pres_inth_int now contains the model interface press. on the obs hor. grid
      !
      !-------------------------------------------
      ! 4) Compute the avg kernel matrix for total columns (cubic),
      !    then convert it to the CSR format
      !
      allocate(col_kernel(n_lvls,nlev_m,NPRF))
      call transf_totcol_avg_kernel(nlev_m,n_lvls,NPRF,pres_inth, &
         pres_inth_int,pres_bounds,col_kernel)
      !
      deallocate(pres_inth)
      deallocate(pres_inth_int)
      !
      Hc_csr%NPRF = NPRF
      Hc_csr%ncol = nlev_m
      Hc_csr%nrow = n_lvls
      nmax=Hc_csr%ncol*Hc_csr%nrow
      Hc_csr%nmax = nmax
      allocate(Hc_csr%H(nmax,NPRF))    !for now the same for all profiles
      allocate(Hc_csr%j(nmax,NPRF))
      allocate(Hc_csr%i(Hc_csr%nrow+1,NPRF))
      do np=1,NPRF
         call dnscsr(Hc_csr%nrow,Hc_csr%ncol,nmax,  &
            col_kernel(:,:,np),Hc_csr%nrow,Hc_csr%H(:,np),  &
            Hc_csr%j(:,np),Hc_csr%i(:,np),ierr)
      end do
      deallocate(col_kernel)
      !
      !---
      !-------------------------------------------
      ! 5) the Averaging matrix is already calculated, convert to the CSR format

      if (flag_averaging_kernel) then
         Ha_csr%NPRF = NPRF
         Ha_csr%ncol = n_lvls
         Ha_csr%nrow = n_lvls_out
         nmax=Ha_csr%ncol*Ha_csr%nrow
         Ha_csr%nmax = nmax
         allocate(Ha_csr%H(nmax,NPRF))    !for now the same for all profiles
         allocate(Ha_csr%j(nmax,NPRF))
         allocate(Ha_csr%i(Ha_csr%nrow+1,NPRF))
         do np=1,NPRF
            call dnscsr(Ha_csr%nrow,Ha_csr%ncol,nmax,  &
               avg_kernel(:,:,np),Ha_csr%nrow,Ha_csr%H(:,np),  &
               Ha_csr%j(:,np),Ha_csr%i(:,np),ierr)
         end do
      end if

      !---
      !
      ! calculate H_c*H_res where H_res is : H_t*H_h
      !-------------------------------------------
      !
      lvls_out = nlev_m/2+1
      nmax = lvls_out * 16 * Hc_csr%nrow

      allocate(iw(Hres2_csr%ncol))    !work array needed for routine amub
      iw=0
      call multiply_matrices_csr(iw,            &
         Hc_csr,        &    !matrix A
         Hres2_csr,     &    !matrix B
         Hres3_csr,      &    !out:A*B
         nmax)
      call deallocate_operator(Hc_csr)
      deallocate(iw)
      call deallocate_operator(Hres2_csr)
      !---
      !
      ! and finally calculate H_a*H_res where H_res is : H_c*H_t*H_h
      !-------------------------------------------
      !
      if ( flag_averaging_kernel ) then
         nmax = get_number_col(Hres3_csr,NPRF,discard_obs)*Ha_csr%nrow
         allocate(iw(Hres3_csr%ncol))    !work array needed for routine amub
         iw=0
         call multiply_matrices_csr(iw      &
            ,Ha_csr        &    !matrix A
            ,Hres3_csr     &    !matrix B
            ,Hres_csr      &    !out:A*B
            ,nmax)
         call deallocate_operator(Ha_csr)
         deallocate(iw)
      else
         call copy_matrix_csr(Hres3_csr,Hres_csr)
      end if

      call deallocate_operator(Hres3_csr)
      !
      !-------------------------------------------
      !
      ! finished calculating the H matrix,
      !-------------------------------------------
      !
      return

   end subroutine observ_operator_tc

   subroutine observ_operator_vert  (                 &
      ! INPUT MODEL DATA:
      n_prof           & !dimensions:
      ,nlev_m          &
      ,pres_m          &
      ! INPUT OBS DATA:
      ,n_lvls          &
      ,n_lvls_out   &
      ,pres_o          &
      ! EXECUTION OPTIONS:
      ,flag_averaging_kernel &
      ! OUTPUT DATA:
      ,Hres_csr        & !Observational Op. in CSR format
      ! INPUT COMPLEMENTARY DATA:
      ,avg_kernel      &
      ) !
      !
      use floats
      use interp_matrix_structure_mod
      use transformations_matrix
      use matrix_manipulations
      implicit none
      !-------------------------------------------
      ! DECLARATIONS:
      !-------------------------------------------
      ! INPUT DATA:
      !numb. of model/obs profiles
      integer*4,intent(in)                      :: n_prof
      !numb. of model veritcal levels
      integer*4,intent(in)                      :: nlev_m
      ! model fields:
      real(std), dimension(nlev_m,n_prof),intent(in) :: pres_m
      !maximum numb. of obs vertical levels
      integer*4,intent(in)                      :: n_lvls
      !this variable dimensions the AvKernel(n_lvls_out,n_lvls),=1 in case
      ! of integrated column data
      integer*4,intent(in)                      :: n_lvls_out

      real(std),dimension(n_lvls,n_prof),intent(in):: pres_o

      logical,intent(in)                      :: flag_averaging_kernel

      !-------------------------------------------
      ! OUTPUT DATA:
      type(CSR_FORMAT),intent(out)              :: Hres_csr

      real(std), dimension(n_lvls_out,n_lvls,n_prof),intent(in), optional :: avg_kernel

      !-------------------------------------------
      ! LOCAL VARIABLES:
      integer, pointer, dimension(:)        :: iw  ! work array for the routine amub
      ! THE FOUR MATRICES CORRESPONDING TO THE 4 OPERATORS, MATRICES IN CSR FORMAT:
      type(CSR_FORMAT)                          :: Hh_csr  !horizontal operator
      type(CSR_FORMAT)                          :: Hv_csr  !vertical operator
      type(CSR_FORMAT)                          :: Ha_csr  !smoothing operator
      ! TEMPORARY MATRICES FOR STORING THE RESULTS OF MULTIPLICATIONS OF MATRICES:
      type(CSR_FORMAT)                          :: Hres1_csr  !result of Hv*Hh

      logical, dimension(n_prof)      :: discard_obs

      !# of model levels used in vertical interpolation, if the averaging kernel is used,
      ! this number is used to store the  number of non-zero elements in the final
      !  matrix
      integer*4                                 :: np   !indexing variable
      ! number of non-zero elements in the CSR format:
      integer*4                                 :: nmax
      integer*4                                 :: ierr !error index

      integer*4                                 :: ALLOC_ERR

      ! loop counters and indexing variables
      integer*4 :: ip, i, k, n, indk

      allocate(Hres1_csr%H(1,1))
      allocate(Hres1_csr%j(1,1))
      allocate(Hres1_csr%i(1,1))

      !******************************************************************
      ! END OF DECLARATIONS
      ! Build a fake horizontal interpolation matrix, which is indeed not interpolating
      ! anything but only used on the right side of Hv to ensure the correct array
      ! multiplications without need to rewrite specific matrix_multiply routines.
      ! The construction is inspired by transf_horiz_interp_matrix but without search
      ! of locations.
      ! define the dimensions of the CRS arrays
      Hh_csr%NPRF  = n_prof ! The destination number of profiles
      Hh_csr%nrow  = nlev_m ! The destination size (vertical on grid model)
      Hh_csr%ncol  = n_prof*nlev_m ! The input size (the number of input 1D profiles)
      nmax         = Hh_csr%nrow*1  ! One to one correspondence between input profiles and output ones
      Hh_csr%nmax  = nmax

      allocate(Hh_csr%H(nmax,n_prof),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_horiz_interp_matrix")
      Hh_csr%H(:,:) = 0.
      allocate(Hh_csr%j(nmax,n_prof),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_horiz_interp_matrix")
      Hh_csr%j(:,:) = 0
      allocate(Hh_csr%i(Hh_csr%nrow+1,n_prof),STAT=ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"transf_horiz_interp_matrix")
      Hh_csr%i(:,:) = 0

      do n=1,n_prof
         ! construct the operator matrix: (for explanation see the
         ! corresponding document)
         ! the dimensions of the full matrix are :
         ! Rows=nlev ,  Columns=nprof*nlev
         !
         do k=1,nlev_m
            ! assign the real array with non-zero weights
            Hh_csr%H(k,n)=1.
            ! assign the integer array with the column positions for non-zero weights
            indk=(k-1)*n_prof
            Hh_csr%j(k,n)=indk+n
            !
         end do
         !assign the integer array with the positions of the first non-zero element
         ! in each row of the full matrix
         do i=1,nlev_m
            Hh_csr%i(i,n)=i
         end do

         Hh_csr%i(nlev_m+1,n)=Hh_csr%i(1,n)+Hh_csr%nmax

         ! matrix construction finished
      end do  !n=1,n_prof
      !-------------------------------------------
      ! 2) PREPARE THE MATRIX TO INTERPOLATE TO OBS VERTICAL LEVELS (coarse):
      !
      discard_obs(:) = .false.

      call transf_vert_interp_matrix(pres_m,pres_o,    &
         discard_obs,     &
         n_prof,nlev_m,n_lvls,1,       &
         Hv_csr)

      !-------------------------------------------
      ! 4) the Averaging matrix is already calculated, convert to the CSR format
      if (flag_averaging_kernel) then
         Ha_csr%NPRF = n_prof
         Ha_csr%ncol = n_lvls
         Ha_csr%nrow = n_lvls_out
         nmax=Ha_csr%ncol*Ha_csr%nrow
         Ha_csr%nmax = nmax
         allocate(Ha_csr%H(nmax,n_prof))    !for now the same for all profiles
         allocate(Ha_csr%j(nmax,n_prof))
         allocate(Ha_csr%i(Ha_csr%nrow+1,n_prof))
         do np=1,n_prof
            call dnscsr(Ha_csr%nrow,Ha_csr%ncol,nmax,  &
               avg_kernel(:,:,np),Ha_csr%nrow,Ha_csr%H(:,np),  &
               Ha_csr%j(:,np),Ha_csr%i(:,np),ierr)
         end do
      end if
      !---
      !-------------------------------------------
      !-------------------------------------------
      ! FINISHED CONSTRUCTING THE MATRICES, NOW MULTIPLY THEM TO OBTAIN THE FINAL
      ! OBSERVATIONAL OPERATOR H=Ha*Ht*Hv*Hh
      ! FOR THE INDIVIDUAL MATRICES IT WAS STRAIGHTFORWARD TO FIND THE NUMBER OF
      ! NON ZERO ELEMENTS. FOR THE PRODUCTS OF TWO MATRICES, THE NUMBER OF NON ZERO
      ! ELEMENTS CAN BE ALSO ESTIMATED ANALYTICALLY. THE FORMULAS ARE GIVEN
      ! BELOW.
      !-------------------------------------------
      !-------------------------------------------
      ! the multplication of Hv*Hh gives 2 non zero elements in each row.

      nmax = 2*Hv_csr%nrow

      allocate(iw(Hh_csr%ncol))           !work array needed for routine amub
      iw=0

      call multiply_matrices_csr(iw &
         ,Hv_csr    & !matrix A
         ,Hh_csr    & !matrix B
         ,Hres1_csr & !out:A*B
         ,nmax)

      call deallocate_operator(Hv_csr)
      deallocate(iw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( flag_averaging_kernel ) then
         ! and finaly calculate H_a*H_res where H_res is : H_t*H_v*H_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nmax = get_number_col(Hres1_csr,n_prof,discard_obs)*Ha_csr%nrow
         allocate(iw(Hres1_csr%ncol))    !work array needed for routine amub
         iw=0
         call multiply_matrices_csr(iw      &
            ,Ha_csr        &    !matrix A
            ,Hres1_csr     &    !matrix B
            ,Hres_csr      &    !out:A*B
            ,nmax)
         call deallocate_operator(Ha_csr)
         deallocate(iw)
      else
         call copy_matrix_csr(Hres1_csr,Hres_csr)
      end if

      call deallocate_operator(Hres1_csr)
      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! finished calculating the H matrix,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
   end subroutine observ_operator_vert
   !
   !
   subroutine TRANSF_TOTCOL_AVG_KERNEL(NIV,NIV_OBS,NPRF,RP,RPLT,RP_O,AVK)
      !
      !CCCC
      !
      !   VERTICAL INTEGRATION (CUBIC)
      !   __  __   _____   _____   _____   _____   _____
      !  |  \/  | /     \ /     \ /  _  \ /  ___\ |  ___|
      !  |      | |  |  | |  <--< |  _  | |  \_ \ |  ___|
      !  |_|\/|_| \_____/ \_____/ \_/ \_/ \_____/ |_____|
      !
      !   MODELE DE CHIMIE ATMOSPHERIQUE A GRANDE ECHELLE
      !   VERSION_RELEASE_BUGFIX : 0_1_1
      !
      !CC DATE : 09/06/2005
      !
      !   Type      : Chemistry
      !   Purpose   : Compute a total column
      !   Called by :
      !
      !CC
      !
      !  Original : B. Josse & V.-H. Peuch, 04/2000
      !
      !  Modified :
      !
      !CCCC
      !
      implicit none
      !
      !
      !     Subroutine arguments.
      !     ~~~~~~~~~~~~~~~~~~~~
      !
      integer :: NIV, NIV_OBS, NPRF
      real (KIND=8) , dimension (NIV,NPRF) :: RP,RPLT
      real (KIND=8) , dimension (2,NIV_OBS,NPRF) :: RP_O
      real (KIND=8) , dimension (NIV_OBS,NIV,NPRF) :: AVK
      !
      ! LOOPS
      integer :: I,L,K
      ! INDEXES
      integer :: LTOP,LBOT
      ! COEFFICIENTS
      real(KIND=8) :: rl_TOP,rl_BOT
      ! CONSTANTS
      real(KIND=8) :: MAIR,GG,AVO,DOBSON
      ! WEIGHTS FOR NEIGHBORING LEVELS
      real (KIND=8), allocatable   :: ALFZ1MAS(:,:),ALFZ2MAS(:,:)
      real (KIND=8), allocatable   :: ALFZ3MAS(:,:),ALFZ4MAS(:,:)
      real (KIND=8), allocatable   :: ALFZ5MAS(:,:)
      ! WORK ARRAY
      real (KIND=8), dimension (NIV) :: AVK_COEFF
      real (KIND=8), allocatable :: RP_SWAP(:,:)
      real (KIND=8) :: RP_O_TOP, RP_O_BOT
      logical :: LL_SWAP, LL_OUT
      !
      allocate (ALFZ1MAS(NPRF,NIV))
      allocate (ALFZ2MAS(NPRF,NIV))
      allocate (ALFZ3MAS(NPRF,NIV))
      allocate (ALFZ4MAS(NPRF,NIV))
      allocate (ALFZ5MAS(NPRF,NIV))
      !
      ! 1. INITIALIZATIONS
      !
      do I = 1,NPRF
         do L=1,NIV
            do K=1,NIV_OBS
               AVK(K,L,I)=0.
            enddo
         enddo
      enddo
      !
      LL_SWAP = .false.
      if (RP(2,1) .lt. RP(1,1)) then
         allocate(RP_SWAP(NIV,NPRF))
         RP_SWAP(:,:)=RP(:,:)
         do L=1,NIV
            RP(L,:)=RP_SWAP(NIV-L+1,:)
         end do
         !
         RP_SWAP(:,:)=RPLT(:,:)
         do L=1,NIV
            RPLT(L,:)=RP_SWAP(NIV-L+1,:)
         end do
         LL_SWAP = .true.
         deallocate(RP_SWAP)
         allocate(RP_SWAP(NIV_OBS,NIV))
      end if

      !
      MAIR=0.029
      GG=9.80665
      AVO=6.02E23
      DOBSON=2.685E20
      !
      ! 2. COMPUTE WEIGHTS
      !
      call TRANSF_TOTCOL_WEIGHTS(NIV,NPRF,RP,RPLT &
         ,ALFZ1MAS,ALFZ2MAS &
         ,ALFZ3MAS,ALFZ4MAS &
         ,ALFZ5MAS)
      !
      ! 3. COMPUTE NUMBER OF MOLES IN EACH MODEL
      !    CELL (CUBIC INTERPOLATION)
      !
      do I = 1,NPRF
         AVK_COEFF(1)=ALFZ3MAS(I,1)+ALFZ2MAS(I,2)+ALFZ1MAS(I,3)
         AVK_COEFF(2)=ALFZ4MAS(I,1)+ALFZ3MAS(I,2)+ALFZ2MAS(I,3)+ALFZ1MAS(I,4)
         AVK_COEFF(3)=ALFZ4MAS(I,2)+ALFZ3MAS(I,3)+ALFZ2MAS(I,4)+ALFZ1MAS(I,5)
         do L = 4,NIV-2
            AVK_COEFF(L)=ALFZ5MAS(I,L-2)+ALFZ4MAS(I,L-1)+ALFZ3MAS(I,L)&
               +ALFZ2MAS(I,L+1)+ALFZ1MAS(I,L+2)
         end do
         AVK_COEFF(NIV-1)=ALFZ5MAS(I,NIV-3)+ALFZ4MAS(I,NIV-2)&
            +ALFZ3MAS(I,NIV-1)+ALFZ2MAS(I,NIV)
         AVK_COEFF(NIV)=ALFZ5MAS(I,NIV-2)+ALFZ4MAS(I,NIV-1)&
            +ALFZ3MAS(I,NIV)
         !
         AVK_COEFF(:) = AVK_COEFF(:)/GG/MAIR*AVO/DOBSON
         !
         ! 4. DETERMINE LEVELS/COLUMNS INTERSECTIONS
         !
         do K=1,NIV_OBS ! For each subcolumn
            !
            LL_OUT = .false.
            RP_O_TOP = min(RP_O(1,K,I),RP_O(2,K,I))
            RP_O_BOT = max(RP_O(1,K,I),RP_O(2,K,I))

            if (abs(RP_O_TOP).lt.1.e-6 .or. RPLT(1,I).gt.RP_O_TOP) then
               !***        Bring down the 0. or very small top boundary to the model roof
               LTOP=1
               rl_TOP=1.
               LL_OUT = .true.
            else if (RPLT(NIV,I).lt.RP_O_TOP) then
               !***        Bring up the very large bottom boundary to the model ground
               LTOP=NIV
               rl_TOP=0. ! Put zero since the LBOT index has already weight = 1
               LL_OUT = .true.
            else
               LTOP=1
               do while(RPLT(LTOP,I).lt.RP_O_TOP)
                  LTOP=LTOP+1
               end do
               !
               rl_TOP=(RPLT(LTOP,I)-RP_O_TOP)/ &
                  (RPLT(LTOP,I)-RPLT(LTOP-1,I))

            end if
            !
            if (RPLT(NIV,I).lt.RP_O_BOT) then
               !***        Bring up the very large bottom boundary to the model ground
               LBOT=NIV
               rl_BOT=1.
               RP_O_BOT = RPLT(NIV,I)
            else if (RPLT(1,I).gt.RP_O_BOT) then
               !***        Bring down the 0. or very small bottom boundary to the model roof
               LBOT=1
               rl_BOT=0. ! Put zero because the LTOP index has already weight 1
               LL_OUT = .true.
            else
               LBOT=NIV
               do while(RPLT(LBOT,I).gt.RP_O_BOT)
                  LBOT=LBOT-1
               end do
               LBOT=LBOT+1
               !
               rl_BOT=(RP_O_BOT-RPLT(LBOT-1,I))/ &
                  (RPLT(LBOT,I)-RPLT(LBOT-1,I))
            end if
            !
            ! 5. FILL AVK_KERNEL LINES
            !
            !        IF(.NOT.LL_OUT) THEN
            if (LTOP.eq.LBOT) then
               rl_TOP=(RP_O_BOT-RP_O_TOP)/ &
                  (RPLT(LTOP,I)-RPLT(LTOP-1,I))
               AVK(K,LTOP,I)=rl_TOP*AVK_COEFF(LTOP)
            else
               AVK(K,LTOP,I)=rl_TOP*AVK_COEFF(LTOP)
               do L = LTOP+1,LBOT-1
                  AVK(K,L,I)=AVK_COEFF(L)
               end do
               AVK(K,LBOT,I)=rl_BOT*AVK_COEFF(LBOT)
            end if
            !        END IF
         end do ! K=1,NIV_OBS
         if (LL_SWAP) then
            RP_SWAP(:,:)=AVK(:,:,I)
            do L=1,NIV
               AVK(:,L,I)=RP_SWAP(:,NIV-L+1)
            end do
         end if
      end do ! I = 1,NPRF
      !
      ! CLEAN UP
      !
      if (LL_SWAP) then
         deallocate(RP_SWAP)
         allocate(RP_SWAP(NIV,NPRF))
         RP_SWAP(:,:)=RP(:,:)
         do L=1,NIV
            RP(L,:)=RP_SWAP(NIV-L+1,:)
         end do
         !
         RP_SWAP(:,:)=RPLT(:,:)
         do L=1,NIV
            RPLT(L,:)=RP_SWAP(NIV-L+1,:)
         end do
         LL_SWAP = .true.
         deallocate(RP_SWAP)
      end if
      deallocate (ALFZ1MAS)
      deallocate (ALFZ2MAS)
      deallocate (ALFZ3MAS)
      deallocate (ALFZ4MAS)
      deallocate (ALFZ5MAS)
      !
   end subroutine TRANSF_TOTCOL_AVG_KERNEL
   !
   subroutine TRANSF_TOTCOL_WEIGHTS(NIV,NPRF,PL,PLT  &
      ,ALFZ1MAS,ALFZ2MAS  &
      ,ALFZ3MAS,ALFZ4MAS  &
      ,ALFZ5MAS)
      !
      implicit none
      !
      !
      !     Subroutine arguments.
      !     ~~~~~~~~~~~~~~~~~~~~
      !
      integer :: NIV,NPRF
      real (KIND=8) , dimension (NIV,NPRF) :: PL,PLT
      real (KIND=8) , dimension (NPRF,NIV) :: ALFZ1MAS,ALFZ2MAS
      real (KIND=8) , dimension (NPRF,NIV) :: ALFZ3MAS,ALFZ4MAS
      real (KIND=8) , dimension (NPRF,NIV) :: ALFZ5MAS
      !
      !     Variables locales.
      !     ~~~~~~~~~~~~~~~~~
      !
      integer :: I,J,L,LP1
      real (KIND=8) :: PL1,PL2,PLT0,PLT1,PLT2,COEF,DP2
      real (KIND=8) :: RPN,RPNM1,RPTN,RPTNM1
      real (KIND=8) :: RP1,RP2,RP3,RP4,DP1
      real (KIND=8) :: RP12,RP13,RP14,RP23,RP24,RP34
      real (KIND=8) , allocatable :: INVPROD(:,:),SOMPROD(:,:)
      real (KIND=8) , allocatable :: PROD(:,:),SOM(:,:)
      real (KIND=8) , allocatable :: ZCOEF(:)
      !
      allocate (INVPROD(NIV,4))
      allocate (SOMPROD(NIV,4))
      allocate (PROD(NIV,4))
      allocate (SOM(NIV,4))
      allocate (ZCOEF(4))
      !
      do I=1,NPRF
         !
         ! 1er cas: maille superieure L=1
         ! Valeur prise constante dans la premiere moitie de maille,
         ! puis interpolation lineaire entre L=1 et L~=1
         !
         PL1=PL(1,I)
         PL2=PL(2,I)
         PLT1=PLT(1,I)
         PLT0=0.
         COEF=((PLT1-PL1)**2)/(PL2-PL1)
         COEF=COEF/2.
         ALFZ1MAS(I,1)=0.
         ALFZ2MAS(I,1)=0.
         ALFZ3MAS(I,1)=(PLT1-PLT0)-COEF
         ALFZ4MAS(I,1)=COEF
         ALFZ5MAS(I,1)=0.
         !
         ! 2 eme cas: maille L=2.
         ! Interpolation lineaire dans la premiere moitie de maille
         ! puis cubique dans la seconde moitie (les coefficients
         ! beta correspondant a la seconde moitie sont calcules plus loin)
         !
         COEF=((PL2-PLT1)**2)/2.
         COEF=COEF/(PL2-PL1)
         ALFZ1MAS(I,2)=0.
         ALFZ2MAS(I,2)=COEF
         ALFZ3MAS(I,2)=(PL2-PLT1)-COEF
         ALFZ4MAS(I,2)=0.
         !
         ! Cas general: interpolation cubique dans les deux moities de mailles
         ! coefficients de la moitie superieure (courbe passant par les points
         ! L-2,L-1,L et L+1):ALF
         !
         ! La masse dans la maille sera calculee par:
         !  ALFZ1 X(L-2) + (ALFZ2)X(L-1) + (ALFZ3) X(L)
         ! + (ALFZ4) X(L+1) + ALFZ5 X(L+2)
         !
         ! Calcul de ALF(Gauche) de L=3 a NIV-1
         ! Calcul de ALF(Droite) de L=2 a NIV-2
         !
         do L=3,NIV-1
            RP1=PL(L-2,I)
            RP2=PL(L-1,I)
            RP3=PL(L,I)
            RP4=PL(L+1,I)
            DP1=PLT(L-1,I)
            RP12=RP1*RP2
            RP13=RP1*RP3
            RP14=RP1*RP4
            RP23=RP2*RP3
            RP24=RP2*RP4
            RP34=RP3*RP4
            INVPROD(L,1)=1./((RP1-RP2)*(RP1-RP3)*(RP1-RP4))
            INVPROD(L,2)=1./((RP2-RP1)*(RP2-RP3)*(RP2-RP4))
            INVPROD(L,3)=1./((RP3-RP1)*(RP3-RP2)*(RP3-RP4))
            INVPROD(L,4)=1./((RP4-RP1)*(RP4-RP2)*(RP4-RP3))
            SOM(L,1)=RP2+RP3+RP4
            SOM(L,2)=RP1+RP3+RP4
            SOM(L,3)=RP1+RP2+RP4
            SOM(L,4)=RP1+RP2+RP3
            PROD(L,1)=RP2*RP3*RP4
            PROD(L,2)=RP1*RP3*RP4
            PROD(L,3)=RP1*RP2*RP4
            PROD(L,4)=RP1*RP2*RP3
            SOMPROD(L,1)=RP23+RP24+RP34
            SOMPROD(L,2)=RP13+RP14+RP34
            SOMPROD(L,3)=RP12+RP14+RP24
            SOMPROD(L,4)=RP12+RP13+RP23
            do J=1,4
               ZCOEF(J)=             (RP3**4-DP1**4)/4.&
                  -SOM    (L,J)*(RP3**3-DP1**3)/3.&
                  +SOMPROD(L,J)*(RP3**2-DP1**2)/2.&
                  -PROD   (L,J)*(RP3   -DP1   )
               ZCOEF(J)=INVPROD(L,J)*ZCOEF(J)
            enddo
            ALFZ1MAS(I,L)=ZCOEF(1)
            ALFZ2MAS(I,L)=ZCOEF(2)
            ALFZ3MAS(I,L)=ZCOEF(3)
            ALFZ4MAS(I,L)=ZCOEF(4)
         enddo
         !
         do L=2,NIV-2
            LP1=L+1
            RP3=PL(L,I)
            RP4=PL(L+1,I)
            DP2=PLT(L,I)
            do J=1,4
               ZCOEF(J)=               (DP2**4-RP3**4)/4.&
                  -SOM    (LP1,J)*(DP2**3-RP3**3)/3.&
                  +SOMPROD(LP1,J)*(DP2**2-RP3**2)/2.&
                  -PROD   (LP1,J)*(DP2   -RP3   )
               ZCOEF(J)=INVPROD(LP1,J)*ZCOEF(J)
            enddo
            ALFZ2MAS(I,L)=ALFZ2MAS(I,L)+ZCOEF(1)
            ALFZ3MAS(I,L)=ALFZ3MAS(I,L)+ZCOEF(2)
            ALFZ4MAS(I,L)=ALFZ4MAS(I,L)+ZCOEF(3)
            ALFZ5MAS(I,L)=ZCOEF(4)
         enddo
         !
         ! Avant-derniere maille (pres du sol: L=NIV-1)
         ! Meme principe que pour L=2
         ! Les ALFZ ont ete calcules plus haut
         !
         RPN   =PL (NIV,I)
         RPNM1 =PL (NIV-1,I)
         RPTN  =PLT(NIV,I)
         RPTNM1=PLT(NIV-1,I)
         COEF=(RPTNM1**2-RPNM1**2)/2.-RPN*(RPTNM1-RPNM1)
         COEF=COEF/(RPNM1-RPN)
         ALFZ3MAS(I,NIV-1)=ALFZ3MAS(I,NIV-1)+               COEF
         ALFZ4MAS(I,NIV-1)=ALFZ4MAS(I,NIV-1)+(RPTNM1-RPNM1)-COEF
         ALFZ5MAS(I,NIV-1)=0.
         !
         ! Derniere maille (L=NIV)
         ! meme principe que pour L=1
         !
         COEF=(RPN**2-RPTNM1**2)/2.-RPN*(RPN-RPTNM1)
         COEF=COEF/(RPNM1-RPN)
         !
         ! essai d'interpolation lineaire jusqu'au sol
         ! sans correction en cas de passage sous 0
         !
         ALFZ1MAS(I,NIV)=0.
         ALFZ2MAS(I,NIV)=           -COEF
         ALFZ3MAS(I,NIV)=RPTN-RPTNM1+COEF
         ALFZ4MAS(I,NIV)=0.
         ALFZ5MAS(I,NIV)=0.
         !
      enddo
      !
      deallocate (INVPROD)
      deallocate (SOMPROD)
      deallocate (PROD)
      deallocate (SOM)
      deallocate (ZCOEF)
      !
   end subroutine TRANSF_TOTCOL_WEIGHTS

   subroutine observ_operator_tc2  (                 &
      ! INPUT MODEL DATA:
      nlat            & !dimensions:
      ,nlon            &
      ,nlev_m          &
      ,ntime           &
      ,latmod          & !data:
      ,lonmod          &
      ,pres_m          &
      ,pres_m_int      &
      ,date_MOCAGE     &
      ! INPUT OBS DATA:
      ,NPRF            & !dimensions:
      ,n_lvls_out      &
      ,lat             & !data:
      ,lon             &
      ,periodic_lon    &
      ,date            &
      ,pres_bounds     &
      ! EXECUTION OPTIONS:
      ,flag_temp_interp      &
      ! OUTPUT DATA:
      ,Hres_csr        & !Observational Op. in CSR format
      ,avg_kernel      & !Per prof. integration weights
      ,cell_fractions  & !Per subcolumn cell contribution
      ) !
      !
      use floats
      use interp_matrix_structure_mod
      use transformations_matrix
      use matrix_manipulations
      implicit none
      !-------------------------------------------
      ! DECLARATIONS:
      !-------------------------------------------
      ! INPUT DATA:
      !numb. of model longitudes and latitudes
      integer*4,intent(in)                      :: nlon,nlat
      !numb. of model veritcal levels
      integer*4,intent(in)                      :: nlev_m
      !numb. of temporary MOCAGE files: if 1 no temporal interpolation
      !                                 if 2 - temporal interpolation :
      integer*4,intent(in)                      :: ntime

      ! model lats and lons:
      real(std), dimension(nlat),intent(in):: latmod
      real(std), dimension(nlon),intent(in):: lonmod
      ! model fields:
      real(std), dimension(nlon*nlat*nlev_m*ntime),intent(in) :: pres_m
      real(std), dimension(nlon*nlat*nlev_m*ntime),intent(in) :: pres_m_int
      real(std), dimension(ntime),intent(in)    :: date_MOCAGE

      !numb. of obser. profiles in MOCAGE interval being treated:
      integer,intent(in)                        :: NPRF
      !this variable dimensions the AvKernel(n_lvls_out,n_lvls),=1 in case
      ! of integrated column data
      integer*4,intent(in)                      :: n_lvls_out

      real(std),dimension(NPRF),intent(in)      :: lat      !IASI lat
      real(std),dimension(NPRF),intent(in)      :: lon      !IASI lon
      logical,intent(in)                        :: periodic_lon
      real(std),dimension(NPRF),intent(in)      :: date     !IASI dates
      real(std),dimension(2,n_lvls_out,NPRF),intent(in) :: pres_bounds !col. bounds

      logical,intent(in)                      :: flag_temp_interp

      !-------------------------------------------
      ! OUTPUT DATA:
      type(CSR_FORMAT),intent(out)              :: Hres_csr

      !-------------------------------------------
      ! LOCAL VARIABLES:
      integer, pointer, dimension(:)        :: iw  ! work array for the routine amub
      ! THE FOUR MATRICES CORRESPONDING TO THE 4 OPERATORS, MATRICES IN CSR FORMAT:
      type(CSR_FORMAT)                          :: Hh_csr  !horizontal operator
      type(CSR_FORMAT)                          :: Ht_csr  !temporal operator
      type(CSR_FORMAT)                          :: Ha_csr  !smoothing operator
      ! TEMPORARY MATRICES FOR STORING THE RESULTS OF MULTIPLICATIONS OF MATRICES:
      !  TYPE(CSR_FORMAT)                          :: Hres2_csr  !result of Ht*(Hv*Hh)

      logical, dimension(NPRF)      :: discard_obs

      ! MODEL PRESSURE INTEGRATED TO THE LOCATION OF OBSERVATION
      real(std),pointer,dimension(:,:)      :: pres_inth
      real(std),pointer,dimension(:,:)      :: pres_inth_int

      real(std), dimension(nlev_m,NPRF), intent(out)  :: avg_kernel
      real(std), dimension(n_lvls_out,nlev_m,NPRF), intent(out)  :: cell_fractions
      !maximum numb. of obs vertical levels
      integer*4                                 :: n_lvls

      !# of model levels used in vertical interpolation,
      ! if the averaging kernel is used,
      ! this number is used to store the  number of non-zero elements in the final
      !  matrix
      integer*4                                 :: lvls_out
      integer*4                                 :: np   !indexing variable
      ! number of non-zero elements in the CSR format:
      integer*4                                 :: nmax
      integer*4                                 :: ierr !error index

      integer*4                                 :: ALLOC_ERR


      !******************************************************************
      ! END OF DECLARATIONS

      !-------------------------------------------
      ! 1) PREPARE THE MATRIX TO INTERPOLATE MODEL CO AND PRES FROM MODEL GRID TO
      !    OBS LOCATIONS
      call transf_horiz_interp_matrix(NPRF,nlon,nlat,nlev_m,ntime,      &
         latmod,lonmod,                      &
         lat,lon,                            &
         periodic_lon,                       &
         discard_obs,                        &
         Hh_csr)
      !
      !-------------------------------------------
      ! 2) NO NEED OF VERTICAL INTERPOLATION:
      !
      n_lvls = nlev_m
      !
      !-------------------------------------------
      ! 3) PREPARE THE MATRIX TO INTERPOLATE WITH RESPECT TO TIME

      if (flag_temp_interp) then   !if the temporal interpolation is
         ! necessary:

         call transf_temp_interp_matrix(n_lvls,                &
            NPRF,ntime,            &
            date_MOCAGE(1:2),      &
            date,                  &
            Ht_csr)
      end if
      !
      !-------------------------------------------
      !-------------------------------------------
      ! FINISHED CONSTRUCTING THE MATRICES, NOW MULTIPLY THEM TO OBTAIN THE FINAL
      ! OBSERVATIONAL OPERATOR H=Ha*Ht*Hh
      ! FOR THE INDIVIDUAL MATRICES IT WAS STRAIGHTFORWARD TO FIND THE NUMBER OF
      ! NON ZERO ELEMENTS. FOR THE PRODUCTS OF TWO MATRICES, THE NUMBER OF NON ZERO
      ! ELEMENTS CAN BE ALSO ESTIMATED ANALYTICALLY. THE FORMULAS ARE GIVEN
      ! BELOW.
      !-------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculate: H_t*H_h:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flag_temp_interp) then
         !the multiplication of Ht*Hres gives 16 non zero values in each row
         nmax = 16 * Ht_csr%nrow

         allocate(iw(Hh_csr%ncol))    !work array needed for routine amub
         iw=0
         call multiply_matrices_csr(iw,         &
            Ht_csr,     &    !matrix A
            Hh_csr,     &    !matrix B
            Hres_csr,  &    !out:A*B
            nmax)
         call deallocate_operator(Ht_csr)
         deallocate(iw)
      else
         call copy_matrix_csr(Hh_csr,Hres_csr)
      end if!if flag_temp_interp
      !
      call deallocate_operator(Hh_csr)
      !
      !-------------------------------------------
      ! need the interpolated pressure levels for the total column:
      allocate(pres_inth(Hres_csr%nrow,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"multiply_matrix_csr_vector")
      call multiply_matrix_csr_vector(Hres_csr,pres_m,          &
         Hres_csr%nrow,NPRF,pres_inth)
      !
      allocate(pres_inth_int(Hres_csr%nrow,NPRF),STAT = ALLOC_ERR)
      call check_allocate(ALLOC_ERR,"multiply_matrix_csr_vector")
      call multiply_matrix_csr_vector(Hres_csr,pres_m_int,          &
         Hres_csr%nrow,NPRF,pres_inth_int)
      !
      ! pres_inth now contains the model pressures on the obs horizontal grid
      ! pres_inth_int now contains the model interface press. on the obs hor. grid
      !
      !-------------------------------------------
      ! 4) Compute the avg kernel matrix for total columns (cubic),
      !    then convert it to the CSR format
      !
      call transf_totcol_avg_kernel2(nlev_m,n_lvls_out,NPRF,pres_inth, &
         pres_inth_int,pres_bounds,avg_kernel,&
         cell_fractions)
      !
      deallocate(pres_inth)
      deallocate(pres_inth_int)
      !
      !-------------------------------------------
      !
      ! finished calculating the H matrix,
      !-------------------------------------------
      !
      return
   end subroutine observ_operator_tc2
   !
   subroutine TRANSF_TOTCOL_AVG_KERNEL2(NIV,NIV_OBS,NPRF,RP,RPLT,RP_O,AVK,FRACS)
      !
      !CCCC
      !
      !   VERTICAL INTEGRATION (CUBIC)
      !   __  __   _____   _____   _____   _____   _____
      !  |  \/  | /     \ /     \ /  _  \ /  ___\ |  ___|
      !  |      | |  |  | |  <--< |  _  | |  \_ \ |  ___|
      !  |_|\/|_| \_____/ \_____/ \_/ \_/ \_____/ |_____|
      !
      !   MODELE DE CHIMIE ATMOSPHERIQUE A GRANDE ECHELLE
      !   VERSION_RELEASE_BUGFIX : 0_1_1
      !
      !CC DATE : 09/06/2005
      !
      !   Type      : Chemistry
      !   Purpose   : Compute a total column
      !   Called by :
      !
      !CC
      !
      !  Original : B. Josse & V.-H. Peuch, 04/2000
      !
      !  Modified :
      !
      !CCCC
      !
      implicit none
      !
      !
      !     Subroutine arguments.
      !     ~~~~~~~~~~~~~~~~~~~~
      !
      integer :: NIV, NIV_OBS, NPRF
      real (KIND=8) , dimension (NIV,NPRF) :: RP,RPLT
      real (KIND=8) , dimension (2,NIV_OBS,NPRF) :: RP_O
      real (KIND=8) , dimension (NIV,NPRF) :: AVK
      real (KIND=8) , dimension (NIV_OBS,NIV,NPRF) :: FRACS
      !
      ! LOOPS
      integer :: I,L,K
      ! INDEXES
      integer :: LTOP,LBOT
      ! COEFFICIENTS
      real(KIND=8) :: rl_TOP,rl_BOT
      ! CONSTANTS
      real(KIND=8) :: MAIR,GG,AVO,DOBSON
      ! WEIGHTS FOR NEIGHBORING LEVELS
      real (KIND=8), allocatable   :: ALFZ1MAS(:,:),ALFZ2MAS(:,:)
      real (KIND=8), allocatable   :: ALFZ3MAS(:,:),ALFZ4MAS(:,:)
      real (KIND=8), allocatable   :: ALFZ5MAS(:,:)
      ! WORK ARRAY
      real (KIND=8), allocatable :: RP_SWAP(:,:)
      real (KIND=8) :: RP_O_TOP, RP_O_BOT
      logical :: LL_SWAP, LL_OUT
      !
      allocate (ALFZ1MAS(NPRF,NIV))
      allocate (ALFZ2MAS(NPRF,NIV))
      allocate (ALFZ3MAS(NPRF,NIV))
      allocate (ALFZ4MAS(NPRF,NIV))
      allocate (ALFZ5MAS(NPRF,NIV))
      !
      ! 1. INITIALIZATIONS
      !
      do I = 1,NPRF
         do L=1,NIV
            AVK(L,I)=0.
            do K=1,NIV_OBS
               FRACS(K,L,I)=0.
            enddo
         enddo
      enddo
      !
      LL_SWAP = .false.
      if (RP(2,1) .lt. RP(1,1)) then
         allocate(RP_SWAP(NIV,NPRF))
         RP_SWAP(:,:)=RP(:,:)
         do L=1,NIV
            RP(L,:)=RP_SWAP(NIV-L+1,:)
         end do
         !
         RP_SWAP(:,:)=RPLT(:,:)
         do L=1,NIV
            RPLT(L,:)=RP_SWAP(NIV-L+1,:)
         end do
         LL_SWAP = .true.
         deallocate(RP_SWAP)
         allocate(RP_SWAP(NIV_OBS,NIV))
      end if

      !
      MAIR=0.029
      GG=9.80665
      AVO=6.02E23
      DOBSON=2.685E20
      !
      ! 2. COMPUTE WEIGHTS
      !
      call TRANSF_TOTCOL_WEIGHTS(NIV,NPRF,RP,RPLT &
         ,ALFZ1MAS,ALFZ2MAS &
         ,ALFZ3MAS,ALFZ4MAS &
         ,ALFZ5MAS)
      !
      ! 3. COMPUTE NUMBER OF MOLES IN EACH MODEL
      !    CELL (CUBIC INTERPOLATION)
      !
      do I = 1,NPRF
         AVK(1,I)=ALFZ3MAS(I,1)+ALFZ2MAS(I,2)+ALFZ1MAS(I,3)
         AVK(2,I)=ALFZ4MAS(I,1)+ALFZ3MAS(I,2)+ALFZ2MAS(I,3)+ALFZ1MAS(I,4)
         AVK(3,I)=ALFZ4MAS(I,2)+ALFZ3MAS(I,3)+ALFZ2MAS(I,4)+ALFZ1MAS(I,5)
         do L = 4,NIV-2
            AVK(L,I)=ALFZ5MAS(I,L-2)+ALFZ4MAS(I,L-1)+ALFZ3MAS(I,L)&
               +ALFZ2MAS(I,L+1)+ALFZ1MAS(I,L+2)
         end do
         AVK(NIV-1,I)=ALFZ5MAS(I,NIV-3)+ALFZ4MAS(I,NIV-2)&
            +ALFZ3MAS(I,NIV-1)+ALFZ2MAS(I,NIV)
         AVK(NIV,I)=ALFZ5MAS(I,NIV-2)+ALFZ4MAS(I,NIV-1)&
            +ALFZ3MAS(I,NIV)
         !
         AVK(:,I) = AVK(:,I)/GG/MAIR*AVO/DOBSON
         !
         ! 4. DETERMINE LEVELS/COLUMNS INTERSECTIONS
         !
         do K=1,NIV_OBS ! For each subcolumn
            !
            LL_OUT = .false.
            RP_O_TOP = min(RP_O(1,K,I),RP_O(2,K,I))
            RP_O_BOT = max(RP_O(1,K,I),RP_O(2,K,I))

            if (abs(RP_O_TOP).lt.1.e-6 .or. RPLT(1,I).gt.RP_O_TOP) then
               !***        Bring down the 0. or very small top boundary to the model roof
               LTOP=1
               rl_TOP=1.
               LL_OUT = .true.
            else if (RPLT(NIV,I).lt.RP_O_TOP) then
               !***        Bring up the very large bottom boundary to the model ground
               LTOP=NIV
               rl_TOP=0. ! Put zero since the LBOT index has already weight = 1
               LL_OUT = .true.
            else
               LTOP=1
               do while(RPLT(LTOP,I).lt.RP_O_TOP)
                  LTOP=LTOP+1
               end do
               !
               rl_TOP=(RPLT(LTOP,I)-RP_O_TOP)/ &
                  (RPLT(LTOP,I)-RPLT(LTOP-1,I))

            end if
            !
            if (RPLT(NIV,I).lt.RP_O_BOT) then
               !***        Bring up the very large bottom boundary to the model ground
               LBOT=NIV
               rl_BOT=1.
               RP_O_BOT = RPLT(NIV,I)
            else if (RPLT(1,I).gt.RP_O_BOT) then
               !***        Bring down the 0. or very small bottom boundary to the model roof
               LBOT=1
               rl_BOT=0. ! Put zero because the LTOP index has already weight 1
               LL_OUT = .true.
            else
               LBOT=NIV
               do while(RPLT(LBOT,I).gt.RP_O_BOT)
                  LBOT=LBOT-1
               end do
               LBOT=LBOT+1
               !
               rl_BOT=(RP_O_BOT-RPLT(LBOT-1,I))/ &
                  (RPLT(LBOT,I)-RPLT(LBOT-1,I))
            end if
            !
            ! 5. FILL AVK_KERNEL LINES
            !
            if(.not.LL_OUT) then
               if (LTOP.eq.LBOT) then
                  rl_TOP=(RP_O_BOT-RP_O_TOP)/ &
                     (RPLT(LTOP,I)-RPLT(LTOP-1,I))
                  FRACS(K,LTOP,I)=rl_TOP
               else
                  FRACS(K,LTOP,I)=rl_TOP
                  do L = LTOP+1,LBOT-1
                     FRACS(K,L,I)=1.
                  end do
                  FRACS(K,LBOT,I)=rl_BOT
               end if
            end if
         end do! K=1,NIV_OBS
         if (LL_SWAP) then
            RP_SWAP(:,:)=FRACS(:,:,I)
            do L=1,NIV
               FRACS(:,L,I)=RP_SWAP(:,NIV-L+1)
            end do
            RP_SWAP(:,:)=0
            RP_SWAP(1,:)=AVK(:,I)
            do L=1,NIV
               AVK(L,I)=RP_SWAP(1,NIV-L+1)
            end do
         end if
      end do! I = 1,NPRF
      !
      ! CLEAN UP
      !
      if (LL_SWAP) then
         deallocate(RP_SWAP)
         allocate(RP_SWAP(NIV,NPRF))
         RP_SWAP(:,:)=RP(:,:)
         do L=1,NIV
            RP(L,:)=RP_SWAP(NIV-L+1,:)
         end do
         !
         RP_SWAP(:,:)=RPLT(:,:)
         do L=1,NIV
            RPLT(L,:)=RP_SWAP(NIV-L+1,:)
         end do
         LL_SWAP = .true.
         deallocate(RP_SWAP)
      end if
      deallocate (ALFZ1MAS)
      deallocate (ALFZ2MAS)
      deallocate (ALFZ3MAS)
      deallocate (ALFZ4MAS)
      deallocate (ALFZ5MAS)
      !
   end subroutine TRANSF_TOTCOL_AVG_KERNEL2

   function get_number_col(Hv_csr,NPRF,discard_obs)
      ! This function returns  the number of model levels used in vertical interpolation
      ! to pressure levels.  This number is needed to determine the number of non-zero
      ! element in a matrix containing the result of the multiplication of the averaging kernel
      ! and the matrix containing the result of Ht*Hv*Hh.

      use interp_matrix_structure_mod

      implicit none

      ! input variables:
      type(CSR_FORMAT),intent(in)               :: Hv_csr  !vertical operator
      !maximum number of obs vertical levels:
      integer, intent(in)                       :: NPRF
      logical, dimension(NPRF),intent(in)       :: discard_obs


      integer*4                              :: jcol,col_in,col_out
      integer*4                              :: k,i,np
      integer*4, dimension(:), allocatable   :: iw
      integer*4 :: get_number_col


      allocate(iw(Hv_csr%ncol))

      col_out=0
      do np=1,Hv_csr%NPRF
         if (.not.discard_obs(np)) then
            iw(:) = 0
            col_in = 0
            do i=1,Hv_csr%nrow
               do k = Hv_csr%i(i,np),Hv_csr%i(i+1,np)-1
                  jcol = Hv_csr%j(k,np)
                  if (iw(jcol) .eq. 0) then
                     col_in = col_in + 1
                     iw(jcol) = 1
                  endif
               end do
            end do
            col_out=max(col_out,col_in)
         end if
      end do

      deallocate(iw)

      get_number_col=col_out
      return

   end function get_number_col

   subroutine dnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)

      implicit none

      real(KIND=8) :: dns(ndns,*),a(*)
      integer      :: ia(*),ja(*)
      integer      :: nrow,ncol,nzmax,ndns,ierr
      !-----------------------------------------------------------------------
      ! Dense           to    Compressed Row Sparse
      !-----------------------------------------------------------------------
      !
      ! converts a densely stored matrix into a row orientied
      ! compactly sparse matrix. ( reverse of csrdns )
      ! Note: this routine does not check whether an element
      ! is small. It considers that a(i,j) is zero if it is exactly
      ! equal to zero: see test below.
      !-----------------------------------------------------------------------
      ! on entry:
      !---------
      !
      ! nrow    = row-dimension of a
      ! ncol    = column dimension of a
      ! nzmax = maximum number of nonzero elements allowed. This
      !         should be set to be the lengths of the arrays a and ja.
      ! dns   = input nrow x ncol (dense) matrix.
      ! ndns    = first dimension of dns.
      !
      ! on return:
      !----------
      !
      ! a, ja, ia = value, column, pointer  arrays for output matrix
      !
      ! ierr    = integer error indicator:
      !         ierr .eq. 0 means normal retur
      !         ierr .eq. i means that the the code stopped while
      !         processing row number i, because there was no space left in
      !         a, and ja (as defined by parameter nzmax).
      !-----------------------------------------------------------------------
      integer :: next, i, j
      ierr = 0
      next = 1
      ia(1) = 1
      do i=1,nrow
         do j=1, ncol
            if (dns(i,j) .eq. 0.0d0) cycle
            if (next .gt. nzmax) then
               ierr = i
               return
            endif
            ja(next) = j
            a(next) = dns(i,j)
            next = next+1
         end do
         ia(i+1) = next
      end do
      return
      !---- end of dnscsr ----------------------------------------------------
      !-----------------------------------------------------------------------
   end subroutine dnscsr

end module space_time_operator_mod
