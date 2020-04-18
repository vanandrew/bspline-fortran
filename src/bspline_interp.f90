!*****************************************************************************************
!> author: Andrew Van
!  license: TBD
!
!### Description
!
!  Calls db3val on array of values
!
!### Notes
!   None

    subroutine db3interp(x,y,z,idx,idy,idz,&
                                    tx,ty,tz,ntx,nty,ntz,&
                                    nx,ny,nz,kx,ky,kz,bcoef,f,iflag,&
                                    inbvx,inbvy,inbvz,iloy,iloz,extrap,a,gs,gt,mask)
    use bspline_module
    use bspline_kinds_module, only: wp

    implicit none

    integer,intent(in)                      :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: ntx      !! the number of points to interpolate over in \(x\).
    integer,intent(in)                      :: nty      !! the number of points to interpolate over in \(y\).
    integer,intent(in)                      :: ntz      !! the number of points to interpolate over in \(z\).
    integer,intent(in)                      :: nx       !! the number of interpolation points in \(x\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: ny       !! the number of interpolation points in \(y\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: nz       !! the number of interpolation points in \(z\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: kx       !! order of polynomial pieces in \(z\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: ky       !! order of polynomial pieces in \(y\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: kz       !! order of polynomial pieces in \(z\).
                                                        !! (same as in last call to [[db3ink]])
    real(wp),dimension(nx),intent(in)       :: x        !! \(x\) coordinates of evaluation point.
    real(wp),dimension(ny),intent(in)       :: y        !! \(y\) coordinates of evaluation point.
    real(wp),dimension(nz),intent(in)       :: z        !! \(z\) coordinates of evaluation point.
    real(wp),dimension(nx+kx),intent(in)    :: tx       !! sequence of knots defining the piecewise polynomial
                                                        !! in the \(x\) direction. (same as in last call to [[db3ink]])
    real(wp),dimension(ny+ky),intent(in)    :: ty       !! sequence of knots defining the piecewise polynomial
                                                        !! in the \(y\) direction. (same as in last call to [[db3ink]])
    real(wp),dimension(nz+kz),intent(in)    :: tz       !! sequence of knots defining the piecewise polynomial
                                                        !! in the \(z\) direction. (same as in last call to [[db3ink]])
    real(wp),dimension(nx,ny,nz),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db3ink]].
    real(wp),dimension(4,4),intent(in)      :: a        !! affine matrix to change positions
    real(wp),dimension(ntx,nty,ntz),intent(out):: f     !! interpolated values
    integer,intent(out)                     :: iflag    !! status flag:
                                                        !!
                                                        !! * \( = 0 \)   : no errors
                                                        !! * \( \ne 0 \) : error
    integer,intent(inout)                   :: inbvx    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer,intent(inout)                   :: inbvy    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer,intent(inout)                   :: inbvz    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer,intent(inout)                   :: iloy     !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer,intent(inout)                   :: iloz     !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    logical,intent(in),optional             :: extrap   !! if extrapolation is allowed
                                                        !! (if not present, default is False)
    real(wp),dimension(4,4),intent(in)      :: gs       !! changes grid orientation (target)
    real(wp),dimension(4,4),intent(in)      :: gt       !! changes grid orientation (target)
    real(wp),dimension(ntx,nty,ntz),intent(in) :: mask     !! target mask

    integer :: i,j,k
    integer, save :: tflag,ivx,ivy,ivz,ily,ilz
    real(wp), save :: iv
    !$omp threadprivate(iv,tflag,ivx,ivy,ivz,ily,ilz)

    !! Initialize values of the changing parameters for threads
    ivx = inbvx
    ivy = inbvy
    ivz = inbvz
    ily = iloy
    ilz = iloz

    !$omp parallel copyin(ivx,ivy,ivz,ily,ilz)
    !$omp do
    do k=1,ntz
        do j=1,nty
            do i=1,ntx
                if (mask(i,j,k) > 0.5) then
                    call db3val((x(i)*gt(1,1)*a(1,1)+y(j)*gt(2,2)*a(1,2)+z(k)*gt(3,3)*a(1,3)+a(1,4))*gs(1,1),&
                        (x(i)*gt(1,1)*a(2,1)+y(j)*gt(2,2)*a(2,2)+z(k)*gt(3,3)*a(2,3)+a(2,4))*gs(2,2),&
                        (x(i)*gt(1,1)*a(3,1)+y(j)*gt(2,2)*a(3,2)+z(k)*gt(3,3)*a(3,3)+a(3,4))*gs(3,3),&
                        idx,idy,idz,tx,ty,tz,&
                        nx,ny,nz,kx,ky,kz,bcoef,iv,tflag,&
                        ivx,ivy,ivz,ily,ilz,extrap)
                else
                    iv = 0
                end if
                f(i,j,k) = iv
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    end subroutine db3interp