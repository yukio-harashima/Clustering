!****************************************************************
!
! This program is based on Snippets / src / PotMecFaultPara.f90
!
!****************************************************************
program fgenpdtdis
  implicit none
  integer :: m, n, mn, nn, icmn, m0, n0, jtn, i, ns
  real :: ylen, rtime, rake
  real :: strike, dip, xx, yy, epilat, epilon, epidepth
  real, parameter   :: dx = 0.5d0, rad = 3.1415927/180.
  real, allocatable :: lat(:, :), lon(:, :), depth(:, :)
  real, allocatable :: slip(:, :)           !cumulate slip on each knot
  real, allocatable :: tr(:, :)             !start time of each knot
  real, allocatable :: i_slip(:, :, :)      !cumulate slip at each knot and slip comp.
  real, allocatable :: sliprate(:, :, :, :) !sliprate at each knot,timestep and slip comp.
  real, allocatable :: mij(:, :), tmp_mrf(:, :), xmo(:), potency(:, :), mrf(:, :, :, :)
  real, allocatable :: str1(:, :), dip1(:, :), rake1(:, :) ! conjugated faultparameter 1
  real, allocatable :: str2(:, :), dip2(:, :), rake2(:, :) ! conjugated faultparameter 2
  integer :: nodal_p
  !---
  integer :: nx, ny, nexpo
  real    :: dm, point1, point2, tmp_potency, lat_tmp, lon_tmp, slip_xy
  real    :: lonmax, lonmin, latmax, latmin, rad_dip, tmp_depth
  real    :: lat_w, lat_s, lat_n, lat_e, w_depth, s_depth, n_depth, e_depth
  real    :: lon_w, lon_s, lon_n, lon_e
  real    :: ax(3)
  !-- For FPSPACK.FOR
  integer, parameter :: FPS = 1 ! 1 uses FPSPACK. 0 does not use it.
  real am,ama,am0,am1,e,am0b,slipa,&
       slipb,eta
  real, allocatable :: trendp(:, :),plungp(:, :),trendt(:, :)
  real, allocatable :: plungt(:, :),trendb(:, :),plungb(:, :)
  integer ierr
  !-- For correct fault depth
  real    :: min_depth, max_depth
  !---
  integer :: dum
  real    :: Mxy(1:3, 1:3), NDC, eigenval(1:3), EV(1:3, 1:3)
  !---
  open(42,file='pdtdis.dat',status='replace')
  open(43,file='pddis.dat',status='replace')
  open(44,file='faultline.dat',status='replace')
  !---
  read(40, *); read(40, *)dm, dm, dm, epilat, epilon, epidepth, dm, ns
  read(40, *); read(40, *)strike, dip, dm, rake, dm
  read(40, *); read(40, *)xx, yy, mn, nn, m0, n0, rtime, jtn, icmn, ylen
  close(40)
  allocate(lat(1:mn, 1:nn), lon(1:mn, 1:nn), depth(1:mn, 1:nn))
  allocate(slip(1:mn, 1:nn), tr(1:mn, 1:nn), i_slip(1:mn, 1:nn, 1:icmn))
  allocate(sliprate(1:mn, 1:nn, 1:icmn, 1:jtn), mrf(1:mn, 1:nn, 1:3, 1:3))
  allocate(str1(1:mn, 1:nn), str2(1:mn, 1:nn), dip1(1:mn, 1:nn), dip2(1:mn, 1:nn))
  allocate(rake1(1:mn, 1:nn), rake2(1:mn, 1:nn), potency(0:mn+2, 0:nn+2))
  allocate(trendp(1:mn, 1:nn),plungp(1:mn, 1:nn),trendt(1:mn, 1:nn))
  allocate(plungt(1:mn, 1:nn),trendb(1:mn, 1:nn),plungb(1:mn, 1:nn))
  potency = 0. ; mrf = 0.
  call readfrom40(mn, nn, icmn, jtn, slip, i_slip, tr, sliprate)
  call getknot(strike, dip, mn, nn, m0, n0, xx, yy,  &
       epilat, epilon, epidepth, lat, lon, depth)
  !---
  do n=1,nn
    do m=1,mn
      allocate(xmo(1:6), mij(1:3, 1:3), tmp_mrf(1:3, 1:3))
      xmo = 0.
      do i = 1, icmn
        xmo(i) = i_slip(m, n, i)
      enddo
      call mtrx(xmo, mij)
      if( icmn .ne. 5 )call fslip2mij(xmo,icmn,strike,dip,rake, mij)
      call meca_rotation(mij, rad) ! add for meca rotation
      call mxy_mrf(mij, tmp_mrf)
      call d_cp(xmo, str1(m, n), dip1(m, n), rake1(m, n), potency(m, n), dm, ax)
      call conj(str1(m, n), dip1(m, n), rake1(m, n), str2(m, n), dip2(m, n), rake2(m, n))
      if ( FPS .eq. 1 ) then ! add for FPSPACK.FOR
        call mtrx(xmo, mij)
        call AR2PLP(mij, am0,am1,e,am0b,str1(m, n),dip1(m, n),rake1(m, n),slipa,&
                    str2(m, n),dip2(m, n),rake2(m, n),slipb,&
                    trendp(m, n),plungp(m, n),trendt(m, n),plungt(m, n),&
                    trendb(m, n),plungb(m, n),eta,ierr)
      end if
      !---
      if (icmn.eq.2) then
        potency(m, n) = sqrt(xmo(1)**2 + xmo(2)**2)
      endif
      nexpo = int(log10(maxval(abs(tmp_mrf))))
      !      tmp_mrf = tmp_mrf / (10.**nexpo)
      mrf(m, n, 1:3, 1:3) = tmp_mrf(1:3, 1:3)
      deallocate(xmo, mij, tmp_mrf)
    end do
  end do
  do n=1,nn
    do m=1,mn
      if (potency(m,n) > 0.) then
        !---
        call mrf_mxy(Mrf(m, n, 1:3, 1:3), Mxy)
        call EIG1(Mxy, 3, 3, eigenval, EV, dum)
        ! based on Dziewonski+1981JGR
        NDC = - eigenval(2) / max(abs(eigenval(1)), abs(eigenval(3)))
        NDC = 100 * NDC / 0.5 ! convert to parcentage
        !---
        !--Only works well if the mechanism is not rotating.
        call selection_nodal (str1(m,n), dip1(m,n), str2(m,n), dip2(m,n), strike, dip, rad ,nodal_p )
        call swap_fault (str1(m,n), dip1(m,n), rake1(m,n), str2(m,n), dip2(m,n), rake2(m,n) ,nodal_p )
        !---
        write(42, '(2i5, 2f9.2, 2f12.5, f10.3, f11.4, 6f16.9, i5, 4f8.2, 2f9.2, 6f8.2, f9.2)') n, m, &
            yy * real(n-n0), xx * real(m-m0), &
            lat(m, n), lon(m, n), depth(m, n), potency(m, n), &
            mrf(m, n, 1, 1), mrf(m, n, 2, 2), mrf(m, n, 3, 3), &
            mrf(m, n, 1, 2), mrf(m, n, 1, 3), mrf(m, n, 2, 3), 20, &
            str1(m, n),str2(m, n),dip1(m, n),dip2(m, n),rake1(m, n),rake2(m, n), &
            trendp(m, n),trendt(m, n),trendb(m, n),&
            plungp(m, n),plungt(m, n),plungb(m, n),NDC
      end if
    end do
  end do
  !---
  point1 = -xx * real(m0)
  point2 =  xx * (real(mn - m0) + 1.)
  nx = int((point2 - point1) * 2. + 1)
  ! nx = int((point2 - point1) * 2. + 2) to prevent from NaN value
  point1 = -yy * real(n0)
  if (ns == 1 .or. ns == 2)then
     point2 = yy * real(nn - n0) + ylen
  else
     point2 = yy * (real(nn - n0) + 1.)
  endif
  ny = int((point2 - point1) * 2. + 1)
  ! ny = int((point2 - point1) * 2. + 2) to prevent from NaN value
  latmin =  1000
  latmax = -1000
  lonmin =  1000
  lonmax = -1000
  !
  min_depth = 1000
  max_depth = -1000
  !
  rad_dip = dip * rad
  do n=1,ny
    do m=1,nx
      point1 = int( (-xx * real(m0)) * 5)/5 + dx * (m - 1)
      point2 = int( (-yy * real(n0)) * 5)/5 + dx * (n - 1)
      tmp_depth = epidepth - (point2 * sin(rad_dip))
      tmp_potency = slip_xy(potency(0, 0), mn, nn, m0, n0, xx, yy, point1, point2)
      call xy2geo(epilat, epilon, point1, point2, strike, dip, lat_tmp, lon_tmp)
      write(43, '(4f9.3, f13.6, f9.3)') lat_tmp, lon_tmp, point1, point2, tmp_potency, tmp_depth
      if (lon_tmp < 0.) then
        lon_tmp=360.+lon_tmp
      end if
      if (latmin > lat_tmp) then
         latmin = lat_tmp
         lat_s  = lat_tmp
         lon_s  = lon_tmp
         s_depth = tmp_depth
      elseif (latmax < lat_tmp) then
         latmax = lat_tmp
         lat_n  = lat_tmp
         lon_n  = lon_tmp
         n_depth = tmp_depth
      elseif (lonmin > lon_tmp) then
         lonmin = lon_tmp
         lat_w  = lat_tmp
         lon_w  = lon_tmp
         w_depth = tmp_depth
      elseif (lonmax < lon_tmp) then
         lonmax = lon_tmp
         lat_e  = lat_tmp
         lon_e  = lon_tmp
         e_depth = tmp_depth
      endif
      !
      if (tmp_depth < min_depth) min_depth=tmp_depth
      if (tmp_depth > max_depth) max_depth=tmp_depth
      !
    end do
  end do

  if ( mod(strike, 90.) == 0 )then
     write(44, '(3f10.4)')lat_n, lon_w, min_depth !w_depth
     write(44, '(3f10.4)')lat_n, lon_e, max_depth !e_depth
     write(44, '(3f10.4)')lat_s, lon_e, max_depth !e_depth
     write(44, '(3f10.4)')lat_s, lon_w, min_depth !w_depth
     write(44, '(3f10.4)')lat_n, lon_w, min_depth !w_depth
  else
     write(44, '(3f10.4)')lat_w, lon_w, min_depth !w_depth
     write(44, '(3f10.4)')lat_n, lon_n, max_depth !n_depth
     write(44, '(3f10.4)')lat_e, lon_e, max_depth !e_depth
     write(44, '(3f10.4)')lat_s, lon_s, min_depth !s_depth
     write(44, '(3f10.4)')lat_w, lon_w, min_depth !w_depth
  endif
  !---
  close(42)
  close(43)
  close(44)
  !---
  call faultline(strike,dip,epilat,epilon,epidepth,xx,yy,mn,nn,m0,n0,1.0)
end program
