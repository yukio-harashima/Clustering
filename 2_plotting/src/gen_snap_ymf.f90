!****************************************************************
!
! Program: fgensnap_y_mf
! 
! _y is generate n v vector
! 
! Description: 
!   Generates time-dependent snapshot data for multi-segment 
!   fault models defined in Multiple_fault.dat.
!
! Key Features:
!   - Reads Multiple_fault.dat (3rd=dy, 4th=dx).
!   - Output snap_y.dat using local coordinates/parameters.
!   - Calculates spatial centroid using direct lat/lon weighting.
!   - Interpolates snap2_y_*.dat per segment with Virtual Origin.
!
!****************************************************************
program fgensnap_y_mf
   implicit none
   integer :: AveSnapFlag
   real    :: twrange, duration
   !---
   integer :: mn, nn, icmn, m, n, icm, jtn, jt, num_timewindow, tw, tcount
   integer :: m0, n0, ns, tmax, t, nx, ny, num_t_in_tw, cum_t_bf_tw, nexpo
   real    :: dm, xx, yy, rtime, strike0, dip0, rake0
   real    :: ylen, x1, x2, y1, y2, dx, dy, x, y, wk, slip_xy
   real    :: shift, strike, a_strike, dip, a_dip, rake, a_rake, moth, ch
   real    :: epilat, epilon, epidep, lat, lon, dep
   
   ! --- Multiple Fault Variables ---
   real, allocatable :: mf_dx(:,:), mf_dy(:,:)
   real, allocatable :: mf_lat(:,:), mf_lon(:,:), mf_depth(:,:)
   real, allocatable :: mf_stk(:,:), mf_dip(:,:), mf_rake(:,:)
   real, allocatable :: mf_atime(:,:)
   integer, allocatable :: mf_no(:,:)
   
   ! --- Interpolation & Loop Variables ---
   integer :: fn, max_fn, tm, tn, ierr, tno
   real    :: tdx, tdy, tlat, tlon, tdep, tstk, tdip, trake, tatime
   character(len=30) :: filename
   real    :: local_dx_min, local_dx_max, local_dy_min, local_dy_max
   real    :: r_dx_min, r_dy_min
   real    :: ref_lat, ref_lon, ref_dx, ref_dy, ref_stk, ref_dip
   integer :: vm0, vn0, ref_m_idx, ref_n_idx
   real    :: point1, point2
 
   ! --- Centroid Variables ---
   real    :: sum_moment, sum_lat, sum_lon, sum_dep
 
   ! --- Standard Arrays ---
   real, parameter :: rad = 3.1415927/180., dt = 0.1
   real, allocatable :: disp1(:, :, :), mij(:, :), mrf(:, :), mtx(:), mty(:)
   real, allocatable :: xmo(:), slip(:, :), disp2(:, :, :, :), tr(:, :)
   real, allocatable :: potency(:, :, :), sliprate(:, :, :)
   real, allocatable :: stf(:, :, :, :), istf(:, :, :, :)
   integer :: nodal_p
   
   !--- dummy matrix ---
   real, allocatable :: w(:), w_1(:)
   
   !-- For FPSPACK.FOR
   integer, parameter :: FPS = 1 
   real am,ama,am0,am1,e,am0b,slipa,&
        slipb,trendp,plungp,trendt,plungt,trendb,plungb,eta
   
   !--- Faultline vars ---
   integer :: min_m, max_m, min_n, max_n
   logical :: found_plane
 
   !---
   integer :: dum
   real    :: Mxy(1:3, 1:3), NDC, eigenval(1:3), EV(1:3, 1:3)
   real    :: n1_vector(3), n2_vector(3)
 
   !_________________________________________________________________
   ! 1. Read Inputs
   !_________________________________________________________________
   read(5,*,end=100)twrange, duration, AveSnapFlag
   
   ! Open static output files
   open(30, file = 'snap_y.dat', status = 'replace')
   open(32, file = 'tw_mec_xy.dat', status = 'replace')
   open(34, file = 'n_vector.dat', status = 'replace')
   
   !-- read fort.40 --
   open(40, file='fort.40', status='old')
   read(40, *); read(40, *)dm, dm, dm, epilat, epilon, epidep, dm, ns
   read(40, *); read(40, *)strike0, dip0, dm, rake0, dm
   read(40, *); read(40, *)xx, yy, mn, nn, m0, n0, rtime, jtn, icmn, ylen
   close(40)
   
   allocate(slip(mn, nn), tr(mn, nn))
   allocate(disp1(mn, nn, icmn), disp2(mn, nn, icmn, jtn))
   
   call readfrom40(mn, nn, icmn, jtn, slip, disp1, tr, disp2)
 
   ! --- Read Multiple_fault.dat ---
   allocate(mf_dx(1:mn, 1:nn), mf_dy(1:mn, 1:nn))
   allocate(mf_lat(1:mn, 1:nn), mf_lon(1:mn, 1:nn), mf_depth(1:mn, 1:nn))
   allocate(mf_stk(1:mn, 1:nn), mf_dip(1:mn, 1:nn), mf_rake(1:mn, 1:nn))
   allocate(mf_atime(1:mn, 1:nn), mf_no(1:mn, 1:nn))
   mf_no = -1
 
   open(50, file='Multiple_fault.dat', status='old')
   read(50, *, iostat=ierr) ! Skip Header
   max_fn = 0
   do
      ! Read: n, m, dy(3rd), dx(4th), lat, lon, depth, strike, dip, rake, ...
      read(50, *, iostat=ierr) tn, tm, tdy, tdx, tlat, tlon, tdep, tstk, tdip, trake, tatime, tno
      if (ierr /= 0) exit
      if (tm >= 1 .and. tm <= mn .and. tn >= 1 .and. tn <= nn) then
         mf_dx(tm, tn)    = tdx
         mf_dy(tm, tn)    = tdy
         mf_lat(tm, tn)   = tlat
         mf_lon(tm, tn)   = tlon
         mf_depth(tm, tn) = tdep
         mf_stk(tm, tn)   = tstk
         mf_dip(tm, tn)   = tdip
         mf_rake(tm, tn)  = trake
         mf_no(tm, tn)    = tno
         if (tno > max_fn) max_fn = tno
      end if
   end do
   close(50)
 
   if (AveSnapFlag /= 0 .and. AveSnapFlag /= 1) then
      write(6, *)" AveSnapFlag error. [0: cumulative] [1: moment] [else: stop]"
      stop
   end if
   num_timewindow = int( duration / twrange )
   write(6, *)" # of timewindow  :", num_timewindow
 
   !_________________________________________________________________
   ! 2. Reconstruct Source Time Function (STF)
   !_________________________________________________________________
   tmax = int(duration / dt) + 1 
   allocate(stf(0:mn+2, 0:nn+2, 0:tmax, icmn), w(0:tmax), w_1(0:tmax))
   stf = 0. ; w = 0. ; w_1 = 0.
   call stime1(w_1, tmax, dt, rtime, rtime*2.)
   do icm = 1, icmn
      do n = 1, nn
         do m = 1, mn
            do jt = 1, jtn
               shift = rtime * (jt - 1) + tr(m, n)
               call resample_shift(w_1, tmax, dt, w, tmax, dt, shift)
               stf(m, n, 0:tmax, icm) = stf(m, n, 0:tmax, icm) + disp2(m, n, icm, jt) * w(0:tmax)
            end do
         end do
      end do
   end do
   deallocate(w_1)
 
   !_________________________________________________________________
   ! 3. Process Snapshots (snap_y.dat)
   !_________________________________________________________________
   allocate(istf(0:mn+2, 0:nn+2, 1:num_timewindow, 1:icmn))
   allocate(potency(0:mn+2, 0:nn+2, 1:num_timewindow), sliprate(0:mn+2, 0:nn+2, 1:num_timewindow))
   allocate(w_1(1:3)) 
   istf = 0. ;  potency = 0. ; sliprate = 0.
 
   do n = 1, nn
      do m = 1, mn
         ! Skip if not in Multiple_fault
         if (mf_no(m,n) == -1) cycle
 
         do tw = 1, num_timewindow
            allocate(xmo(1:6), mij(1:3, 1:3), mrf(1:3, 1:3))
            xmo = 0.
            do icm = 1, icmn
               num_t_in_tw = int(twrange / dt)
               cum_t_bf_tw = int(((tw - 1) * twrange) / dt)
               tcount = 0
               ! Average or Moment selection
               if (AveSnapFlag == 0)then
                  do t = cum_t_bf_tw + 1, cum_t_bf_tw + num_t_in_tw
                     if (t >= tmax) exit
                     istf(m, n, tw, icm) = istf(m, n, tw, icm) + stf(m, n, t, icm)
                     tcount = tcount + 1
                  end do
                  istf(m, n, tw, icm) = istf(m, n, tw, icm) / tcount
               elseif (AveSnapFlag == 1)then
                  istf(m, n, tw, icm) = stf(m, n, cum_t_bf_tw + num_t_in_tw, icm)
               end if
               xmo(icm) = istf(m, n, tw, icm)
            end do
            
            ! --- Mechanism Calculation using LOCAL parameters ---
            call mtrx(xmo, mij)
            if( icmn .ne. 5 ) call fslip2mij(xmo,icmn,mf_stk(m,n),mf_dip(m,n),mf_rake(m,n), mij)
            call meca_rotation(mij, rad) 
            call mxy_mrf(mij, mrf)
            
            if (icmn == 5)then
               call d_cp(xmo, strike, dip, rake, sliprate(m, n, tw), dm, w_1)
               call conj(strike, dip, rake, a_strike, a_dip, a_rake)
            else
               call d_cp(xmo, strike, dip, rake, dm, dm, w_1)
               call conj(strike, dip, rake, a_strike, a_dip, a_rake)
               sliprate(m, n, tw) = sqrt(xmo(1)**2 + xmo(2)**2)
            end if
            
            if ( FPS .eq. 1 ) then 
              call mtrx(xmo, mij)
              call AR2PLP(mij, am0,am1,e,am0b,strike,dip,rake,slipa,&
                          a_strike,a_dip,a_rake,slipb,&
                          trendp,plungp,trendt,plungt,trendb,plungb,eta,ierr)
            end if
            
            potency(m, n, tw) = sliprate(m, n, tw) * xx * yy * 1.e6
 
            ! --- Output snap_y.dat ---
            if (sliprate(m, n, tw) .gt. 0.0000001)then
               call mrf_mxy(Mrf, Mxy)
               call EIG1(Mxy, 3, 3, eigenval, EV, dum)
               NDC = - eigenval(2) / max(abs(eigenval(1)), abs(eigenval(3)))
               NDC = 100 * NDC / 0.5 
               
               ! Use MF coordinates directly
               lat = mf_lat(m, n)
               lon = mf_lon(m, n)
               dep = mf_depth(m, n)
               
               call selection_nodal (strike, dip, a_strike, a_dip, mf_stk(m,n), mf_dip(m,n), rad, nodal_p)
               call swap_fault( strike, dip, rake, a_strike, a_dip, a_rake, nodal_p) 
               
               write(30, '(3i6, 2f9.2, 7f16.9, 4f8.2, 2f9.2, 3f9.3, 6f8.2, f9.2)') n, m, tw, &
                    mf_dx(m,n), mf_dy(m,n), & ! Local coords
                    sliprate(m, n, tw), mrf(1, 1), mrf(2, 2), mrf(3, 3), mrf(1, 2), &
                    mrf(1, 3), mrf(2, 3), strike, a_strike, dip, a_dip, rake, a_rake, &
                    lat, lon, dep, &
                    trendp,trendt,trendb,plungp,plungt,plungb,NDC
               
               ! n_vector calculation
               call get_n_vector (strike, dip, rad, n1_vector ) 
               call get_n_vector (a_strike, a_dip, rad, n2_vector ) 
               write(34,'(3i5,11f12.5)') n, m, tw, real(tw)/num_timewindow, & 
                    mf_dx(m,n), mf_dy(m,n), &
                    n1_vector(1), n1_vector(2), n1_vector(3), &
                    n2_vector(1), n2_vector(2), n2_vector(3), &
                    sliprate(m, n, tw)
            end if
            deallocate(xmo, mij, mrf)
         end do
      end do
   end do
   close(30)
   close(34)
 
   !_________________________________________________________________
   ! 4. Centroid (tw_mec_xy.dat) - Modified for Spatial Averaging
   !_________________________________________________________________
   do tw = 1, num_timewindow
      allocate(xmo(1:6), mij(1:3,1:3), mrf(1:3,1:3))
      xmo = 0
      
      sum_moment = 0.
      sum_lat = 0.
      sum_lon = 0.
      sum_dep = 0.
 
      do m = 1, mn
         do n = 1, nn
            if (mf_no(m, n) == -1) cycle
            
            ! Sum Moment Tensor
            do icm = 1, icmn
               xmo(icm) = xmo(icm) + istf(m, n, tw, icm)
            end do
            
            ! Sum for Spatial Centroid (Direct weighted average of lat/lon/dep)
            wk = sliprate(m, n, tw)
            sum_moment = sum_moment + wk
            sum_lat = sum_lat + wk * mf_lat(m, n)
            sum_lon = sum_lon + wk * mf_lon(m, n)
            sum_dep = sum_dep + wk * mf_depth(m, n)
         end do
      end do
      
      ! Mechanism Calculation (Same as original)
      call mtrx(xmo, mij)
      call meca_rotation(mij, rad) 
      call mxy_mrf(mij, mrf)
      call d_cp(xmo, strike, dip, rake, dm, dm, w_1)
      call conj(strike, dip, rake, a_strike, a_dip, a_rake)
      if ( FPS .eq. 1 ) then
        call mtrx(xmo, mij)
        call AR2PLP(mij, am0,am1,e,am0b,strike,dip,rake,slipa,&
                    a_strike,a_dip,a_rake,slipb,&
                    trendp,plungp,trendt,plungt,trendb,plungb,eta,ierr)
      end if
      
      ! Calculate Centroid
      if (sum_moment > 0.) then
         lat = sum_lat / sum_moment
         lon = sum_lon / sum_moment
         dep = sum_dep / sum_moment
         ! Use dummy x,y for output compatibility (or set to 0)
         x = 0. ; y = 0. 
      else
         lat = epilat; lon = epilon; dep = epidep
         x = 0. ; y = 0.
      end if
 
      nexpo = int(log10(maxval(abs(mrf))))
      mrf = mrf / (10.**nexpo)
      
      ! Output with calculated lat/lon/dep
      write(32,'(i5,8f12.6, 6f9.2, 3f9.3)')tw, x, y, &
           mrf(1,1), mrf(2,2), mrf(3,3), mrf(1,2), mrf(1,3),mrf(2,3), &
           strike, a_strike, dip, a_dip, rake, a_rake, lat, lon, dep
      deallocate(xmo, mij, mrf)
   end do
   close(32)
 
   !_________________________________________________________________
   ! 5. Snap2_y (Interpolated Contours) - Per Segment
   !_________________________________________________________________
   dx = 0.5 ; dy = 0.5
 
   do fn = 1, max_fn
      write(filename, '("snap2_y_", i0, ".dat")') fn
      open(31, file=trim(filename), status='replace')
      
      ! --- Determine Boundaries ---
      local_dx_min = 1.e10; local_dx_max = -1.e10
      local_dy_min = 1.e10; local_dy_max = -1.e10
      ref_lat = -999.
 
      do n=1,nn
         do m=1,mn
            if (mf_no(m,n) == fn) then
               if (mf_dx(m,n) < local_dx_min) local_dx_min = mf_dx(m,n)
               if (mf_dx(m,n) > local_dx_max) local_dx_max = mf_dx(m,n)
               if (mf_dy(m,n) < local_dy_min) local_dy_min = mf_dy(m,n)
               if (mf_dy(m,n) > local_dy_max) local_dy_max = mf_dy(m,n)
               
               if (ref_lat == -999.) then
                  ref_lat = mf_lat(m,n)
                  ref_lon = mf_lon(m,n)
                  ref_dx  = mf_dx(m,n)
                  ref_dy  = mf_dy(m,n)
                  ref_stk = mf_stk(m,n)
                  ref_dip = mf_dip(m,n)
                  ref_m_idx = m
                  ref_n_idx = n
               endif
            endif
         enddo
      enddo
      
      if (ref_lat == -999.) then
          close(31)
          cycle
      endif
 
      ! --- Virtual Origin Setup ---
      vm0 = ref_m_idx - nint(ref_dx / xx)
      vn0 = ref_n_idx - nint(ref_dy / yy)
 
      ! --- Rounding Start Point (NINT & 0.1km step) ---
      ! fgenpdtdis_mfと同様に、データ最小値から1グリッド分戻した位置を基準とする
      r_dx_min = real(nint((local_dx_min - xx) * 10.)) / 10.
      r_dy_min = real(nint((local_dy_min - yy) * 10.)) / 10.
 
      ! --- Grid Size Calculation (Modified Safety Margin) ---
      ! 安全マージンを fgenpdtdis_mf と同じ + 1.0 に変更
      ! 範囲: (最大値 + xx) - 開始点(r_dx_min)
      nx = int((local_dx_max + xx - r_dx_min) * 2. + 1.0)
      ny = int((local_dy_max + yy - r_dy_min) * 2. + 1.0)
 
      ! --- Time Loop ---
      do tw = 1, num_timewindow
         do n = 1, ny
            do m = 1, nx
               point1 = r_dx_min + dx * (m - 1)
               point2 = r_dy_min + dx * (n - 1)
               
               wk = slip_xy(sliprate(0:mn+2, 0:nn+2, tw), mn, nn, vm0, vn0, xx, yy, point1, point2)
               
               if (wk >= 0.) then
                  call xy2geo(ref_lat, ref_lon, point1 - ref_dx, point2 - ref_dy, &
                              ref_stk, ref_dip, lat, lon)
                  dep = epidep - (point2 * sin(ref_dip * rad))
                  
                  write(31,'(3i5,2f11.3,f11.6, 3f9.3)') n, m, tw, point1, point2, wk, lat, lon, dep
               end if
            end do
         end do
         write(31, *) 
         write(31, *) 
      end do
      close(31)
   end do
 
   !_________________________________________________________________
   ! 6. Faultline (Per Segment)
   !_________________________________________________________________
   open(33, file='faultline.dat', status='replace')
   do fn = 1, max_fn
      min_m = 100000 ; max_m = -100000
      min_n = 100000 ; max_n = -100000
      found_plane = .false.
      do n = 1, nn
         do m = 1, mn
            if (mf_no(m, n) == fn) then
               if (m < min_m) min_m = m
               if (m > max_m) max_m = m
               if (n < min_n) min_n = n
               if (n > max_n) max_n = n
               found_plane = .true.
            end if
         end do
      end do
      if (found_plane) then
         write(33, '(3f12.4)') mf_lat(min_m, min_n), mf_lon(min_m, min_n), mf_depth(min_m, min_n)
         write(33, '(3f12.4)') mf_lat(max_m, min_n), mf_lon(max_m, min_n), mf_depth(max_m, min_n)
         write(33, '(3f12.4)') mf_lat(max_m, max_n), mf_lon(max_m, max_n), mf_depth(max_m, max_n)
         write(33, '(3f12.4)') mf_lat(min_m, max_n), mf_lon(min_m, max_n), mf_depth(min_m, max_n)
         write(33, '(3f12.4)') mf_lat(min_m, min_n), mf_lon(min_m, min_n), mf_depth(min_m, min_n)
         write(33, *) ' > ' 
         write(33, *) 
      end if
   end do
   close(33)
 
   stop
 100 stop "failed to read arguments. stop."
 end program fgensnap_y_mf
 
 ! Subroutine: get_n_vector
 subroutine get_n_vector (str, dip, rad ,n1 )
   implicit none 
   real,intent(in):: str, dip, rad
   real,intent(out):: n1(3)
   n1(1) = -sin(dip*rad)*sin(str*rad)  
   n1(2) =  sin(dip*rad)*cos(str*rad)  
   n1(3) = -cos(dip*rad)
   return
 end subroutine get_n_vector