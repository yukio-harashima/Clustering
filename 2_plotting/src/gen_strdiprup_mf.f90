!****************************************************************
!
! Program: fgenstrdiprup_mf (Hybrid Unfolding Version)
! Description: 
!   Generates Spatio-temporal rupture propagation data (VR).
!
!   Modifications:
!     - [HYBRID UNFOLDING]
!       Calculates distance based on physical coordinate differences (mf_dx/mf_dy).
!       1. If diff > threshold: Uses physical distance (Handles gaps/jumps).
!       2. If diff ~ 0: Uses grid spacing xx/yy (Handles identical coords).
!
!****************************************************************
program fgenstrdiprup_mf
   implicit none
   real    :: duration
   !---
   integer :: mn, nn, icmn, m, n, i, icm, jtn, jt, m0, n0, ns, tmax, t, nx, ny, tmax2
   real  :: dm, xx, yy, rtime, x1, x2, y1, y2
   real  :: ylen, shift, loc, val, strike0, dip0
   real, parameter   :: dt = 0.1, d = 0.1
   real, parameter   :: rad = 3.1415927/180.
   
   ! --- Arrays ---
   real, allocatable :: disp1(:, :, :), mij(:, :), mrf(:, :)
   real, allocatable :: xmo(:), slip(:, :), disp2(:, :, :, :), tr(:, :)
   real, allocatable :: sliprate(:, :, :), i_sliprate(:, :, :, :), w(:), w_1(:)
   real, allocatable :: str_sliprate(:, :), dip_sliprate(:, :)
   
   ! --- Multiple_fault Variables ---
   real, allocatable :: mf_dx(:,:), mf_dy(:,:)
   real, allocatable :: mf_lat(:,:), mf_lon(:,:), mf_depth(:,:)
   real, allocatable :: mf_stk(:,:), mf_dip(:,:), mf_rake(:,:)
   real, allocatable :: mf_atime(:,:)
   integer, allocatable :: mf_no(:,:)
   
   ! --- Unfolding Coordinates ---
   real, allocatable :: unfolded_dx(:,:), unfolded_dy(:,:)
   
   ! --- Loop & Temp Variables ---
   integer :: fn, max_fn, tm, tn, ierr, tno
   real    :: tdx, tdy, tlat, tlon, tdep, tstk, tdip, trake, tatime
   character(len=30) :: filename1, filename2
   real    :: local_dx_min, local_dx_max, local_dy_min, local_dy_max
   integer :: min_m, max_m, min_n, max_n
   real    :: epilat, epilon, epidepth
   
   ! --- Interpolation Variables ---
   real :: x_start, x_end, y_start, y_end, rate
   integer :: n_step
   logical :: valid_start, valid_end
   real :: loc_phys, loc_comb
 
   ! --- Variables for Combined Output ---
   integer, parameter :: n_bin = 400000 
   integer :: offset_bin = 200000
   real, allocatable :: vr_str_all(:,:), vr_dip_all(:,:)
   integer :: idx_x, idx_y
   integer :: min_valid_x, max_valid_x, min_valid_y, max_valid_y
   
   ! --- Offset Variables ---
   real, allocatable :: offset_x(:), offset_y(:)
   real, allocatable :: seg_min_x(:), seg_max_x(:)
   real, allocatable :: seg_min_y(:), seg_max_y(:)
   real :: current_max_x, current_max_y, shift_x, shift_y, dist_step
   real :: base_dx, base_dy
   
   ! --- Reference Plane Parameters ---
   real :: ref_strike, ref_dip
   logical :: file_exist
   character(len=256) :: bash_file, line_buf
   integer :: pos_eq
 
   !- - - - - - - - - - - - - - - - - - - - - - - - - - -
   read(5,*,end=100) duration
   !- - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   ! Read fort.40
   open(40, file='fort.40', status='old')
   read(40, *) ; read(40, *)dm, dm, dm, epilat, epilon, epidepth, dm, ns
   read(40, *) ; read(40, *)strike0, dip0, dm, dm, dm
   read(40, *) ; read(40, *)xx, yy, mn, nn, m0, n0, rtime, jtn, icmn, ylen
   close(40)
   
   ! Allocate Standard Arrays
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
   ! Removed header skip line
   ! read(50, *, iostat=ierr) 
   max_fn = 0
   do
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
 
   !--GET SOURCE TIME FUNCTION--
   tmax = int(duration / dt) + 1
   allocate(i_sliprate(0:mn+2, 0:nn+2, 0:tmax, icmn), w(0:tmax), w_1(0:tmax))
   i_sliprate = 0. ; w = 0. ; w_1 = 0.
   
   call stime1(w_1, tmax, dt, rtime, rtime*2.)
   do icm = 1, icmn
      do n = 1, nn
         do m = 1, mn
            do jt = 1, jtn
               shift = rtime * (jt - 1) + tr(m, n)
               call resample_shift(w_1, tmax, dt, w, tmax, dt, shift)
               i_sliprate(m,n,0:tmax,icm) = i_sliprate(m,n,0:tmax,icm) + disp2(m,n,icm,jt) * w(0:tmax)
            end do
         end do
      end do
   end do
   deallocate(w_1)
   
   ! -- Calculate Scalar Slip Rate --
   allocate(xmo(1:6), mij(1:3, 1:3), mrf(1:3, 1:3), sliprate(0:mn+2, 0:nn+2, 0:tmax))
   allocate(w_1(1:3))
   sliprate = 0.
   do t = 0, tmax
      do m = 1, mn
         do n = 1, nn
            xmo = 0.
            if (icmn == 5)then
               do icm = 1, icmn
                  xmo(icm) = i_sliprate(m, n, t, icm)
               end do
               call D_CP(xmo, dm, dm, dm, sliprate(m, n, t), dm, w_1)
            elseif (icmn == 2)then
               sliprate(m,n,t) = sqrt(i_sliprate(m,n,t,1)**2 + i_sliprate(m,n,t,2)**2)
            elseif (icmn == 1)then
               sliprate(m,n,t) = i_sliprate(m,n,t,1)
            end if
         end do
      end do
   end do
 
   !=================================================================
   ! 0. INTERNAL UNFOLDING LOGIC (Hybrid: Coords & Grid)
   !=================================================================
   allocate(unfolded_dx(mn, nn), unfolded_dy(mn, nn))
   unfolded_dx = 0. ; unfolded_dy = 0.
   
   allocate(seg_min_x(max_fn), seg_max_x(max_fn))
   allocate(seg_min_y(max_fn), seg_max_y(max_fn))
   
   do fn = 1, max_fn
      ! Find range for this fault
      min_m = 10000; max_m = -10000
      min_n = 10000; max_n = -10000
      local_dx_min = 1.e10 ; local_dy_min = 1.e10
      
      do n = 1, nn
         do m = 1, mn
            if (mf_no(m,n) == fn) then
               if (m < min_m) min_m = m
               if (m > max_m) max_m = m
               if (n < min_n) min_n = n
               if (n > max_n) max_n = n
               ! Capture absolute physical min
               if (mf_dx(m,n) < local_dx_min) local_dx_min = mf_dx(m,n)
               if (mf_dy(m,n) < local_dy_min) local_dy_min = mf_dy(m,n)
            end if
         end do
      end do
      if (min_m > max_m) cycle
      
      ! Save local minimums
      seg_min_x(fn) = local_dx_min
      seg_min_y(fn) = local_dy_min
 
      ! --- DEBUG: Print Physical Range ---
      print *, "--- Fault:", fn, " Range ---"
      print *, "  Grid m:", min_m, " to ", max_m
      print *, "  Phys Start X:", local_dx_min
 
      ! --- Unfold Strike (Along m) ---
      unfolded_dx(min_m, min_n:max_n) = 0.0
      
      do m = min_m, max_m - 1
         do n = min_n, max_n
            if (mf_no(m,n) == fn .and. mf_no(m+1,n) == fn) then
               ! 1. Calculate physical difference
               dist_step = abs(mf_dx(m+1, n) - mf_dx(m, n))
               
               ! 2. HYBRID CHECK:
               !    If dist is tiny (duplicate coords), use Grid Spacing (xx).
               !    If dist is large (jump/gap), use Physical Distance (dist_step).
               if (dist_step < 0.001) then
                   dist_step = xx
                   if (fn == 2 .and. n == min_n) &
                      print *, "  m=", m, " Duplicate dx. Using xx."
               endif
               
               unfolded_dx(m+1, n) = unfolded_dx(m, n) + dist_step
               
            else
               ! Gap in fault index (should not happen in valid segment)
               if (mf_no(m+1,n) == fn) then
                  unfolded_dx(m+1, n) = unfolded_dx(m, n) + xx
               endif
            endif
         end do
      end do
      seg_max_x(fn) = maxval(unfolded_dx(min_m:max_m, min_n:max_n))
      print *, "  Unfolded Max Length:", seg_max_x(fn)
      
      ! --- Unfold Dip (Along n) ---
      unfolded_dy(min_m:max_m, min_n) = 0.0
      do n = min_n, max_n - 1
         do m = min_m, max_m
            if (mf_no(m,n) == fn .and. mf_no(m,n+1) == fn) then
               ! 1. Calculate physical difference
               dist_step = abs(mf_dy(m, n+1) - mf_dy(m, n))
               
               ! 2. HYBRID CHECK for Dip
               if (dist_step < 0.001) dist_step = yy
               
               unfolded_dy(m, n+1) = unfolded_dy(m, n) + dist_step
            else
               if (mf_no(m,n+1) == fn) unfolded_dy(m, n+1) = unfolded_dy(m, n) + yy
            endif
         end do
      end do
      seg_max_y(fn) = maxval(unfolded_dy(min_m:max_m, min_n:max_n))
   end do
 
   !=================================================================
   ! PRE-CALCULATION: Determine Offsets for Global Layout
   !=================================================================
   allocate(offset_x(max_fn), offset_y(max_fn))
   offset_x = 0. ; offset_y = 0.
   
   current_max_x = seg_max_x(1) + seg_min_x(1) 
   current_max_y = seg_max_y(1) + seg_min_y(1)
   
   do fn = 2, max_fn
      base_dx = seg_min_x(fn)
      base_dy = seg_min_y(fn)
      
      if (base_dx < current_max_x) then
         shift_x = current_max_x - base_dx + d
         offset_x(fn) = shift_x
         current_max_x = base_dx + seg_max_x(fn) + shift_x
      else
         offset_x(fn) = 0.
         current_max_x = max(current_max_x, base_dx + seg_max_x(fn))
      endif
      
      if (base_dy < current_max_y) then
         shift_y = current_max_y - base_dy + d
         offset_y(fn) = shift_y
         current_max_y = base_dy + seg_max_y(fn) + shift_y
      else
         offset_y(fn) = 0.
         current_max_y = max(current_max_y, base_dy + seg_max_y(fn))
      endif
   end do
 
   print *, "--- Unfolding Offsets ---"
   do fn = 1, max_fn
      print *, "Fault", fn, " Offset X:", offset_x(fn), " Y:", offset_y(fn)
   end do
   print *, "-------------------------"
 
   !=================================================================
   ! Initialize Combined Arrays
   !=================================================================
   allocate(vr_str_all(n_bin, 0:tmax), vr_dip_all(n_bin, 0:tmax))
   vr_str_all = 0. ; vr_dip_all = 0.
   min_valid_x = n_bin; max_valid_x = 1
   min_valid_y = n_bin; max_valid_y = 1
 
   allocate(str_sliprate(0:mn+2, 0:tmax), dip_sliprate(0:nn+2, 0:tmax))
 
   ! --- Info Display ---
   ref_strike = strike0
   ref_dip    = dip0
   bash_file  = "none"
   inquire(file='./run_PDTI_TA.bash', exist=file_exist)
   if (file_exist) then
      bash_file = './run_PDTI_TA.bash'
   else
      inquire(file='./const/run_PDTI_TA.bash', exist=file_exist)
      if (file_exist) then
         bash_file = './const/run_PDTI_TA.bash'
      endif
   endif
   if (trim(bash_file) /= "none") then
      open(60, file=trim(bash_file), status='old', action='read')
      do
         read(60, '(A)', iostat=ierr) line_buf
         if (ierr /= 0) exit
         pos_eq = index(line_buf, 'rstike=')
         if (pos_eq > 0) read(line_buf(pos_eq+7:), *, iostat=ierr) ref_strike
         pos_eq = index(line_buf, 'rstrike=')
         if (pos_eq > 0) read(line_buf(pos_eq+8:), *, iostat=ierr) ref_strike
         pos_eq = index(line_buf, 'rdip=')
         if (pos_eq > 0) read(line_buf(pos_eq+5:), *, iostat=ierr) ref_dip
      end do
      close(60)
      print *, 'Read reference info from: ', trim(bash_file)
   endif
   print *, "Reference Strike:", ref_strike
   print *, "Reference Dip   :", ref_dip
 
   !=================================================================
   ! Main Processing Loop
   !=================================================================
   do fn = 1, max_fn
      write(filename1, '("vr_str_", i0, ".dat")') fn
      write(filename2, '("vr_dip_", i0, ".dat")') fn
      open(11, file=trim(filename1), status='replace')
      open(12, file=trim(filename2), status='replace')
 
      ! Find Index Boundaries for this segment
      min_m = 10000; max_m = -10000
      min_n = 10000; max_n = -10000
      
      do n = 1, nn
         do m = 1, mn
            if (mf_no(m,n) == fn) then
               if (m < min_m) min_m = m
               if (m > max_m) max_m = m
               if (n < min_n) min_n = n
               if (n > max_n) max_n = n
            end if
         end do
      end do
      
      if (min_m > max_m) then
         close(11); close(12)
         cycle
      endif
 
      ! --- Projection (Extract Max Sliprate) ---
      str_sliprate = 0. ; dip_sliprate = 0.
      do t = 0, tmax
         do m = min_m, max_m
            val = 0.
            do n = min_n, max_n
               if (mf_no(m,n) == fn) then
                  if (sliprate(m,n,t) > val) val = sliprate(m,n,t)
               endif
            end do
            str_sliprate(m, t) = val
         end do
         do n = min_n, max_n
            val = 0.
            do m = min_m, max_m
               if (mf_no(m,n) == fn) then
                  if (sliprate(m,n,t) > val) val = sliprate(m,n,t)
               endif
            end do
            dip_sliprate(n, t) = val
         end do
      end do
 
      ! --- Find tmax2 (Backward Scan) ---
      tmax2 = 0
      do t = tmax, 0, -1
        if (sum(str_sliprate(min_m:max_m, t)) > 1.e-6 .or. &
            sum(dip_sliprate(min_n:max_n, t)) > 1.e-6) then
          tmax2 = t 
          exit
        end if
      enddo
      if (tmax2 == 0) tmax2 = tmax 
 
      ! --- Interpolation and Output ---
      do t = 0, tmax2
         
         ! === Strike Interpolation ===
         do m = min_m, max_m - 1
            valid_start = .false.
            valid_end   = .false.
            do n = min_n, max_n
               if (mf_no(m, n) == fn) then
                  x_start = unfolded_dx(m, n) 
                  valid_start = .true.
                  exit
               endif
            end do
            do n = min_n, max_n
               if (mf_no(m+1, n) == fn) then
                  x_end = unfolded_dx(m+1, n) 
                  valid_end = .true.
                  exit
               endif
            end do
 
            if (valid_start .and. valid_end) then
                ! Determine number of steps based on actual unfolded distance
                n_step = int(abs(x_end - x_start) / d)
                if (n_step == 0) n_step = 1
                
                do i = 0, n_step
                   rate = real(i) / real(n_step)
                   loc = x_start * (1.-rate) + x_end * rate
                   val = str_sliprate(m, t) * (1.-rate) + str_sliprate(m+1, t) * rate
                   
                   ! Output Individual
                   loc_phys = loc + seg_min_x(fn)
                   write(11,'(3f15.4)') loc_phys, t*dt, val
                   
                   ! Accumulate Combined
                   loc_comb = loc_phys + offset_x(fn)
                   idx_x = nint(loc_comb / d) + offset_bin
                   if (idx_x >= 1 .and. idx_x <= n_bin) then
                      if (val > vr_str_all(idx_x, t)) vr_str_all(idx_x, t) = val
                      if (idx_x < min_valid_x) min_valid_x = idx_x
                      if (idx_x > max_valid_x) max_valid_x = idx_x
                   endif
                end do
            endif
         end do
         
         ! === Dip Interpolation ===
         do n = min_n, max_n - 1
            valid_start = .false.
            valid_end   = .false.
            do m = min_m, max_m
               if (mf_no(m, n) == fn) then
                  y_start = unfolded_dy(m, n) 
                  valid_start = .true.
                  exit
               endif
            end do
            do m = min_m, max_m
               if (mf_no(m, n+1) == fn) then
                  y_end = unfolded_dy(m, n+1) 
                  valid_end = .true.
                  exit
               endif
            end do
 
            if (valid_start .and. valid_end) then
                n_step = int(abs(y_end - y_start) / d)
                if (n_step == 0) n_step = 1
                
                do i = 0, n_step
                   rate = real(i) / real(n_step)
                   loc = y_start * (1.-rate) + y_end * rate
                   val = dip_sliprate(n, t) * (1.-rate) + dip_sliprate(n+1, t) * rate
                   
                   ! Output Individual
                   loc_phys = loc + seg_min_y(fn)
                   write(12,'(3f15.4)') loc_phys, t*dt, val
                   
                   ! Accumulate Combined
                   loc_comb = loc_phys + offset_y(fn)
                   idx_y = nint(loc_comb / d) + offset_bin
                   if (idx_y >= 1 .and. idx_y <= n_bin) then
                      if (val > vr_dip_all(idx_y, t)) vr_dip_all(idx_y, t) = val
                      if (idx_y < min_valid_y) min_valid_y = idx_y
                      if (idx_y > max_valid_y) max_valid_y = idx_y
                   endif
                end do
            endif
         end do
      end do
      close(11); close(12)
   end do
 
   !=================================================================
   ! PART 3: Output Combined Files
   !=================================================================
   open(11, file='vr_str.dat', status='replace')
   open(12, file='vr_dip.dat', status='replace')
   
   tmax2 = 0
   do t = tmax, 0, -1
     if (sum(vr_str_all(min_valid_x:max_valid_x, t)) > 1.e-6 .or. &
         sum(vr_dip_all(min_valid_y:max_valid_y, t)) > 1.e-6) then
       tmax2 = t
       exit
     endif
   enddo
   if (tmax2 == 0) tmax2 = tmax 
 
   do t = 0, tmax2
      ! Strike
      do i = min_valid_x, max_valid_x
         val = vr_str_all(i, t)
         loc = (real(i) - real(offset_bin)) * d
         write(11, '(3f15.4)') loc, t*dt, val
      end do
      
      ! Dip
      do i = min_valid_y, max_valid_y
         val = vr_dip_all(i, t)
         loc = (real(i) - real(offset_bin)) * d
         write(12, '(3f15.4)') loc, t*dt, val
      end do
   end do
   
   close(11); close(12)
 
   !- - -
 100 stop 
 end program fgenstrdiprup_mf