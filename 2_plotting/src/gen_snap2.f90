!****************************************************************
!
! This program is based on GenerateSnapShot / src / SnapShot.f90
!
!****************************************************************
program fgensnap2
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
  real    :: lat_w, lat_s, lat_n, lat_e, w_depth, s_depth, n_depth, e_depth
  real    :: lon_w, lon_s, lon_n, lon_e, lonmax, lonmin, latmax, latmin
  real, parameter :: rad = 3.1415927/180., dt = 0.1
  real, allocatable :: disp1(:, :, :), mij(:, :), mrf(:, :), mtx(:), mty(:)
  real, allocatable :: xmo(:), slip(:, :), disp2(:, :, :, :), tr(:, :)
  real, allocatable :: potency(:, :, :), sliprate(:, :, :)
  real, allocatable :: stf(:, :, :, :), istf(:, :, :, :)
  integer :: nodal_p
  !--- dummy matrix ---
  real, allocatable :: w(:), w_1(:)
  !-- For FPSPACK.FOR
  integer, parameter :: FPS = 1 ! 1 uses FPSPACK. 0 does not use it.
  real am,ama,am0,am1,e,am0b,slipa,&
       slipb,trendp,plungp,trendt,plungt,trendb,plungb,eta
  integer ierr
  !---
  integer :: dum
  real    :: Mxy(1:3, 1:3), NDC, eigenval(1:3), EV(1:3, 1:3)
  !____
  read(5,*,end=100)twrange, duration, AveSnapFlag
  !____ output file ___
  open(30, file = 'snap.dat', status = 'replace')
  open(31, file = 'snap2.dat', status = 'replace')
  open(32, file='tw_mec_xy.dat', status='replace')
  open(33, file = 'faultline.dat', status = 'replace')
  !-- read fort.40 --
  read(40, *); read(40, *)dm, dm, dm, epilat, epilon, epidep, dm, ns
  read(40, *); read(40, *)strike0, dip0, dm, rake0, dm
  read(40, *); read(40, *)xx, yy, mn, nn, m0, n0, rtime, jtn, icmn, ylen
  close(40)
  allocate(slip(mn, nn), tr(mn, nn))
  allocate(disp1(mn, nn, icmn), disp2(mn, nn, icmn, jtn))
  !- - -
  call readfrom40(mn, nn, icmn, jtn, slip, disp1, tr, disp2)
  !- - -
  ! read(5, *)twrange, duration, AveSnapFlag
  !duration=maxval(tr)+(jtn+2)*rtime+5.
  ! if (twrange >= duration)then
  !    twrange = duration
  ! end if
  if (AveSnapFlag /= 0 .and. AveSnapFlag /= 1) then
     write(6, *)" AveSnapFlag error. [0: cumulative] [1: moment] [else: stop]"
     stop
  end if
  num_timewindow = int( duration / twrange ) !+ 1
  write(6, *)" # of timewindow  :", num_timewindow
  !- - -
  allocate(mtx(0:mn+1), mty(0:nn+1))
  call get_xy_location(xx, yy, mn, nn, m0, n0, mtx, mty)
  !- - - - - - - - - - - - - - - - - - - - - -
  !  get source time function at every knot   |
  !- - - - - - - - - - - - - - - - - - - - - -
  tmax = int(duration / dt) + 1 !  dt = 0.1
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !    stf: slip-rate function of m1 ~ m5 at every knot
  ! w, w_1: work space
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  w = 0.
  allocate(istf(0:mn+2, 0:nn+2, 1:num_timewindow, 1:icmn))
  allocate(potency(0:mn+2, 0:nn+2, 1:num_timewindow), sliprate(0:mn+2, 0:nn+2, 1:num_timewindow))
  allocate(w_1(1:3)) ! re-used dummy vector
  
  !   ////////////////////////////////////////////////////////////////////////////////////////////////////////
  !   snap.datの生成
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! istf: slip-rate function of m1 ~ m5 for each time window at every knot
  ! potency  : potency-rate function at every knot
  ! sliprate : slip-rate function at every knot
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  istf = 0. ;  potency = 0. ; sliprate = 0.
  !- -
  do n = 1, nn
     do m = 1, mn
        do tw = 1, num_timewindow
           allocate(xmo(1:6), mij(1:3, 1:3), mrf(1:3, 1:3))
           xmo = 0.
           do icm = 1, icmn
              num_t_in_tw = int(twrange / dt)
              cum_t_bf_tw = int(((tw - 1) * twrange) / dt)
              tcount = 0
              !==== average value in each time window ====
              if (AveSnapFlag == 0)then
                 do t = cum_t_bf_tw + 1, cum_t_bf_tw + num_t_in_tw
                    if (t >= tmax) exit
                    istf(m, n, tw, icm) = istf(m, n, tw, icm) + stf(m, n, t, icm)
                    tcount = tcount + 1
                 end do
                 istf(m, n, tw, icm) = istf(m, n, tw, icm) / tcount
              !==== value at each time step ====
              elseif (AveSnapFlag == 1)then
                 istf(m, n, tw, icm) = stf(m, n, cum_t_bf_tw + num_t_in_tw, icm)
              end if
              !=================================
              xmo(icm) = istf(m, n, tw, icm)
           end do
           call mtrx(xmo, mij)
           if( icmn .ne. 5 )call fslip2mij(xmo,icmn,strike0,dip0,rake0, mij)
           call meca_rotation(mij, rad) ! add for meca rotation
           call mxy_mrf(mij, mrf)
           !- - - - - - - - - - - - - - - - - - - - - -
           if (icmn == 5)then
              call d_cp(xmo, strike, dip, rake, sliprate(m, n, tw), dm, w_1)
              call conj(strike, dip, rake, a_strike, a_dip, a_rake)
           else
              call d_cp(xmo, strike, dip, rake, dm, dm, w_1)
              call conj(strike, dip, rake, a_strike, a_dip, a_rake)
              !- - - - - - - - - - - - - - - - - - - - - -
              sliprate(m, n, tw) = sqrt(xmo(1)**2 + xmo(2)**2)
           end if
           if ( FPS .eq. 1 ) then ! add for FPSPACK.FOR
             call mtrx(xmo, mij)
             call AR2PLP(mij, am0,am1,e,am0b,strike,dip,rake,slipa,&
                         a_strike,a_dip,a_rake,slipb,&
                         trendp,plungp,trendt,plungt,trendb,plungb,eta,ierr)
           end if
           if (ns /= 0 .and. n == nn) then
              potency(m, n, tw) = sliprate(m, n, tw) * xx * (yy + ylen) * 1.e6
           else
              potency(m, n, tw) = sliprate(m, n, tw) * xx * yy * 1.e6
           endif
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           nexpo = int(log10(maxval(abs(mrf))))
            !           mrf = mrf / (10.**nexpo)
           !- -
            !     if (sliprate(m, n, tw) .ne. 0)then
           if (sliprate(m, n, tw) .gt. 0.0000001)then
              !---
              call mrf_mxy(Mrf, Mxy)
              call EIG1(Mxy, 3, 3, eigenval, EV, dum)
              ! based on Dziewonski+1981JGR
              NDC = - eigenval(2) / max(abs(eigenval(1)), abs(eigenval(3)))
              NDC = 100 * NDC / 0.5 ! convert to parcentage
              !---
              call xy2geo(epilat, epilon, mtx(m), mty(n), strike0, dip0, lat, lon)
              dep = epidep - mty(n) * sin(dip0 * rad)
              call selection_nodal (strike, dip, a_strike, a_dip, strike0, dip0, rad, nodal_p)
              call swap_fault( strike, dip, rake, a_strike, a_dip, a_rake, nodal_p) 
              write(30, '(3i6, 2f9.2, 7f16.9, 4f8.2, 2f9.2, 3f9.3, 6f8.2, f9.2)')n, m, tw, mtx(m), mty(n), &
                   sliprate(m, n, tw), mrf(1, 1), mrf(2, 2), mrf(3, 3), mrf(1, 2), &
                   mrf(1, 3), mrf(2, 3), strike, a_strike, dip, a_dip, rake, a_rake, &
                   lat, lon, dep, &
                   trendp,trendt,trendb,plungp,plungt,plungb,NDC
           end if
           deallocate(xmo, mij, mrf)
        end do
     end do
  end do
  close(30)
  
  !---
  call faultline(strike0,dip0,epilat,epilon,epidep,xx,yy,mn,nn,m0,n0,1.0)
  !---
  stop
100 stop "failed to read arguments. stop."
end program
