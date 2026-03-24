!****************************************************************
!
! Program: fgenpdtdis_mf
! Description: 
!   Generates detailed fault analysis data (pdtdis.dat) and 
!   interpolated slip distribution (pddis_*.dat) for multi-segment 
!   fault models defined in Multiple_fault.dat.
!
! Input: 
!   - fort.40 (Slip distribution data)
!   - Multiple_fault.dat (Multi-segment geometry & parameters)
!
! Output:
!   - pdtdis.dat (Detailed mechanism data)
!   - pddis_{No}.dat (Interpolated slip for each segment)
!   - faultline.dat (Boundaries of fault segments)
!
!　Compile
!　 gfortran gen_pdtdis_multi.f90 sub.gen.f90 -o fgenpdtdis_mf
!****************************************************************
program fgenpdtdis_mf
    implicit none
  
    ! --- Grid and Dimensions ---
    integer :: m, n, mn, nn, icmn, m0, n0, jtn, i, ns
    real :: ylen, rtime, rake
    real :: strike, dip, xx, yy, epilat, epilon, epidepth
    real, parameter   :: dx = 0.5d0, rad = 3.1415927/180.
  
    ! --- Standard Arrays ---
    real, allocatable :: lat(:, :), lon(:, :), depth(:, :)
    real, allocatable :: slip(:, :), tr(:, :)
    real, allocatable :: i_slip(:, :, :), sliprate(:, :, :, :)
    
    ! --- Mechanism Calculation Arrays ---
    real, allocatable :: mij(:, :), tmp_mrf(:, :), xmo(:), potency(:, :), mrf(:, :, :, :)
    real, allocatable :: str1(:, :), dip1(:, :), rake1(:, :) 
    real, allocatable :: str2(:, :), dip2(:, :), rake2(:, :) 
    integer :: nodal_p
    real :: ax(3), dm
    
    ! --- Multiple Fault Data Arrays ---
    real, allocatable :: mf_dx(:,:), mf_dy(:,:)
    real, allocatable :: mf_lat(:,:), mf_lon(:,:), mf_depth(:,:)
    real, allocatable :: mf_stk(:,:), mf_dip(:,:), mf_rake(:,:)
    real, allocatable :: mf_atime(:,:)
    integer, allocatable :: mf_no(:,:)
    
    ! --- Loop & Temp Variables for Reading ---
    integer :: fn, max_fn, tm, tn, ierr, tno
    real :: tdx, tdy, tlat, tlon, tdep, tstk, tdip, trake, tatime
    character(len=30) :: filename
  
    ! --- Interpolation Variables (pddis) ---
    integer :: nx, ny, vm0, vn0, ref_m_idx, ref_n_idx
    real    :: point1, point2, tmp_potency, lat_tmp, lon_tmp, slip_xy, tmp_depth
    real    :: local_dx_min, local_dx_max, local_dy_min, local_dy_max
    real    :: r_dx_min, r_dy_min
    real    :: ref_lat, ref_lon, ref_dx, ref_dy, ref_stk, ref_dip
  
    ! --- Faultline Variables ---
    integer :: min_m, max_m, min_n, max_n
    logical :: found_plane
  
    ! --- FPSPACK Variables ---
    integer, parameter :: FPS = 1 ! 1 uses FPSPACK
    real :: am0, am1, e, am0b, slipa, slipb, eta
    real, allocatable :: trendp(:, :), plungp(:, :), trendt(:, :)
    real, allocatable :: plungt(:, :), trendb(:, :), plungb(:, :)
    integer :: dum
    real    :: Mxy(1:3, 1:3), NDC, eigenval(1:3), EV(1:3, 1:3)
  
    !=================================================================
    ! 1. fort.40 から基本パラメータの読み込み
    !=================================================================
    open(40, file='fort.40', status='old')
    read(40, *); read(40, *)dm, dm, dm, epilat, epilon, epidepth, dm, ns
    read(40, *); read(40, *)strike, dip, dm, rake, dm
    read(40, *); read(40, *)xx, yy, mn, nn, m0, n0, rtime, jtn, icmn, ylen
    close(40)
  
    ! --- メモリ確保 ---
    allocate(lat(1:mn, 1:nn), lon(1:mn, 1:nn), depth(1:mn, 1:nn))
    allocate(slip(1:mn, 1:nn), tr(1:mn, 1:nn), i_slip(1:mn, 1:nn, 1:icmn))
    allocate(sliprate(1:mn, 1:nn, 1:icmn, 1:jtn), mrf(1:mn, 1:nn, 1:3, 1:3))
    allocate(str1(1:mn, 1:nn), str2(1:mn, 1:nn), dip1(1:mn, 1:nn), dip2(1:mn, 1:nn))
    allocate(rake1(1:mn, 1:nn), rake2(1:mn, 1:nn), potency(0:mn+2, 0:nn+2))
    
    allocate(trendp(1:mn, 1:nn), plungp(1:mn, 1:nn), trendt(1:mn, 1:nn))
    allocate(plungt(1:mn, 1:nn), trendb(1:mn, 1:nn), plungb(1:mn, 1:nn))
  
    ! --- Multiple Fault用配列の確保 ---
    allocate(mf_dx(1:mn, 1:nn), mf_dy(1:mn, 1:nn))
    allocate(mf_lat(1:mn, 1:nn), mf_lon(1:mn, 1:nn), mf_depth(1:mn, 1:nn))
    allocate(mf_stk(1:mn, 1:nn), mf_dip(1:mn, 1:nn), mf_rake(1:mn, 1:nn))
    allocate(mf_atime(1:mn, 1:nn), mf_no(1:mn, 1:nn))
  
    ! mf_no を -1 で初期化 (データ無し判定用)
    mf_no = -1 
    potency = 0. ; mrf = 0.
  
    ! 詳細なすべりデータの読み込み (readfrom40)
    call readfrom40(mn, nn, icmn, jtn, slip, i_slip, tr, sliprate)
  
    !=================================================================
    ! 2. Multiple_fault.dat の読み込みとマッピング
    !=================================================================
    open(50, file='Multiple_fault.dat', status='old')
    
    ! ヘッダーがある場合は読み飛ばす (必要に応じて有効化)
    ! read(50, *, iostat=ierr) 
  
    max_fn = 0
    do
       ! 読み込み: n, m, dx, dy, lat, lon, depth, strike, dip, rake, time, fault_no
       read(50, *, iostat=ierr) tn, tm, tdy, tdx, tlat, tlon, tdep, tstk, tdip, trake, tatime, tno
       
       if (ierr /= 0) exit
       
       ! m, n のインデックスに基づいて配列にマッピング (ファイルの行順序に依存しない)
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
  
    !=================================================================
    ! 3. pdtdis.dat の生成 (平面ごとのパラメータを使用)
    !=================================================================
    open(42, file='pdtdis.dat', status='replace')
  
    do n=1,nn
      do m=1,mn
        ! Multiple_fault.datにデータが存在するグリッドのみ処理
        if (mf_no(m, n) /= -1) then
        
          allocate(xmo(1:6), mij(1:3, 1:3), tmp_mrf(1:3, 1:3))
          xmo = 0.
          do i = 1, icmn
            xmo(i) = i_slip(m, n, i)
          enddo
  
          ! --- モーメントテンソルの計算 ---
          call mtrx(xmo, mij)
          
          ! 局所的な断層パラメータ (mf_stk, mf_dip, mf_rake) を使用
          if( icmn .ne. 5 ) call fslip2mij(xmo, icmn, mf_stk(m,n), mf_dip(m,n), mf_rake(m,n), mij)
          
          call meca_rotation(mij, rad)
          call mxy_mrf(mij, tmp_mrf)
          
          ! --- 断層メカニズム解の導出 (d_cp) ---
          call d_cp(xmo, str1(m, n), dip1(m, n), rake1(m, n), potency(m, n), dm, ax)
          call conj(str1(m, n), dip1(m, n), rake1(m, n), str2(m, n), dip2(m, n), rake2(m, n))
  
          ! --- P/T/B軸の計算 (FPSPACK) ---
          if ( FPS .eq. 1 ) then
             call mtrx(xmo, mij)
             call AR2PLP(mij, am0,am1,e,am0b,str1(m, n),dip1(m, n),rake1(m, n),slipa,&
                  str2(m, n),dip2(m, n),rake2(m, n),slipb,&
                  trendp(m, n),plungp(m, n),trendt(m, n),plungt(m, n),&
                  trendb(m, n),plungb(m, n),eta,ierr)
          end if
  
          if (icmn.eq.2) potency(m, n) = sqrt(xmo(1)**2 + xmo(2)**2)
          mrf(m, n, 1:3, 1:3) = tmp_mrf(1:3, 1:3)
          deallocate(xmo, mij, tmp_mrf)
  
          ! --- pdtdis.dat への出力 ---
          if (potency(m,n) > 0.) then
             ! NDC成分の計算
             call mrf_mxy(Mrf(m, n, 1:3, 1:3), Mxy)
             call EIG1(Mxy, 3, 3, eigenval, EV, dum)
             NDC = - eigenval(2) / max(abs(eigenval(1)), abs(eigenval(3)))
             NDC = 100 * NDC / 0.5
  
             ! 局所的な走向・傾斜に基づいて節面を選択
             call selection_nodal (str1(m,n), dip1(m,n), str2(m,n), dip2(m,n), &
                                   mf_stk(m,n), mf_dip(m,n), rad ,nodal_p )
             call swap_fault (str1(m,n), dip1(m,n), rake1(m,n), str2(m,n), dip2(m,n), rake2(m,n) ,nodal_p )
  
             ! Multiple_fault.datの座標を使用して書き出し
             write(42, '(2i5, 2f9.2, 2f12.5, f10.3, f11.4, 6f16.9, i5, 4f8.2, 2f9.2, 6f8.2, f9.2)') &
                 n, m, &
                 mf_dy(m, n), mf_dx(m, n), &                   ! dx, dy
                 mf_lat(m, n), mf_lon(m, n), mf_depth(m, n), & ! lat, lon, depth
                 potency(m, n), &
                 mrf(m, n, 1, 1), mrf(m, n, 2, 2), mrf(m, n, 3, 3), &
                 mrf(m, n, 1, 2), mrf(m, n, 1, 3), mrf(m, n, 2, 3), &
                 20, &
                 str1(m, n), str2(m, n), dip1(m, n), dip2(m, n), rake1(m, n), rake2(m, n), &
                 trendp(m, n), trendt(m, n), trendb(m, n), &
                 plungp(m, n), plungt(m, n), plungb(m, n), NDC
          end if
        end if
      end do
    end do
    close(42)
  
    !=================================================================
    ! 4. pddis_{fn}.dat の生成 (平面ごとの補間処理)
    !=================================================================
    
    do fn = 1, max_fn
      
        ! ファイル名の生成 (例: pddis_1.dat)
        write(filename, '("pddis_", i0, ".dat")') fn
        open(43, file=trim(filename), status='replace')
    
        ! この断層面の範囲と基準点の探索
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
                 
                 ! 基準点の取得 (最初に見つかった有効な点)
                 if (ref_lat == -999.) then
                    ref_lat = mf_lat(m,n)
                    ref_lon = mf_lon(m,n)
                    ref_dx  = mf_dx(m,n)
                    ref_dy  = mf_dy(m,n)
                    ref_stk = mf_stk(m,n)
                    ref_dip = mf_dip(m,n)
                    ! [修正] 仮想原点計算用にインデックスを保存
                    ref_m_idx = m
                    ref_n_idx = n
                 endif
              endif
           enddo
        enddo
        
        ! このfnに対応する点が見つからない場合はスキップ
        if (ref_lat == -999.) then
            close(43)
            cycle
        endif
    
        ! [修正] 仮想原点 (Virtual Origin) の計算
        ! slip_xy が正しいインデックスを参照できるようにオフセットを逆算
        vm0 = ref_m_idx - nint(ref_dx / xx)
        vn0 = ref_n_idx - nint(ref_dy / yy)
  
        ! [修正] 丸め処理 (NINT & 0.1km刻み) と 開始点の決定
        ! データの最小値から1グリッド分(xx)戻った位置を基準とする
        r_dx_min = real(nint((local_dx_min - xx) * 10.)) / 10.
        r_dy_min = real(nint((local_dy_min - yy) * 10.)) / 10.
  
        ! [修正] 補間用グリッドサイズの決定
        ! 丸められた開始点から、最大値+1グリッド分(xx)までの範囲をカバー
        ! 安全マージン + 1.0 (fgensnap_y_mfと統一)
        nx = int((local_dx_max + xx - r_dx_min) * 2. + 1.0)
        ny = int((local_dy_max + yy - r_dy_min) * 2. + 1.0)
        
        ! 補間ループ
        do n = 1, ny
           do m = 1, nx
              ! [修正] 丸め済みの基準点からステップ移動
              point1 = r_dx_min + dx * (m - 1)
              point2 = r_dy_min + dx * (n - 1)
              
              ! [修正] 仮想原点 vm0, vn0 を使用して補間
              tmp_potency = slip_xy(potency(0, 0), mn, nn, vm0, vn0, xx, yy, point1, point2)
              
              ! [修正] 出力条件: 0.0 も含める (>= 0.)
              if (tmp_potency >= 0.) then
                 ! 基準点からの相対距離を用いて地理座標を計算
                 ! Delta X = point1 - ref_dx
                 ! Delta Y = point2 - ref_dy
                 call xy2geo(ref_lat, ref_lon, point1 - ref_dx, point2 - ref_dy, &
                             ref_stk, ref_dip, lat_tmp, lon_tmp)
                 
                 ! 深さの計算 (基準dipを使用)
                 ! 注意: point2 は深さ方向(down-dip)に増加する
                 tmp_depth = epidepth - (point2 * sin(ref_dip * rad))
                 
                 write(43, '(4f9.3, f13.6, f9.3)') lat_tmp, lon_tmp, point1, point2, &
                      tmp_potency, tmp_depth
              end if
           end do
        end do
        close(43)
      end do
  
    !=================================================================
    ! 5. faultline.dat の生成 (断層面ごとの外形線)
    !=================================================================
    open(44, file='faultline.dat', status='replace')
  
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
          ! 四隅の座標を出力 (左下 -> 右下 -> 右上 -> 左上 -> 左下)
          write(44, '(3f12.4)') mf_lat(min_m, min_n), mf_lon(min_m, min_n), mf_depth(min_m, min_n)
          write(44, '(3f12.4)') mf_lat(max_m, min_n), mf_lon(max_m, min_n), mf_depth(max_m, min_n)
          write(44, '(3f12.4)') mf_lat(max_m, max_n), mf_lon(max_m, max_n), mf_depth(max_m, max_n)
          write(44, '(3f12.4)') mf_lat(min_m, max_n), mf_lon(min_m, max_n), mf_depth(min_m, max_n)
          write(44, '(3f12.4)') mf_lat(min_m, min_n), mf_lon(min_m, min_n), mf_depth(min_m, min_n)
          
          ! GMT/Gnuplot用のセパレータ
          write(44, *) '>' 
        !   write(44, *) 
       end if
    end do
  
    close(44)
  
  end program fgenpdtdis_mf