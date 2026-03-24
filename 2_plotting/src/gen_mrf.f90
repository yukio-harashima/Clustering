!****************************************************************
!
! This program is based on Snippets / src / plotslip.f90
!
!****************************************************************
program fgenmrf
  implicit none
  ! real, parameter :: dt=0.1
  real              :: dt
  real, allocatable :: disp(:,:), disp1(:,:,:), disp2(:,:,:,:),&
                       tr(:,:), slip(:,:)
  real, allocatable :: rigidd(:,:), amp_sf(:,:)
  real, allocatable :: Ro(:), R3(:,:,:,:), Ro1(:), Ro2(:)
  real, allocatable :: xmo(:),mij(:,:),mrf(:,:)
  real, allocatable :: st(:), potencyrate(:)
  integer :: nsurface, mn, nn, m0, n0, jtn, icmn
  real :: moment, Mw, rigid, lat, lon, depth, vr, &
          strike, dip, slip0, slipi, rslip, &
          xx, yy, rtime, ylength, &
          varicance, ABIC
  integer :: i, n, m, jt, icm, imax, imax2
  real :: tmax, shift, SM, DSM, wk, dmm
  real :: AX(3)
  !---
  open(12,file='mrf.dat',status='replace')
  !---
  read(40,*); read(40,*)moment,Mw,rigid ,lat,lon,depth,vr,nsurface
  read(40,*); read(40,*)strike,dip,slip0,slipi,rslip
  read(40,*); read(40,*)xx,yy,mn,nn,m0,n0,rtime,jtn,icmn,ylength
  read(40,*); read(40,*)varicance, ABIC
  read(40,*)
  allocate (disp(0:mn+2,0:nn+2),disp1(0:mn+2,0:nn+2,1:icmn), &
       disp2(1:mn,1:nn,1:jtn,1:icmn),tr(1:mn,1:nn),slip(0:mn+1,0:nn+1), &
       rigidd(1:mn,1:nn),amp_sf(1:mn,1:nn))
  disp = 0.; disp1 = 0.; slip = 0.
  do n=nn,1,-1 ; read(40,*)(disp(m,n),m=1,mn) ; end do
  read(40,*)
  do n=nn,1,-1 ; read(40,*)(slip(m,n),m=1,mn) ; end do
  read(40,*)
  do n=nn,1,-1 ; read(40,*)(tr(m,n),m=1,mn)   ; end do
  read(40,*)
  do icm=1,icmn ; read(40,*)
    do n=nn,1,-1 ; read(40,*)(disp1(m,n,icm),m=1,mn) ; end do
  end do
  do jt=1,jtn ; do icm=1,icmn ;read(40,*)
    do n=nn,1,-1 ; read(40,*)(disp2(m,n,jt,icm),m=1,mn) ; end do
  end do ;  end do
  !---
  !
  ! Prevent resampling from making the nonlinearity of seismic moment
  ! derived from eigenvalue more pronounced.
  ! Yagi & Yamashita 2021/12/14
  !
  dt = rtime
  !---
  open(41,file='rigid_amp.info',status='old',err=999)
  do i = 1, mn*nn
    read(41,*,end=999)n,m,rigidd(m,n), amp_sf(m,n)
  end do
  goto 998
999 write(6,*) 'failed to read rigid_amp.info'
  rigidd = rigid
  amp_sf = 1.
  do n=1,nn
    rigidd(m,n)=rigid
  end do
998 close(41)
  !---
  tmax = maxval(tr) + (jtn+2) * rtime  + 5.
  imax = nint(tmax / dt)
  allocate (st(imax),potencyrate(1:imax))
  allocate (Ro(imax),R3(0:mn+2,0:nn+2,imax,icmn),Ro1(imax),Ro2(imax))
  Ro=0.; R3=0.; Ro1=0.; Ro2=0.
  call stime1(st,imax,dt,rtime,rtime*2.)
  do icm=1,icmn
    do m=1,mn
      do n=1,nn
        do jt=1,jtn
          shift = rtime * (jt-1) + tr(m,n)
          call resample_shift(st, imax, dt, Ro, imax, dt, shift)
          R3(m,n,1:imax,icm) = R3(m,n,1:imax,icm) +  disp2(m,n,jt,icm)*Ro(1:imax)
        end do
      end do
    end do
  end do
  deallocate(st, disp2)
  !---
  Ro=0.; Ro1=0.; Ro2=0.; potencyrate=0.
  if (icmn.eq.2) then
    do i=1,imax
      do m=1,mn
        do n=1,nn
          Ro1(i) = rigidd(m,n) * amp_sf(m,n) * R3(m,n,i,1) + Ro1(i)
          Ro2(i) = rigidd(m,n) * amp_sf(m,n) * R3(m,n,i,2) + Ro2(i)
        end do
      end do
      Ro(i) = sqrt(Ro1(i)**2 + Ro2(i)**2)
    end do
  else if (icmn.eq.5) then
    do i=1, imax
      allocate(xmo(1:6), mij(1:3, 1:3), mrf(1:3, 1:3))
      xmo=0.
      do icm=1,icmn
        do m=1,mn
          do n=1,nn
            xmo(icm) = xmo(icm) + R3(m,n,i,icm) * rigidd(m,n) * amp_sf(m,n)
          end do
        end do
      end do
      call d_cp(xmo,dmm,dmm,dmm,SM,DSM,AX)
      Ro(i) = SM
      deallocate(xmo,mij,mrf)
      allocate(xmo(1:6), mij(1:3, 1:3), mrf(1:3, 1:3))
      xmo=0.
      do m=1,mn
        do n=1,nn
          do icm=1,icmn
            xmo(icm) = xmo(icm) + R3(m,n,i,icm) * amp_sf(m,n)
          end do
        end do
      end do
      call d_cp(xmo,dmm,dmm,dmm,SM,DSM,AX)
      potencyrate(i) = SM
      deallocate(xmo,mij,mrf)
    end do
  else
    write(6, *)" !!!STOP!!! ICMN is not 2 or 5"
    stop
  end if
  !---
  wk = xx * yy * (10.**15.) / (10.**18.)
  Ro = Ro * wk
  !------
  do i = 1, imax-3
     if( sum(Ro(i:i+3)) .eq. 0. ) then
        imax2 = i
        exit
     end if
  end do
!  write(6,*) imax, imax2
  if (icmn.eq.2) then
    write(12,'(f8.3,e15.4)')(real(i-1)*dt,Ro(i),i=1,imax2)
  else if (icmn.eq.5) then
    do i=1,imax2
      allocate(xmo(1:6), mij(1:3, 1:3), mrf(1:3, 1:3))
      xmo=0.; mij=0.; mrf=0.
      do m=1,mn
        do n=1,nn
          do icm=1,icmn
            xmo(icm) = xmo(icm) + R3(m, n, i, icm) * rigidd(m,n) * amp_sf(m,n)
          end do
        end do
      end do
      call mtrx(xmo, mij)
      call mxy_mrf(mij, mrf)
      write(12,'(f8.3, e15.4, 6f15.6)') real(i-1) * dt, Ro(i), &
           mrf(1, 1), mrf(2, 2), mrf(3, 3), mrf(1, 2), mrf(1, 3), mrf(2, 3)
      deallocate(xmo, mij, mrf)
    end do
  end if
  !---
  close(12)
end program
