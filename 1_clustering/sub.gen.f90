!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine readfrom40(mn, nn, icmn, jtn, slip, ic_slip, tr, sliprate)
  implicit none
  integer, intent(in)  :: mn, nn, icmn, jtn
  real,    intent(out) :: slip(1:mn, 1:nn), ic_slip(1:mn, 1:nn, 1:icmn)
  real,    intent(out) :: tr(1:mn, 1:nn), sliprate(1:mn, 1:nn, 1:icmn, 1:jtn)
  integer :: m, n, icm, jt
  !- - - -
  open(430, file = 'fort.40', status = 'old')
  do n = 1, 9
     read(430, *)
  enddo
  do n = nn, 1, -1
     if (mn > 30)then
        read(430, *)(slip(m, n), m = 1, 30)
        read(430, *)(slip(m, n), m = 31, mn)
     else
        read(430, *)(slip(m, n), m = 1, mn)
     endif
  enddo
  read(430, *) !slip angel for each sub-fault
  if (mn > 30)then
     do n = 1, 2 * nn
        read(430, *)
     enddo
  else
     do n = 1,  nn
        read(430, *)
     enddo
  endif
  read(430, *) !start time for each sub-fault
  do n = nn, 1, -1
     if (mn > 30)then
        read(430, *)(tr(m, n), m = 1, 30)
        read(430, *)(tr(m, n), m = 31, mn)
     else
        read(430, *)(tr(m, n), m = 1, mn)
     endif
  enddo
  read(430, *) !total slip vecoter for each sub-fault
  do icm = 1, icmn
     read(430, *)
     do n = nn, 1, -1
        if (mn > 30)then
           read(430, *)(ic_slip(m, n, icm), m = 1, 30)
           read(430, *)(ic_slip(m, n, icm), m = 31, mn)
        else
           read(430, *)(ic_slip(m, n, icm), m = 1, mn)
        endif
     enddo
  enddo

  do jt = 1, jtn
     do icm = 1, icmn
        read(430, *)
        do n = nn, 1, -1
           if (mn > 30)then
              read(430, *)(sliprate(m, n, icm, jt), m = 1, 30)
              read(430, *)(sliprate(m, n, icm, jt), m = 31, mn)
           else
              read(430, *)(sliprate(m, n, icm, jt), m = 1, mn)
           endif
        enddo
     enddo
  enddo
  close(430)

end subroutine readfrom40

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! it is calcurate dip and strike.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine d_cp(v,f1,d1,a1,sm,dsm,ax)
  real v(6),ax(3),eg1(3),eg2(3),ev(3,3),mij(3,3)
  pi = 3.1415926
  rad = pi/180
  root2 = 1.41421356
  call mtrx(v, mij)
  call eig1(mij, 3, 3, eg1, ev, icon)
  if(icon /= 0) print *, '  seig1.cond= ', icon
  e1 = eg1(1)
  e3 = eg1(1)
  i1 = 1
  i3 = 1
  do ii = 2, 3
     if(e1 >= eg1(ii)) goto 1
     e1 = eg1(ii)
     i1 = ii
1    if(e3 < eg1(ii)) cycle
     e3 = eg1(ii)
     i3 = ii
  end do
  sm = (e1-e3) / 2
  dsm = sm - e1
  do ii = 1, 3
     ax(ii) = ev(ii,i1)
     eg1(ii) = (ev(ii, i3)-ev(ii, i1)) / root2
     eg2(ii) = (ev(ii, i1)+ev(ii, i3)) / root2
  end do
  d1 = acos(eg2(3))
  f1 = atan2(eg2(2), eg2(1)) + pi / 2
  cf = cos(f1)
  sf = sin(f1)
  c = eg1(1) * cf + eg1(2) * sf
  if(d1 /= 0.) s = eg1(3) / sin(d1)
  if(d1 == 0.) s = eg1(1) * sf - eg1(2) * cf
  f1 = f1 / rad
  d1 = d1 / rad
  a1 = -atan2(s,c) / rad
  if(d1 < 90.) return
  a1 = -a1
  d1 = 180 - d1
  f1 = f1 + 180
  if(f1 > 360) f1 = f1 - 360
end subroutine d_cp

subroutine eig1(a, k, n, e, u, icon)
  real a(k,n),u(k,n),e(n)
  call jacobi(a, n, k, e, u, m)
  call eigsrt(e, u, n, k)
  icon=0
end subroutine eig1

subroutine jacobi(a, n, np, d, v, nrot)
  parameter (nmax = 100)
  dimension a(np, np), d(np), v(np, np), b(nmax), z(nmax)
  do ip = 1, n
     do iq = 1, n
        v(ip, iq) = 0.
     end do
     v(ip, ip) = 1.
  end do
  do ip = 1, n
     b(ip) = a(ip, ip)
     d(ip) = b(ip)
     z(ip) = 0.
  end do
  nrot = 0
  do i = 1, 50
     sm = 0.
     do ip = 1, n-1
        do iq = ip + 1, n
           sm = sm + abs(a(ip, iq))
        end do
     end do
     if(sm == 0.)return
     if(i < 4)then
        tresh = 0.2 * sm / n **2
     else
        tresh = 0.
     endif
     do ip = 1, n - 1
        do iq = ip + 1, n
           g = 100. * abs(a(ip, iq))
           if((i > 4).and.(abs(d(ip)) + g == abs(d(ip)))             &
                &              .and.(abs(d(iq)) + g == abs(d(iq))))then
              a(ip, iq)=0.
           else if(abs(a(ip, iq)) > tresh)then
              h = d(iq) - d(ip)
              if(abs(h) + g == abs(h))then
                 t = a(ip, iq) / h
              else
                 theta = 0.5 * h / a(ip, iq)
                 t = 1. / (abs(theta) + sqrt(1. + theta **2))
                 if(theta < 0.)t = -t
              endif
              c = 1. / sqrt(1 + t ** 2)
              s = t * c
              tau = s / (1.+c)
              h = t * a(ip, iq)
              z(ip) = z(ip) - h
              z(iq) = z(iq) + h
              d(ip) = d(ip) - h
              d(iq) = d(iq) + h
              a(ip,iq) = 0.
              do j = 1, ip-1
                 g = a(j, ip)
                 h = a(j, iq)
                 a(j, ip) = g - s * (h + g * tau)
                 a(j, iq) = h + s * (g - h * tau)
              end do
              do j = ip + 1, iq - 1
                 g = a(ip, j)
                 h = a(j, iq)
                 a(ip, j) = g - s * (h + g * tau)
                 a(j, iq) = h + s * (g - h * tau)
              end do
              do j = iq + 1, n
                 g = a(ip, j)
                 h = a(iq, j)
                 a(ip, j) = g - s * (h + g * tau)
                 a(iq, j) = h + s * (g - h * tau)
              end do
              do j = 1, n
                 g = v(j, ip)
                 h = v(j, iq)
                 v(j, ip) = g - s * (h + g * tau)
                 v(j, iq) = h + s * (g - h * tau)
              end do
              nrot = nrot + 1
           endif
        end do
     end do
     do ip = 1, n
        b(ip) = b(ip) + z(ip)
        d(ip) = b(ip)
        z(ip) = 0.
     end do
  end do
  return
end subroutine jacobi

subroutine eigsrt(d, v, n, np)
  dimension d(np), v(np, np)
  do i = 1, n - 1
     k = i
     p = d(i)
     do j = i + 1, n
        if(d(j) >=  p)then
           k = j
           p = d(j)
        endif
     end do
     if(k /=  i)then
        d(k) = d(i)
        d(i) = p
        do j = 1, n
           p = v(j, i)
           v(j, i) = v(j, k)
           v(j, k) = p
        end do
     endif
  end do
  return
end subroutine eigsrt

subroutine conj(f1, d1, a1, f2, d2, r2)
  parameter(pi = 3.1415926)
  rad = pi/180.
  if(f1 == 0..and.d1 == 0..and.a1 == 0.) stop
  rf1 = rad * f1
  rd1 = rad * d1
  ra1 = rad * a1
  if(a1 >=  0.) then
     rd2 = acos(sin(rd1) * sin(ra1))
     rf2 = rf1 + atan2(cos(ra1), sin(ra1) * cos(rd1)) + pi
     ra2 = pi - acos(cos(ra1) * sin(rd1)/sin(rd2))
  else
     rd2 = pi - acos(sin(rd1) * sin(ra1))
     rf2 = rf1 + atan2(cos(ra1), sin(ra1) * cos(rd1))
     ra2 = - pi + acos(cos(ra1) * sin(rd1) / sin(rd2))
  end if
  if(rf2 > 2 * pi) rf2 = rf2 - 2 * pi
  f2 = rf2 / rad
  d2 = rd2 / rad
  r2 = ra2 / rad
  if(f2 < 0.) f2 = 360. + f2
  return
end subroutine conj

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine mtrx(v,m)
  real v(6),m(3,3)
  m(1,1) = v(2) - v(5) + v(6)
  m(1,2) = v(1)
  m(2,1) = m(1,2)
  m(1,3) = v(4)
  m(3,1) = m(1,3)
  m(2,2) = -v(2) + v(6)
  m(2,3) = v(3)
  m(3,2) = m(2,3)
  m(3,3) = v(5) + v(6)
end subroutine mtrx

subroutine mxy_mrf(mxy, mrf)
  real mxy(3, 3), mrf(3, 3)
  mrf(1, 1) =  mxy(3, 3)
  mrf(1, 2) =  mxy(3, 1)
  mrf(1, 3) = - mxy(3, 2)
  mrf(2, 1) =  mxy(1, 3)
  mrf(2, 2) =  mxy(1, 1)
  mrf(2, 3) = - mxy(1, 2)
  mrf(3, 1) = - mxy(2, 3)
  mrf(3, 2) = - mxy(2, 1)
  mrf(3, 3) =  mxy(2, 2)
  return
end subroutine mxy_mrf

subroutine mrf_mxy(mrf, mxy)
  implicit none
  real, intent(in)  :: mrf(1:3, 1:3)
  real, intent(out) :: mxy(1:3, 1:3)
  mxy(3, 3) =   mrf(1, 1)
  mxy(3, 1) =   mrf(1, 2)
  mxy(3, 2) = - mrf(1, 3)
  mxy(1, 3) =   mrf(2, 1)
  mxy(1, 1) =   mrf(2, 2)
  mxy(1, 2) = - mrf(2, 3)
  mxy(2, 3) = - mrf(3, 1)
  mxy(2, 1) = - mrf(3, 2)
  mxy(2, 2) =   mrf(3, 3)

end subroutine mrf_mxy

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine stime1(x, n, dt, t1, t2)
  real x(n)
  do i = 1, n
     x(i) = 0
  end do
  it1 = t1 / dt + 1.5
  it2 = t2 / dt + 1.5
  it3 = (t1 + t2) / dt + 1.5
  do i = 2, it1
     x(i) = 2 * (i - 1) * dt / (t1 * t2)
  end do
  do i = it1 + 1, it2
     x(i) = 2 * (t2 - (i - 1) * dt) / (t2 - t1) / t2
  end do
  return
end subroutine stime1

subroutine resample_shift(w1, imax1, dt1, w2, imax2, dt2, shift)
  implicit none
  integer,  intent(in) :: imax1, imax2
  real,     intent(in) :: dt1, dt2
  real,     intent(in) :: shift
  real,     intent(in) :: w1(imax1)
  real,     intent(out):: w2(imax2)
  real :: trend, t1, t2
  integer :: i, i1
  w2 = 0.
  do i = 1, imax2
     t2 = dt2 * (i - 1) - shift
     i1 = int(t2/dt1) + 1
     t1 = dt1 * (i1 - 1)
     if(i1 <= 0  .or. i1 >= imax1) cycle
     if(t2 < 0.) then
        trend = (w1(1) - 0.) / dt1
        w2(i) = w1(1) + trend * (t2)
     else
        trend = (w1(i1 + 1) - w1(i1)) / dt1
        w2(i) = w1(i1) + trend * (t2 - t1)
     endif
  end do
  return
end subroutine resample_shift

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine get_xy_location(xx, yy, mn, nn, m0, n0, x, y)
  implicit none
  integer,intent(in) :: mn, nn, m0, n0
  real,   intent(in) :: xx, yy
  real               :: x(0:mn+1), y(0:nn+1)
  !- - -
  integer :: i
  !- - -
  x = 0.
  y = 0.
  x(0) = - xx * m0
  do i = 0, mn
     x(i+1) = x(i) + xx
  end do
  y(0) = - yy * n0
  do i = 0, nn
     y(i+1) = y(i) + yy
  end do
end subroutine get_xy_location

function slip_xy(d, mn, nn, m0, n0, xx, yy, xl, yl) ! bilinear - inc.
  implicit none
  integer            m1,  m2,  n1,  n2,  mn,  nn,  m0,  n0,  m,  n
  real               xx,  yy,  xl,  yl,  x1,  x2,  y1,  y2,  slip_xy
  real               d(0:mn + 2, 0:nn + 2)
  real, allocatable:: x(:) ,  y(:)
  allocate(x(0:mn + 2), y(0:nn + 2))
  do m = 0, mn + 1
     x(m) = xx * (m - m0)
  enddo
  do n = 0, nn + 2
     y(n) = yy * (n - n0)
  enddo
  do m = 1, mn + 1
     m1 = m - 1
     m2 = m
     if(x(m).gt.xl) exit
  enddo
  do n = 1, nn + 2
     n1 = n - 1
     n2 = n
     if(y(n).gt.yl) exit
  enddo
  x1 = 1. - abs(x(m1) - xl) / xx
  x2 = 1. - abs(x(m2) - xl) / xx
  y1 = 1. - abs(y(n1) - yl) / yy
  y2 = 1. - abs(y(n2) - yl) / yy
  slip_xy =  d(m1, n1) * x1 * y1 + d(m2, n1) * x2 * y1   &
       + d(m1, n2) * x1 * y2 + d(m2, n2) * x2 * y2
  return
end function slip_xy

!------------------------------------------------                       
function getdeglat2 (alat, alon)   
  !-
  !-  partial(latitude)/partial(km)
  !-
  parameter (rad = 0.017453292) 
  parameter (delta = 0.1) 
  !------ GRS80
  a =6378.137
  b = 6356.752
  e = sqrt(1. - (b/a)**2)
  rc = a*(1.-e**2)/sqrt( (1.-(e*sin( (alat+delta/2.)*rad))**2)**3 )
  dist = rc * delta * rad
  getdeglat2 = delta/dist 
  return 
end function getdeglat2
!------------------------------------------------                       
function getdeglon2 (alat, alon) 
  !-
  !-  partial(longitude)/partial(km)
  !-
  parameter (rad = 0.017453292) 
  parameter (delta = 0.1)
  !------ GRS80
  a =6378.137
  b =6356.752
  e = sqrt(1. - (b/a)**2)
  rc = (a/sqrt( 1.- (e* sin(alat*rad))**2)) * cos(alat*rad)
  dist = rc * delta * rad
  getdeglon2 = delta / dist
  return 
end function getdeglon2
!---------------------------------------------                          
subroutine getLoc1 (dx, dy, strike, dip, rlat, rlon, wlat, wlon)
  parameter (rad = 0.017453292)
  dy = dy * cos (dip * rad)
  str1 = (strike-90.) * rad
  work1 = dx * cos (str1) + dy * sin (str1)
  work2 = - dx * sin (str1) + dy * cos (str1)
  deglat = getdeglat2 (rlat, rlon)
  deglon = getdeglon2 (rlat, rlon)
  wlat = rlat + deglat * work2
  wlon = rlon + deglon * work1
  RETURN
end subroutine getLoc1
!--------------------------------------
subroutine xy2geo(olat,olon,xd,yd,str,dip,wlat,wlon)
!	IF You Put Strike, Dip, Epicenter and center of 
!	Subfaulte location, To output Center of Subfaulte
!	Longitude and Latitude.
  parameter (rad = 0.017453292)
  deglat = getdeglat2 (olat, olon) 
  deglon = getdeglon2 (olat, olon) 
  dip1=dip*rad
  str1=(str-90.)*rad
     work1 =  xd*cos(str1) + yd*sin(str1)*cos(dip1)
     work2 = -xd*sin(str1) + yd*cos(str1)*cos(dip1)
     wlon = olon + deglon * work1
     wlat = olat + deglat * work2
  if(wlon > 180. ) wlon = wlon - 360.
  return
end subroutine xy2geo

subroutine getknot(strike, dip, nl, nk, l0, k0, dl, dk, alat, alon, &
     adepth, latitude, longitude, depth)
  implicit none
  !- - - - - - - - - - - - - - - - - - - - -
  integer,intent(in)  :: nl, nk, l0, k0
  real,   intent(in)  :: strike, dip, alat, alon, adepth, dl, dk
  real,   intent(out) :: latitude(1:nl, 1:nk), longitude(1:nl, 1:nk), depth(1:nl, 1:nk)
  integer :: l, k
  real    :: wx, wy
  real    :: rad = 0.0174533
  !- - - - - - - - - - - - - - - - - - - - -
  do k = 1, nk
     do l = 1, nl
        depth(l, k) = adepth + (k0 - k) * dk * sin(dip * rad)
        wx = dl * (l - l0)
        wy = dk * (k - k0)
        call xy2geo(alat, alon, wx, wy, strike, dip, latitude(l, k), longitude(l, k))
     end do
  end do
end subroutine getknot

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine meca_rotation(m, rad)
  implicit none
  real, intent(in)    :: rad
  real, intent(inout) :: m(3,3)
  real                :: strike, dip
  real                :: rad_theta
  real                :: Rz(3,3),Ry(3,3),R(3,3)
  !
  ! Right-hand coordinate system with Euler angles rotation
  ! X -> North direction
  ! Y -> East direction
  ! Z -> Down direction
  !
  ! positive Z rotation X -> Y rotation (strike)
  ! positive X rotation Y -> Z rotation (dip)
  !
  open(11, file="meca_rotation.info", status="old", err=999)
  read(11,*) strike, dip
  close(11)
  rad_theta=strike*rad
  Rz(1,1)=cos(rad_theta); Rz(1,2)=-1.*sin(rad_theta); Rz(1,3)=0.
  Rz(2,1)=sin(rad_theta); Rz(2,2)=cos(rad_theta);     Rz(2,3)=0.
  Rz(3,1)=0.;             Rz(3,2)=0.;                 Rz(3,3)=1.
  rad_theta=dip*rad
  Ry(1,1)=cos(rad_theta);     Ry(1,2)=0.; Ry(1,3)=sin(rad_theta)
  Ry(2,1)=0.;                 Ry(2,2)=1.; Ry(2,3)=0.
  Ry(3,1)=-1.*sin(rad_theta); Ry(3,2)=0.; Ry(3,3)=cos(rad_theta)
  R=matmul(Ry,Rz) ! strike -> dip
  m=matmul(matmul(R,m),transpose(R))
  return
999 m=m
  return
end subroutine

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine faultline(strike,dip,alat,alon,adepth,xx,yy,mn,nn,m0,n0,inter)
  implicit none
  integer, intent(in) :: mn,nn,m0,n0
  real, intent(in)    :: strike,dip,alat,alon,adepth,xx,yy,inter
  integer             :: m,n
  real                :: dx,dy,bdepth,udepth,wx,wy,&
                         lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4
  real, parameter     :: rad = 3.1415927/180.
  ! bottom left
  m=1; n=1
  bdepth=adepth+( (n0-n)*yy + yy/inter )*sin(dip*rad) ! bottom
  wy=yy*(n-n0); wx=xx*(m-m0); dx=wx-xx/inter; dy=wy-yy/inter; call xy2geo(alat,alon,dx,dy,strike,dip,lat1,lon1)
  ! upper left
  m=1; n=nn
  udepth=adepth+( (n0-n)*yy - yy/inter )*sin(dip*rad) ! upper
  wy=yy*(n-n0); wx=xx*(m-m0); dx=wx-xx/inter; dy=wy+yy/inter; call xy2geo(alat,alon,dx,dy,strike,dip,lat2,lon2)
  ! upper right
  m=mn; n=nn
  wy=yy*(n-n0); wx=xx*(m-m0); dx=wx+xx/inter; dy=wy+yy/inter; call xy2geo(alat,alon,dx,dy,strike,dip,lat3,lon3)
  ! bottom right
  m=mn; n=1
  wy=yy*(n-n0); wx=xx*(m-m0); dx=wx+xx/inter; dy=wy-yy/inter; call xy2geo(alat,alon,dx,dy,strike,dip,lat4,lon4)
  open(13,file='faultline.dat',status='replace')
  write(13,'(3f9.3)')lat1,lon1,bdepth
  write(13,'(3f9.3)')lat2,lon2,udepth
  write(13,'(3f9.3)')lat3,lon3,udepth
  write(13,'(3f9.3)')lat4,lon4,bdepth
  write(13,'(3f9.3)')lat1,lon1,bdepth
  close(13)
  return
end subroutine
!-----------------------------------------------------
subroutine fslip2mij(xmo,icmn,strike,dip,rake, mij)
  implicit none
  integer, intent(in) :: icmn
  real,    intent(in) :: strike, dip, rake
  real,    intent(in) :: xmo(1:6)
  real,    intent(out):: mij(3,3)
  real                :: w1_mij(3,3),w2_mij(3,3)
  !----------
  if(icmn .eq. 1 ) then
    call sdr2mij(strike,dip,rake, w1_mij)
    mij = xmo(1) * w1_mij
  else
    call sdr2mij(strike,dip,rake-45, w1_mij)
    call sdr2mij(strike,dip,rake+45, w2_mij)
    mij = xmo(1) * w1_mij  + xmo(2) *  w2_mij
  end if
  return
end subroutine fslip2mij
!-----------------------------------------------------
subroutine sdr2mij(strike,dip,rake, mij)
  !--- 
  !   get mij  - xyz  Aki & Richards
  !---
  implicit none
  real,parameter   :: rad=0.017453292
  real,intent(in)  :: strike,dip,rake
  real,intent(out) :: mij(3,3)
  mij = 0.
  mij(1,1) = -(sin(dip*rad)*cos(rake*rad)*sin(2.*strike*rad) +     sin(2.*dip*rad)*sin(rake*rad)*(sin(strike*rad)**2))
  mij(1,2) =   sin(dip*rad)*cos(rake*rad)*cos(2.*strike*rad) + 0.5*sin(2.*dip*rad)*sin(rake*rad)*sin(2.*strike*rad)
  mij(1,3) = -(cos(dip*rad)*cos(rake*rad)*cos(strike*rad)    +     cos(2.*dip*rad)*sin(rake*rad)*sin(strike*rad))
  mij(2,2) =   sin(dip*rad)*cos(rake*rad)*sin(2.*strike*rad) -     sin(2.*dip*rad)*sin(rake*rad)*(cos(strike*rad)**2)
  mij(2,3) = -(cos(dip*rad)*cos(rake*rad)*sin(strike*rad)    -     cos(2.*dip*rad)*sin(rake*rad)*cos(strike*rad))
  mij(3,3) =   sin(2.*dip*rad)*sin(rake*rad)
  mij(2,1) = mij(1,2)
  mij(3,1) = mij(1,3)
  mij(3,2) = mij(2,3)
  return
end subroutine sdr2mij
!-----------------------------------------------------
subroutine selection_nodal (str1, dip1, str2, dip2, str0, dip0, rad ,nodal_p )
  implicit none
  real,intent(in):: str1, dip1, str2, dip2, str0, dip0, rad
  integer,intent(out):: nodal_p
  real :: n1(3), n2(3), n0(3)
  !-----
  n1(1) = -sin(dip1*rad)*sin(str1*rad)
  n1(2) =  sin(dip1*rad)*cos(str1*rad)
  n1(3) = -cos(dip1*rad)
  n2(1) = -sin(dip2*rad)*sin(str2*rad)
  n2(2) =  sin(dip2*rad)*cos(str2*rad)
  n2(3) = -cos(dip2*rad)
  n0(1) = -sin(dip0*rad)*sin(str0*rad)
  n0(2) =  sin(dip0*rad)*cos(str0*rad)
  n0(3) = -cos(dip0*rad)
  if( abs(dot_product( n1, n0)) >= abs(dot_product( n2, n0)) ) then
    nodal_p = 1
  else
    nodal_p = 2
  endif
end subroutine selection_nodal
!-----------------------------------------------------
subroutine  swap_fault (str1, dip1, rake1, str2, dip2,rake2, nodal_p )
  implicit none
  real,intent(inout):: str1, dip1, rake1, str2, dip2,rake2
  integer,intent(in):: nodal_p
  real:: str0, dip0, rake0
  if(nodal_p .eq. 2) then 
    str0 = str1
    dip0 = dip1
    rake0 = rake1
    str1 = str2
    dip1 = dip2
    rake1 = rake2
    str2 = str0
    dip2 = dip0
    rake2 = rake0
  endif
end subroutine swap_fault
