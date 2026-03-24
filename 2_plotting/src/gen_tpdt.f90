program fgentpdt
  implicit none
  real, allocatable :: disp(:,:), disp1(:,:,:), disp2(:,:,:,:),&
                       tr(:,:), slip(:,:)
  real :: xmo(1:6),mij(1:3,1:3),mrf(1:3,1:3)
  real :: moment, Mw, rigid, lat, lon, depth, vr, &
          strike, dip, slip0, slipi, rslip, &
          xx, yy, rtime, ylength, &
          varicance, ABIC
  integer :: nsurface, mn, nn, m0, n0, jtn, icmn
  integer :: n, m, jt, icm
  real    :: nexpo, maxcomp, m1, m2, m3, m4, m5, m6
  real    :: strike1, dip1, rake1, strike2, dip2, rake2, sm, NDC
  real, parameter :: rad=3.1415927/180.
  !---
  real am,ama,am0,am1,e,am0b,strika,dipa,rakea,slipa,strikb,dipb,&
       rakeb,slipb,trendp,plungp,trendt,plungt,trendb,plungb,eta
  integer ierr
  real :: eigenval(1:3),mxy(1:3,1:3)
  !---
  open(24,file='tpdt.dat',status='replace')
  !---
  read(40,*); read(40,*)moment,Mw,rigid ,lat,lon,depth,vr,nsurface
  read(40,*); read(40,*)strike,dip,slip0,slipi,rslip
  read(40,*); read(40,*)xx,yy,mn,nn,m0,n0,rtime,jtn,icmn,ylength
  read(40,*); read(40,*)varicance, ABIC
  read(40,*)
  allocate (disp(0:mn+2,0:nn+2),disp1(0:mn+2,0:nn+2,1:icmn), &
       disp2(1:mn,1:nn,1:jtn,1:icmn),tr(1:mn,1:nn),slip(0:mn+1,0:nn+1))
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
  xmo = 0. ; mij = 0. ; mrf = 0.
  do icm = 1, icmn
    do n = 1, nn
      do m = 1, mn
        xmo(icm) = xmo(icm) + disp1(m,n,icm)! * xx * yy * 1.e6 ! potency
      end do
    end do
  end do
  call MTRX(xmo,mij)
  if( icmn .ne. 5 )call fslip2mij(xmo,icmn,strike,dip,slipi, mij)
  call d_cp(mij,strika,dipa,rakea,sm,am0,mrf(1,1:3)) ! mrf & am0 are dummy returns
  !---
  mxy=mij
  call AR2PLP(mxy, am0,am1,e,am0b,strika,dipa,rakea,slipa,&
              strikb,dipb,rakeb,slipb,trendp,plungp,trendt,plungt,trendb,&
              plungb,eta,ierr)
  mxy=mij
  call EIG1(mxy, 3, 3, eigenval, mrf, am0) ! mrf & am0 are dummy returns
  NDC = - eigenval(2) / max(abs(eigenval(1)), abs(eigenval(3)))
  NDC = 100 * NDC / 0.5 ! convert to parcentage
  !---
  call Mxy_Mrf(mij,mrf)
  nexpo = int(log10(maxval(abs(mrf))))
  maxcomp = maxval(abs(mrf(1:3, 1:3)))
  mrf(1:3, 1:3) = mrf(1:3, 1:3) / maxcomp
  !---
  write(24,'(6f16.9)') mrf(1,1),mrf(2,2),mrf(3,3),mrf(1,2),mrf(1,3),mrf(2,3)
  write(24,'(6f10.3)') strika,dipa,rakea,strikb,dipb,rakeb
  write(24,'(6f10.3)') sm, NDC, 0., 0., 0., 0.
  close(24)
  !---
end program
