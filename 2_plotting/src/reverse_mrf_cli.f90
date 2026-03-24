!***********************************************************************
!
! reverse_mrf_cli.f90
!
! This version is modified to accept 6 moment tensor components
! as command-line arguments and print the results in a parseable format.
!
! [How to Compile and Link]
! This program MUST be linked with the pre-compiled object files
! of its dependencies.
!
! 1. Create object files (if not already done):
!    gfortran -c sub.gen.f90
!    gfortran -c FPSPACK.FOR
!
! 2. Compile this command-line version:
!    gfortran -c reverse_mrf_cli.f90
!
! 3. Link all object files together:
!    gfortran reverse_mrf_cli.o sub.gen.o FPSPACK.o -o reverse_mrf_v2
!
! [How to Run]
! ./reverse_mrf_v2 <Mrr> <Mss> <Mee> <Mrs> <Mre> <Mse>
!
!***********************************************************************
program reverse_calculation_cli
    implicit none
  
    !----- Input Moment Tensor (MRF) -----
    real :: mrf(1:3, 1:3)
    real :: mrf_rr, mrf_ss, mrf_ee, mrf_rs, mrf_re, mrf_se
    character(len=20) :: arg
  
    !----- Calculated Variables -----
    real :: xmo(1:6)
    real :: strike, dip, rake
    real :: a_strike, a_dip, a_rake
    real :: sliprate
    real :: trendp, trendt, trendb, plungp, plungt, plungb
    real :: NDC
    integer :: i
  
    ! --- Get MRF values from command-line arguments ---
    if (command_argument_count() < 6) then
      write(*,*) "ERROR: 6 moment tensor components are required."
      write(*,*) "USAGE: ./reverse_mrf_v2 Mrr Mss Mee Mrs Mre Mse"
      stop 1
    end if
  
    call get_command_argument(1, arg)
    read(arg, *) mrf_rr
    call get_command_argument(2, arg)
    read(arg, *) mrf_ss
    call get_command_argument(3, arg)
    read(arg, *) mrf_ee
    call get_command_argument(4, arg)
    read(arg, *) mrf_rs
    call get_command_argument(5, arg)
    read(arg, *) mrf_re
    call get_command_argument(6, arg)
    read(arg, *) mrf_se
  
    ! Construct the mrf matrix
    mrf(1,1) = mrf_rr
    mrf(2,2) = mrf_ss
    mrf(3,3) = mrf_ee
    mrf(1,2) = mrf_rs; mrf(2,1) = mrf(1,2) ! Mrt (Mzx)
    mrf(1,3) = mrf_re; mrf(3,1) = mrf(1,3) ! Mrp (Mzy)
    mrf(2,3) = mrf_se; mrf(3,2) = mrf(2,3) ! Mtp (Mxy)
  
    !----- Start Calculation -----
    ! Step 1: Inverse calculation (mrf -> xmo)
    call inverse_mrf_to_xmo(mrf, xmo)
  
    ! Step 2: Forward calculation from the inverted xmo
    call calculate_parameters_from_xmo(xmo, strike, dip, rake, &
      & a_strike, a_dip, a_rake, sliprate, &
      & trendp, trendt, trendb, plungp, plungt, plungb, NDC)
  
    !----- Print Results in Key:Value format for easy parsing -----
    write(*,'(A,F20.10)') "SLIPRATE:", sliprate
    write(*,'(A,F20.10)') "STRIKE1:", strike
    write(*,'(A,F20.10)') "DIP1:", dip
    write(*,'(A,F20.10)') "RAKE1:", rake
    write(*,'(A,F20.10)') "STRIKE2:", a_strike
    write(*,'(A,F20.10)') "DIP2:", a_dip
    write(*,'(A,F20.10)') "RAKE2:", a_rake
    write(*,'(A,F20.10)') "TRENDP:", trendp
    write(*,'(A,F20.10)') "PLUNGP:", plungp
    write(*,'(A,F20.10)') "TRENDT:", trendt
    write(*,'(A,F20.10)') "PLUNGT:", plungt
    write(*,'(A,F20.10)') "TRENDB:", trendb
    write(*,'(A,F20.10)') "PLUNGB:", plungb
    write(*,'(A,F20.10)') "NDC:", NDC
  
  contains
  
  !-----------------------------------------------------------------------
  ! Inverse calculation from MRF to XMO
  !-----------------------------------------------------------------------
  subroutine inverse_mrf_to_xmo(mrf, xmo)
    implicit none
    real, intent(in)  :: mrf(1:3, 1:3)
    real, intent(out) :: xmo(1:6)
    real :: mij(1:3, 1:3)
  
    ! This call requires 'mrf_mxy' from sub.gen.o
    call mrf_mxy(mrf, mij)
  
    ! This call requires 'mtrx_inverse' defined below
    call mtrx_inverse(mij, xmo)
  
  end subroutine inverse_mrf_to_xmo
  
  !-----------------------------------------------------------------------
  ! Main routine to calculate all parameters from XMO
  !-----------------------------------------------------------------------
  subroutine calculate_parameters_from_xmo(xmo, strike, dip, rake, &
      & a_strike, a_dip, a_rake, sliprate, &
      & trendp, trendt, trendb, plungp, plungt, plungb, NDC)
    implicit none
    real,    intent(in)  :: xmo(1:6)
    real,    intent(out) :: strike, dip, rake, a_strike, a_dip, a_rake, sliprate
    real,    intent(out) :: trendp, trendt, trendb, plungp, plungt, plungb, NDC
    
    real    :: mij(1:3, 1:3), mrf(1:3, 1:3), Mxy(1:3, 1:3)
    real    :: w_1(1:3), eigenval(1:3), EV(1:3, 1:3)
    real    :: dm, am0, am1, e, am0b, slipa, slipb, eta
    integer :: ierr, dum
  
    ! The following calls require subroutines from sub.gen.o and FPSPACK.o
    call d_cp(xmo, strike, dip, rake, sliprate, dm, w_1)
    
    call conj(strike, dip, rake, a_strike, a_dip, a_rake)
  
    call mtrx(xmo, mij)
    
    call AR2PLP(mij, am0,am1,e,am0b,strike,dip,rake,slipa,&
          & a_strike,a_dip,a_rake,slipb,&
          & trendp,plungp,trendt,plungt,trendb,plungb,eta,ierr)
    
    call mxy_mrf(mij, mrf)
    call mrf_mxy(mrf, Mxy)
    call EIG1(Mxy, 3, 3, eigenval, EV, dum)
    if (max(abs(eigenval(1)), abs(eigenval(3))) > 1.0e-20) then
        NDC = - eigenval(2) / max(abs(eigenval(1)), abs(eigenval(3)))
        NDC = 100 * NDC / 0.5
    else
        NDC = 0.0
    end if
  
  end subroutine calculate_parameters_from_xmo
  
  !-----------------------------------------------------------------------
  ! Newly created subroutine to perform inverse of mtrx
  !-----------------------------------------------------------------------
  subroutine mtrx_inverse(m, v)
    real, intent(in)  :: m(3,3)
    real, intent(out) :: v(6)
    
    v(1) = m(1,2)
    v(3) = m(2,3)
    v(4) = m(1,3)
    
    v(6) = (m(1,1) + m(2,2) + m(3,3)) / 3.0
    v(2) = v(6) - m(2,2)
    v(5) = m(3,3) - v(6)
  
  end subroutine mtrx_inverse
  
  end program reverse_calculation_cli
  