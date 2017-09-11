! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Hessian File
! 
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : kchem_kpp_Hessian.f90
! Time                 : Thu Sep  7 11:39:19 2017
! Working directory    : /pd/home/khan-b/github/palm-4u/current_version/trunk/KPP_FOR_PALM/tmp_kp4
! Equation file        : kchem_kpp.kpp
! Output root filename : kchem_kpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE kchem_kpp_Hessian

  USE kchem_kpp_Parameters
  USE kchem_kpp_HessianSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Hessian - function for Hessian (Jac derivative w.r.t. variables)
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      HESS      - Hessian of Var (i.e. the 3-tensor d Jac / d Var)
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Hessian ( V, F, RCT, HESS )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(NREACT)
! HESS - Hessian of Var (i.e. the 3-tensor d Jac / d Var)
  REAL(kind=dp) :: HESS(NHESS)

! --------------------------------------------------------
! Note: HESS is represented in coordinate sparse format:
!       HESS(m) = d^2 f_i / dv_j dv_k = d Jac_{i,j} / dv_k
!       where i = IHESS_I(m), j = IHESS_J(m), k = IHESS_K(m).
! --------------------------------------------------------
! Note: d^2 f_i / dv_j dv_k = d^2 f_i / dv_k dv_j,
!       therefore only the terms d^2 f_i / dv_j dv_k
!       with j <= k are computed and stored in HESS.
! --------------------------------------------------------

! Local variables
! D2A - Second derivatives of equation rates
  REAL(kind=dp) :: D2A(8)

! Computation of the second derivatives of equation rates
! D2A(1) = d^2 A(3) / dV(6)dV(12)
  D2A(1) = RCT(3)
! D2A(2) = d^2 A(4) / dV(4)dV(11)
  D2A(2) = RCT(4)
! D2A(3) = d^2 A(5) / dV(9)dV(11)
  D2A(3) = RCT(5)
! D2A(4) = d^2 A(7) / dV(7)dV(12)
  D2A(4) = RCT(7)
! D2A(5) = d^2 A(8) / dV(10)dV(12)
  D2A(5) = RCT(8)
! D2A(6) = d^2 A(9) / dV(8)dV(12)
  D2A(6) = RCT(9)
! D2A(7) = d^2 A(10) / dV(11)dV(13)
  D2A(7) = RCT(10)
! D2A(8) = d^2 A(11) / dV(8)dV(13)
  D2A(8) = RCT(11)

! Computation of the Jacobian derivative
! HESS(1) = d^2 Vdot(1)/{dV(11)dV(13)} = d^2 Vdot(1)/{dV(13)dV(11)}
  HESS(1) = D2A(7)
! HESS(2) = d^2 Vdot(4)/{dV(4)dV(11)} = d^2 Vdot(4)/{dV(11)dV(4)}
  HESS(2) = -D2A(2)
! HESS(3) = d^2 Vdot(5)/{dV(8)dV(13)} = d^2 Vdot(5)/{dV(13)dV(8)}
  HESS(3) = D2A(8)
! HESS(4) = d^2 Vdot(6)/{dV(6)dV(12)} = d^2 Vdot(6)/{dV(12)dV(6)}
  HESS(4) = -D2A(1)
! HESS(5) = d^2 Vdot(7)/{dV(7)dV(12)} = d^2 Vdot(7)/{dV(12)dV(7)}
  HESS(5) = -D2A(4)
! HESS(6) = d^2 Vdot(7)/{dV(10)dV(12)} = d^2 Vdot(7)/{dV(12)dV(10)}
  HESS(6) = D2A(5)
! HESS(7) = d^2 Vdot(8)/{dV(8)dV(12)} = d^2 Vdot(8)/{dV(12)dV(8)}
  HESS(7) = -D2A(6)
! HESS(8) = d^2 Vdot(8)/{dV(8)dV(13)} = d^2 Vdot(8)/{dV(13)dV(8)}
  HESS(8) = -D2A(8)
! HESS(9) = d^2 Vdot(8)/{dV(9)dV(11)} = d^2 Vdot(8)/{dV(11)dV(9)}
  HESS(9) = D2A(3)
! HESS(10) = d^2 Vdot(9)/{dV(9)dV(11)} = d^2 Vdot(9)/{dV(11)dV(9)}
  HESS(10) = -D2A(3)
! HESS(11) = d^2 Vdot(9)/{dV(10)dV(12)} = d^2 Vdot(9)/{dV(12)dV(10)}
  HESS(11) = D2A(5)
! HESS(12) = d^2 Vdot(10)/{dV(4)dV(11)} = d^2 Vdot(10)/{dV(11)dV(4)}
  HESS(12) = D2A(2)
! HESS(13) = d^2 Vdot(10)/{dV(8)dV(12)} = d^2 Vdot(10)/{dV(12)dV(8)}
  HESS(13) = D2A(6)
! HESS(14) = d^2 Vdot(10)/{dV(10)dV(12)} = d^2 Vdot(10)/{dV(12)dV(10)}
  HESS(14) = -D2A(5)
! HESS(15) = d^2 Vdot(11)/{dV(4)dV(11)} = d^2 Vdot(11)/{dV(11)dV(4)}
  HESS(15) = -D2A(2)
! HESS(16) = d^2 Vdot(11)/{dV(7)dV(12)} = d^2 Vdot(11)/{dV(12)dV(7)}
  HESS(16) = D2A(4)
! HESS(17) = d^2 Vdot(11)/{dV(9)dV(11)} = d^2 Vdot(11)/{dV(11)dV(9)}
  HESS(17) = -D2A(3)
! HESS(18) = d^2 Vdot(11)/{dV(11)dV(13)} = d^2 Vdot(11)/{dV(13)dV(11)}
  HESS(18) = -D2A(7)
! HESS(19) = d^2 Vdot(12)/{dV(6)dV(12)} = d^2 Vdot(12)/{dV(12)dV(6)}
  HESS(19) = -D2A(1)
! HESS(20) = d^2 Vdot(12)/{dV(7)dV(12)} = d^2 Vdot(12)/{dV(12)dV(7)}
  HESS(20) = -D2A(4)
! HESS(21) = d^2 Vdot(12)/{dV(8)dV(12)} = d^2 Vdot(12)/{dV(12)dV(8)}
  HESS(21) = -D2A(6)
! HESS(22) = d^2 Vdot(12)/{dV(10)dV(12)} = d^2 Vdot(12)/{dV(12)dV(10)}
  HESS(22) = -D2A(5)
! HESS(23) = d^2 Vdot(13)/{dV(6)dV(12)} = d^2 Vdot(13)/{dV(12)dV(6)}
  HESS(23) = D2A(1)
! HESS(24) = d^2 Vdot(13)/{dV(7)dV(12)} = d^2 Vdot(13)/{dV(12)dV(7)}
  HESS(24) = D2A(4)
! HESS(25) = d^2 Vdot(13)/{dV(8)dV(12)} = d^2 Vdot(13)/{dV(12)dV(8)}
  HESS(25) = D2A(6)
! HESS(26) = d^2 Vdot(13)/{dV(8)dV(13)} = d^2 Vdot(13)/{dV(13)dV(8)}
  HESS(26) = -D2A(8)
! HESS(27) = d^2 Vdot(13)/{dV(10)dV(12)} = d^2 Vdot(13)/{dV(12)dV(10)}
  HESS(27) = D2A(5)
! HESS(28) = d^2 Vdot(13)/{dV(11)dV(13)} = d^2 Vdot(13)/{dV(13)dV(11)}
  HESS(28) = -D2A(7)
      
END SUBROUTINE Hessian

! End of Hessian function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! HessTR_Vec - Hessian transposed times user vectors
!   Arguments :
!      HESS      - Hessian of Var (i.e. the 3-tensor d Jac / d Var)
!      U1        - User vector
!      U2        - User vector
!      HTU       - Transposed Hessian times user vectors: (Hess x U2)^T * U1 = [d (Jac^T*U1)/d Var] * U2
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE HessTR_Vec ( HESS, U1, U2, HTU )

! HESS - Hessian of Var (i.e. the 3-tensor d Jac / d Var)
  REAL(kind=dp) :: HESS(NHESS)
! U1 - User vector
  REAL(kind=dp) :: U1(NVAR)
! U2 - User vector
  REAL(kind=dp) :: U2(NVAR)
! HTU - Transposed Hessian times user vectors: (Hess x U2)^T * U1 = [d (Jac^T*U1)/d Var] * U2
  REAL(kind=dp) :: HTU(NVAR)

! Compute the vector HTU =(Hess x U2)^T * U1 = d (Jac^T*U1)/d Var * U2
  HTU(1) = 0
  HTU(2) = 0
  HTU(3) = 0
  HTU(4) = HESS(2)*(U1(4)*U2(11))+HESS(12)*(U1(10)*U2(11))+HESS(15)*(U1(11)*U2(11))
  HTU(5) = 0
  HTU(6) = HESS(4)*(U1(6)*U2(12))+HESS(19)*(U1(12)*U2(12))+HESS(23)*(U1(13)*U2(12))
  HTU(7) = HESS(5)*(U1(7)*U2(12))+HESS(16)*(U1(11)*U2(12))+HESS(20)*(U1(12)*U2(12))+HESS(24)*(U1(13)*U2(12))
  HTU(8) = HESS(3)*(U1(5)*U2(13))+HESS(7)*(U1(8)*U2(12))+HESS(8)*(U1(8)*U2(13))+HESS(13)*(U1(10)*U2(12))+HESS(21)&
             &*(U1(12)*U2(12))+HESS(25)*(U1(13)*U2(12))+HESS(26)*(U1(13)*U2(13))
  HTU(9) = HESS(9)*(U1(8)*U2(11))+HESS(10)*(U1(9)*U2(11))+HESS(17)*(U1(11)*U2(11))
  HTU(10) = HESS(6)*(U1(7)*U2(12))+HESS(11)*(U1(9)*U2(12))+HESS(14)*(U1(10)*U2(12))+HESS(22)*(U1(12)*U2(12))+HESS(27)&
              &*(U1(13)*U2(12))
  HTU(11) = HESS(1)*(U1(1)*U2(13))+HESS(2)*(U1(4)*U2(4))+HESS(9)*(U1(8)*U2(9))+HESS(10)*(U1(9)*U2(9))+HESS(12)*(U1(10)&
              &*U2(4))+HESS(15)*(U1(11)*U2(4))+HESS(17)*(U1(11)*U2(9))+HESS(18)*(U1(11)*U2(13))+HESS(28)*(U1(13)*U2(13))
  HTU(12) = HESS(4)*(U1(6)*U2(6))+HESS(5)*(U1(7)*U2(7))+HESS(6)*(U1(7)*U2(10))+HESS(7)*(U1(8)*U2(8))+HESS(11)*(U1(9)&
              &*U2(10))+HESS(13)*(U1(10)*U2(8))+HESS(14)*(U1(10)*U2(10))+HESS(16)*(U1(11)*U2(7))+HESS(19)*(U1(12)*U2(6))&
              &+HESS(20)*(U1(12)*U2(7))+HESS(21)*(U1(12)*U2(8))+HESS(22)*(U1(12)*U2(10))+HESS(23)*(U1(13)*U2(6))+HESS(24)&
              &*(U1(13)*U2(7))+HESS(25)*(U1(13)*U2(8))+HESS(27)*(U1(13)*U2(10))
  HTU(13) = HESS(1)*(U1(1)*U2(11))+HESS(3)*(U1(5)*U2(8))+HESS(8)*(U1(8)*U2(8))+HESS(18)*(U1(11)*U2(11))+HESS(26)*(U1(13)&
              &*U2(8))+HESS(28)*(U1(13)*U2(11))
      
END SUBROUTINE HessTR_Vec

! End of HessTR_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Hess_Vec - Hessian times user vectors
!   Arguments :
!      HESS      - Hessian of Var (i.e. the 3-tensor d Jac / d Var)
!      U1        - User vector
!      U2        - User vector
!      HU        - Hessian times user vectors: (Hess x U2) * U1 = [d (Jac*U1)/d Var] * U2
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Hess_Vec ( HESS, U1, U2, HU )

! HESS - Hessian of Var (i.e. the 3-tensor d Jac / d Var)
  REAL(kind=dp) :: HESS(NHESS)
! U1 - User vector
  REAL(kind=dp) :: U1(NVAR)
! U2 - User vector
  REAL(kind=dp) :: U2(NVAR)
! HU - Hessian times user vectors: (Hess x U2) * U1 = [d (Jac*U1)/d Var] * U2
  REAL(kind=dp) :: HU(NVAR)

! Compute the vector HU =(Hess x U2) * U1 = d (Jac*U1)/d Var * U2
  HU(1) = HESS(1)*(U1(11)*U2(13))+HESS(1)*(U1(13)*U2(11))
  HU(2) = 0
  HU(3) = 0
  HU(4) = HESS(2)*(U1(4)*U2(11))+HESS(2)*(U1(11)*U2(4))
  HU(5) = HESS(3)*(U1(8)*U2(13))+HESS(3)*(U1(13)*U2(8))
  HU(6) = HESS(4)*(U1(6)*U2(12))+HESS(4)*(U1(12)*U2(6))
  HU(7) = HESS(5)*(U1(7)*U2(12))+HESS(5)*(U1(12)*U2(7))+HESS(6)*(U1(10)*U2(12))+HESS(6)*(U1(12)*U2(10))
  HU(8) = HESS(7)*(U1(8)*U2(12))+HESS(7)*(U1(12)*U2(8))+HESS(8)*(U1(8)*U2(13))+HESS(8)*(U1(13)*U2(8))+HESS(9)*(U1(9)&
            &*U2(11))+HESS(9)*(U1(11)*U2(9))
  HU(9) = HESS(10)*(U1(9)*U2(11))+HESS(10)*(U1(11)*U2(9))+HESS(11)*(U1(10)*U2(12))+HESS(11)*(U1(12)*U2(10))
  HU(10) = HESS(12)*(U1(4)*U2(11))+HESS(12)*(U1(11)*U2(4))+HESS(13)*(U1(8)*U2(12))+HESS(13)*(U1(12)*U2(8))+HESS(14)&
             &*(U1(10)*U2(12))+HESS(14)*(U1(12)*U2(10))
  HU(11) = HESS(15)*(U1(4)*U2(11))+HESS(15)*(U1(11)*U2(4))+HESS(16)*(U1(7)*U2(12))+HESS(16)*(U1(12)*U2(7))+HESS(17)&
             &*(U1(9)*U2(11))+HESS(17)*(U1(11)*U2(9))+HESS(18)*(U1(11)*U2(13))+HESS(18)*(U1(13)*U2(11))
  HU(12) = HESS(19)*(U1(6)*U2(12))+HESS(19)*(U1(12)*U2(6))+HESS(20)*(U1(7)*U2(12))+HESS(20)*(U1(12)*U2(7))+HESS(21)&
             &*(U1(8)*U2(12))+HESS(21)*(U1(12)*U2(8))+HESS(22)*(U1(10)*U2(12))+HESS(22)*(U1(12)*U2(10))
  HU(13) = HESS(23)*(U1(6)*U2(12))+HESS(23)*(U1(12)*U2(6))+HESS(24)*(U1(7)*U2(12))+HESS(24)*(U1(12)*U2(7))+HESS(25)&
             &*(U1(8)*U2(12))+HESS(25)*(U1(12)*U2(8))+HESS(26)*(U1(8)*U2(13))+HESS(26)*(U1(13)*U2(8))+HESS(27)*(U1(10)&
             &*U2(12))+HESS(27)*(U1(12)*U2(10))+HESS(28)*(U1(11)*U2(13))+HESS(28)*(U1(13)*U2(11))
      
END SUBROUTINE Hess_Vec

! End of Hess_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE kchem_kpp_Hessian
