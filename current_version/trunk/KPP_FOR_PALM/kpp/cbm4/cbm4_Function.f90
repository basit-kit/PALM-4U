! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Function of Chemical Model File
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
! File                 : cbm4_Function.f90
! Time                 : Mon Mar  6 12:48:45 2017
! Working directory    : /pd/home/khan-b/palm/current_version/trunk/KPP_FOR_PALM/kpp/cbm
! Equation file        : cbm4.kpp
! Output root filename : cbm4
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm4_Function

  USE cbm4_Parameters
  IMPLICIT NONE

! A - Rate for each equation
  REAL(kind=dp) :: A(NREACT)

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Fun - time derivatives of variables - Agregate form
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      Vdot      - Time derivative of variable species concentrations
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Fun ( V, F, RCT, Vdot )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(NREACT)
! Vdot - Time derivative of variable species concentrations
  REAL(kind=dp) :: Vdot(NVAR)


! Computation of equation rates
  A(1) = RCT(1)*V(26)
  A(2) = RCT(2)*V(29)
  A(3) = RCT(3)*V(25)*V(31)
  A(4) = 9.3e-12*V(26)*V(29)
  A(5) = RCT(5)*V(26)*V(29)
  A(6) = RCT(6)*V(29)*V(31)
  A(7) = RCT(7)*V(25)*V(26)
  A(8) = RCT(8)*V(25)
  A(9) = RCT(9)*V(25)
  A(10) = RCT(10)*V(1)
  A(11) = 2.2e-10*V(1)*F(1)
  A(12) = RCT(12)*V(25)*V(27)
  A(13) = RCT(13)*V(25)*V(28)
  A(14) = RCT(14)*V(30)
  A(15) = RCT(15)*V(30)*V(31)
  A(16) = RCT(16)*V(26)*V(30)
  A(17) = RCT(17)*V(26)*V(30)
  A(18) = 1.3e-21*V(6)*F(1)
  A(19) = RCT(19)*V(6)
  A(20) = RCT(20)*V(31)*V(31)
  A(21) = 4.39999e-40*V(26)*V(31)*F(1)
  A(22) = RCT(22)*V(27)*V(31)
  A(23) = RCT(23)*V(9)
  A(24) = 6.6e-12*V(9)*V(27)
  A(25) = 1e-20*V(9)*V(9)
  A(26) = RCT(26)*V(26)*V(27)
  A(27) = RCT(27)*V(12)*V(27)
  A(28) = RCT(28)*V(28)*V(31)
  A(29) = RCT(29)*V(26)*V(28)
  A(30) = RCT(30)*V(10)
  A(31) = RCT(31)*V(10)*V(27)
  A(32) = RCT(32)*V(28)*V(28)
  A(33) = RCT(33)*V(28)*V(28)*F(1)
  A(34) = RCT(34)*V(2)
  A(35) = RCT(35)*V(2)*V(27)
  A(36) = 2.2e-13*V(16)*V(27)
  A(37) = 1e-11*V(21)*V(27)
  A(38) = RCT(38)*V(21)
  A(39) = RCT(39)*V(21)
  A(40) = RCT(40)*V(21)*V(29)
  A(41) = 6.3e-16*V(21)*V(30)
  A(42) = RCT(42)*V(24)*V(29)
  A(43) = RCT(43)*V(24)*V(27)
  A(44) = 2.5e-15*V(24)*V(30)
  A(45) = RCT(45)*V(24)
  A(46) = RCT(46)*V(31)*V(32)
  A(47) = RCT(47)*V(26)*V(32)
  A(48) = RCT(48)*V(3)
  A(49) = 2e-12*V(32)*V(32)
  A(50) = 6.5e-12*V(28)*V(32)
  A(51) = RCT(51)*V(27)
  A(52) = 8.1e-13*V(20)*V(27)
  A(53) = RCT(53)*V(13)
  A(54) = 1600*V(13)
  A(55) = 1.5e-11*V(13)*V(26)
  A(56) = RCT(56)*V(23)*V(29)
  A(57) = RCT(57)*V(23)*V(27)
  A(58) = RCT(58)*V(23)*V(25)
  A(59) = 7.7e-15*V(23)*V(30)
  A(60) = RCT(60)*V(17)*V(29)
  A(61) = RCT(61)*V(17)*V(27)
  A(62) = RCT(62)*V(17)*V(25)
  A(63) = RCT(63)*V(5)*V(27)
  A(64) = 8.1e-12*V(11)*V(31)
  A(65) = 4.2*V(11)
  A(66) = 4.1e-11*V(14)*V(27)
  A(67) = 2.2e-11*V(14)*V(30)
  A(68) = 1.4e-11*V(4)*V(26)
  A(69) = RCT(69)*V(7)*V(27)
  A(70) = 3e-11*V(19)*V(27)
  A(71) = RCT(71)*V(19)
  A(72) = RCT(72)*V(19)*V(25)
  A(73) = 1.7e-11*V(15)*V(27)
  A(74) = RCT(74)*V(15)
  A(75) = 1.8e-11*V(22)*V(29)
  A(76) = 9.6e-11*V(22)*V(27)
  A(77) = 1.2e-17*V(22)*V(25)
  A(78) = 3.2e-13*V(22)*V(30)
  A(79) = 8.1e-12*V(18)*V(31)
  A(80) = RCT(80)*V(18)*V(18)
  A(81) = 6.8e-13*V(8)*V(31)

! Aggregate function
  Vdot(1) = A(9)-A(10)-A(11)
  Vdot(2) = A(32)+A(33)-A(34)-A(35)
  Vdot(3) = A(47)-A(48)
  Vdot(4) = 0.4*A(66)+A(67)-A(68)
  Vdot(5) = -A(63)
  Vdot(6) = A(17)-A(18)-A(19)
  Vdot(7) = -A(69)
  Vdot(8) = 0.13*A(52)+0.04*A(53)+0.02*A(56)+0.09*A(59)+0.13*A(76)+A(78)-A(81)
  Vdot(9) = 2*A(21)+A(22)-A(23)-A(24)-2*A(25)
  Vdot(10) = A(29)-A(30)-A(31)
  Vdot(11) = 0.56*A(63)-A(64)-A(65)+0.3*A(69)
  Vdot(12) = 2*A(18)+A(26)-A(27)+A(41)+A(44)+A(67)
  Vdot(13) = 0.76*A(52)-0.98*A(53)-A(54)-A(55)
  Vdot(14) = 0.36*A(63)+A(65)-A(66)-A(67)+0.2*A(69)
  Vdot(15) = 0.8*A(69)+0.2*A(72)-A(73)-A(74)+0.4*A(76)+0.2*A(77)
  Vdot(16) = -A(36)+A(37)+A(38)+A(39)+A(40)+A(41)+A(45)+0.3*A(56)+0.33*A(58)+A(60)+0.42*A(62)+2*A(70)+A(71)+0.69*A(72)&
               &+A(74)+0.5*A(75)+0.06*A(77)
  Vdot(17) = -A(60)-A(61)-A(62)+0.45*A(75)+A(76)+0.55*A(77)
  Vdot(18) = A(45)+A(46)+2*A(49)+0.79*A(50)+A(51)+0.87*A(52)+0.96*A(53)+0.28*A(56)+A(57)+0.22*A(58)+0.91*A(59)+0.7*A(60)&
               &+A(61)+0.08*A(63)+0.6*A(66)+0.5*A(69)+A(70)+0.03*A(72)+A(73)+0.5*A(75)+A(76)-A(79)-2*A(80)
  Vdot(19) = 0.9*A(64)+0.3*A(66)-A(70)-A(71)-A(72)
  Vdot(20) = -1.11*A(52)-2.1*A(53)+0.22*A(56)-A(57)-A(58)-A(59)+1.1*A(69)+0.9*A(75)+0.1*A(77)
  Vdot(21) = -A(37)-A(38)-A(39)-A(40)-A(41)+A(45)+A(46)+2*A(49)+0.79*A(50)+A(51)+0.2*A(56)+A(57)+0.74*A(58)+A(59)+A(60)&
               &+1.56*A(61)+A(62)+A(70)+0.7*A(72)+A(76)+A(77)
  Vdot(22) = -A(75)-A(76)-A(77)-A(78)
  Vdot(23) = -A(56)-A(57)-A(58)-A(59)+0.55*A(75)
  Vdot(24) = -A(42)-A(43)-A(44)-A(45)+0.11*A(52)+1.1*A(53)+0.63*A(56)+A(57)+0.5*A(58)+A(59)+0.22*A(61)+0.03*A(72)+0.8&
               &*A(75)+0.2*A(76)+0.4*A(77)
  Vdot(25) = A(2)-A(3)-A(7)-A(8)-A(9)-A(12)-A(13)-A(58)-A(62)-A(72)-A(77)
  Vdot(26) = -A(1)+A(3)-A(4)-A(5)+A(6)-A(7)+0.89*A(14)+2*A(15)-A(17)+A(19)+2*A(20)-A(21)+A(24)+A(25)-A(26)+A(28)-A(29)&
               &+A(30)+A(31)+A(46)-A(47)+A(48)-A(55)+A(59)+0.9*A(64)-A(68)+A(79)
  Vdot(27) = 2*A(11)-A(12)+A(13)-A(22)+A(23)-A(24)-A(26)-A(27)+A(28)-A(31)+2*A(34)-A(35)-A(36)-A(37)+A(40)+A(42)-A(43)&
               &+0.79*A(50)-A(51)-A(52)+0.2*A(56)-A(57)+0.1*A(58)+0.3*A(60)-A(61)-A(63)-A(66)-A(69)-A(70)+0.08*A(72)-A(73)&
               &-A(76)+0.1*A(77)
  Vdot(28) = A(12)-A(13)-A(28)-A(29)+A(30)-2*A(32)-2*A(33)+A(35)+A(36)+A(37)+2*A(38)+A(40)+A(41)+2*A(45)+A(46)+2*A(49)&
               &-0.21*A(50)+A(51)+0.11*A(52)+0.94*A(53)+A(54)+0.38*A(56)+A(57)+0.44*A(58)+1.7*A(60)+A(61)+0.12*A(62)+0.44&
               &*A(63)+0.9*A(64)+A(65)+0.6*A(66)+0.7*A(69)+2*A(70)+A(71)+0.76*A(72)+A(74)+0.6*A(75)+0.67*A(76)+0.44*A(77)
  Vdot(29) = A(1)-A(2)-A(4)-A(5)-A(6)+A(8)+A(10)+0.89*A(14)-A(40)-A(42)-A(56)-A(60)-A(75)
  Vdot(30) = A(5)+A(7)-A(14)-A(15)-A(16)-A(17)+A(19)+A(27)-A(41)-A(44)-A(59)-A(67)-A(78)
  Vdot(31) = A(1)-A(3)+A(4)-A(6)+0.11*A(14)-A(15)+A(16)-2*A(20)-A(21)-A(22)+A(23)+A(25)-A(28)-A(46)-A(64)-A(79)-A(81)
  Vdot(32) = A(42)+A(43)+A(44)-A(46)-A(47)+A(48)-2*A(49)-A(50)+A(70)+A(71)+0.62*A(72)+A(73)+A(74)+0.2*A(76)
      
END SUBROUTINE Fun

! End of Fun function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE cbm4_Function

