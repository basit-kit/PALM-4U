! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The Stoichiometric Chemical Model File
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
! File                 : cbm4_Stoichiom.f90
! Time                 : Mon Mar  6 12:48:45 2017
! Working directory    : /pd/home/khan-b/palm/current_version/trunk/KPP_FOR_PALM/kpp/cbm
! Equation file        : cbm4.kpp
! Output root filename : cbm4
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm4_Stoichiom

  USE cbm4_Parameters
  USE cbm4_StoichiomSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! ReactantProd - Reactant Products in each equation
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      ARP       - Reactant product in each equation
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE ReactantProd ( V, F, ARP )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! ARP - Reactant product in each equation
  REAL(kind=dp) :: ARP(NREACT)


! Reactant Products in each equation are useful in the
!     stoichiometric formulation of mass action law
  ARP(1) = V(26)
  ARP(2) = V(29)
  ARP(3) = V(25)*V(31)
  ARP(4) = V(26)*V(29)
  ARP(5) = V(26)*V(29)
  ARP(6) = V(29)*V(31)
  ARP(7) = V(25)*V(26)
  ARP(8) = V(25)
  ARP(9) = V(25)
  ARP(10) = V(1)
  ARP(11) = V(1)*F(1)
  ARP(12) = V(25)*V(27)
  ARP(13) = V(25)*V(28)
  ARP(14) = V(30)
  ARP(15) = V(30)*V(31)
  ARP(16) = V(26)*V(30)
  ARP(17) = V(26)*V(30)
  ARP(18) = V(6)*F(1)
  ARP(19) = V(6)
  ARP(20) = V(31)*V(31)
  ARP(21) = V(26)*V(31)*F(1)
  ARP(22) = V(27)*V(31)
  ARP(23) = V(9)
  ARP(24) = V(9)*V(27)
  ARP(25) = V(9)*V(9)
  ARP(26) = V(26)*V(27)
  ARP(27) = V(12)*V(27)
  ARP(28) = V(28)*V(31)
  ARP(29) = V(26)*V(28)
  ARP(30) = V(10)
  ARP(31) = V(10)*V(27)
  ARP(32) = V(28)*V(28)
  ARP(33) = V(28)*V(28)*F(1)
  ARP(34) = V(2)
  ARP(35) = V(2)*V(27)
  ARP(36) = V(16)*V(27)
  ARP(37) = V(21)*V(27)
  ARP(38) = V(21)
  ARP(39) = V(21)
  ARP(40) = V(21)*V(29)
  ARP(41) = V(21)*V(30)
  ARP(42) = V(24)*V(29)
  ARP(43) = V(24)*V(27)
  ARP(44) = V(24)*V(30)
  ARP(45) = V(24)
  ARP(46) = V(31)*V(32)
  ARP(47) = V(26)*V(32)
  ARP(48) = V(3)
  ARP(49) = V(32)*V(32)
  ARP(50) = V(28)*V(32)
  ARP(51) = V(27)
  ARP(52) = V(20)*V(27)
  ARP(53) = V(13)
  ARP(54) = V(13)
  ARP(55) = V(13)*V(26)
  ARP(56) = V(23)*V(29)
  ARP(57) = V(23)*V(27)
  ARP(58) = V(23)*V(25)
  ARP(59) = V(23)*V(30)
  ARP(60) = V(17)*V(29)
  ARP(61) = V(17)*V(27)
  ARP(62) = V(17)*V(25)
  ARP(63) = V(5)*V(27)
  ARP(64) = V(11)*V(31)
  ARP(65) = V(11)
  ARP(66) = V(14)*V(27)
  ARP(67) = V(14)*V(30)
  ARP(68) = V(4)*V(26)
  ARP(69) = V(7)*V(27)
  ARP(70) = V(19)*V(27)
  ARP(71) = V(19)
  ARP(72) = V(19)*V(25)
  ARP(73) = V(15)*V(27)
  ARP(74) = V(15)
  ARP(75) = V(22)*V(29)
  ARP(76) = V(22)*V(27)
  ARP(77) = V(22)*V(25)
  ARP(78) = V(22)*V(30)
  ARP(79) = V(18)*V(31)
  ARP(80) = V(18)*V(18)
  ARP(81) = V(8)*V(31)
      
END SUBROUTINE ReactantProd

! End of ReactantProd function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! JacReactantProd - Jacobian of Reactant Products vector
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      JVRP      - d ARP(1:NREACT)/d VAR (1:NVAR)
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE JacReactantProd ( V, F, JVRP )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! JVRP - d ARP(1:NREACT)/d VAR (1:NVAR)
  REAL(kind=dp) :: JVRP(NJVRP)


! Reactant Products in each equation are useful in the
!    stoichiometric formulation of mass action law
! Below we compute the Jacobian of the Reactant Products vector
!    w.r.t. variable species: d ARP(1:NREACT) / d Var(1:NVAR)

! JVRP(1) = dARP(1)/dV(26)
  JVRP(1) = 1
! JVRP(2) = dARP(2)/dV(29)
  JVRP(2) = 1
! JVRP(3) = dARP(3)/dV(25)
  JVRP(3) = V(31)
! JVRP(4) = dARP(3)/dV(31)
  JVRP(4) = V(25)
! JVRP(5) = dARP(4)/dV(26)
  JVRP(5) = V(29)
! JVRP(6) = dARP(4)/dV(29)
  JVRP(6) = V(26)
! JVRP(7) = dARP(5)/dV(26)
  JVRP(7) = V(29)
! JVRP(8) = dARP(5)/dV(29)
  JVRP(8) = V(26)
! JVRP(9) = dARP(6)/dV(29)
  JVRP(9) = V(31)
! JVRP(10) = dARP(6)/dV(31)
  JVRP(10) = V(29)
! JVRP(11) = dARP(7)/dV(25)
  JVRP(11) = V(26)
! JVRP(12) = dARP(7)/dV(26)
  JVRP(12) = V(25)
! JVRP(13) = dARP(8)/dV(25)
  JVRP(13) = 1
! JVRP(14) = dARP(9)/dV(25)
  JVRP(14) = 1
! JVRP(15) = dARP(10)/dV(1)
  JVRP(15) = 1
! JVRP(16) = dARP(11)/dV(1)
  JVRP(16) = F(1)
! JVRP(17) = dARP(12)/dV(25)
  JVRP(17) = V(27)
! JVRP(18) = dARP(12)/dV(27)
  JVRP(18) = V(25)
! JVRP(19) = dARP(13)/dV(25)
  JVRP(19) = V(28)
! JVRP(20) = dARP(13)/dV(28)
  JVRP(20) = V(25)
! JVRP(21) = dARP(14)/dV(30)
  JVRP(21) = 1
! JVRP(22) = dARP(15)/dV(30)
  JVRP(22) = V(31)
! JVRP(23) = dARP(15)/dV(31)
  JVRP(23) = V(30)
! JVRP(24) = dARP(16)/dV(26)
  JVRP(24) = V(30)
! JVRP(25) = dARP(16)/dV(30)
  JVRP(25) = V(26)
! JVRP(26) = dARP(17)/dV(26)
  JVRP(26) = V(30)
! JVRP(27) = dARP(17)/dV(30)
  JVRP(27) = V(26)
! JVRP(28) = dARP(18)/dV(6)
  JVRP(28) = F(1)
! JVRP(29) = dARP(19)/dV(6)
  JVRP(29) = 1
! JVRP(30) = dARP(20)/dV(31)
  JVRP(30) = 2*V(31)
! JVRP(31) = dARP(21)/dV(26)
  JVRP(31) = V(31)*F(1)
! JVRP(32) = dARP(21)/dV(31)
  JVRP(32) = V(26)*F(1)
! JVRP(33) = dARP(22)/dV(27)
  JVRP(33) = V(31)
! JVRP(34) = dARP(22)/dV(31)
  JVRP(34) = V(27)
! JVRP(35) = dARP(23)/dV(9)
  JVRP(35) = 1
! JVRP(36) = dARP(24)/dV(9)
  JVRP(36) = V(27)
! JVRP(37) = dARP(24)/dV(27)
  JVRP(37) = V(9)
! JVRP(38) = dARP(25)/dV(9)
  JVRP(38) = 2*V(9)
! JVRP(39) = dARP(26)/dV(26)
  JVRP(39) = V(27)
! JVRP(40) = dARP(26)/dV(27)
  JVRP(40) = V(26)
! JVRP(41) = dARP(27)/dV(12)
  JVRP(41) = V(27)
! JVRP(42) = dARP(27)/dV(27)
  JVRP(42) = V(12)
! JVRP(43) = dARP(28)/dV(28)
  JVRP(43) = V(31)
! JVRP(44) = dARP(28)/dV(31)
  JVRP(44) = V(28)
! JVRP(45) = dARP(29)/dV(26)
  JVRP(45) = V(28)
! JVRP(46) = dARP(29)/dV(28)
  JVRP(46) = V(26)
! JVRP(47) = dARP(30)/dV(10)
  JVRP(47) = 1
! JVRP(48) = dARP(31)/dV(10)
  JVRP(48) = V(27)
! JVRP(49) = dARP(31)/dV(27)
  JVRP(49) = V(10)
! JVRP(50) = dARP(32)/dV(28)
  JVRP(50) = 2*V(28)
! JVRP(51) = dARP(33)/dV(28)
  JVRP(51) = 2*V(28)*F(1)
! JVRP(52) = dARP(34)/dV(2)
  JVRP(52) = 1
! JVRP(53) = dARP(35)/dV(2)
  JVRP(53) = V(27)
! JVRP(54) = dARP(35)/dV(27)
  JVRP(54) = V(2)
! JVRP(55) = dARP(36)/dV(16)
  JVRP(55) = V(27)
! JVRP(56) = dARP(36)/dV(27)
  JVRP(56) = V(16)
! JVRP(57) = dARP(37)/dV(21)
  JVRP(57) = V(27)
! JVRP(58) = dARP(37)/dV(27)
  JVRP(58) = V(21)
! JVRP(59) = dARP(38)/dV(21)
  JVRP(59) = 1
! JVRP(60) = dARP(39)/dV(21)
  JVRP(60) = 1
! JVRP(61) = dARP(40)/dV(21)
  JVRP(61) = V(29)
! JVRP(62) = dARP(40)/dV(29)
  JVRP(62) = V(21)
! JVRP(63) = dARP(41)/dV(21)
  JVRP(63) = V(30)
! JVRP(64) = dARP(41)/dV(30)
  JVRP(64) = V(21)
! JVRP(65) = dARP(42)/dV(24)
  JVRP(65) = V(29)
! JVRP(66) = dARP(42)/dV(29)
  JVRP(66) = V(24)
! JVRP(67) = dARP(43)/dV(24)
  JVRP(67) = V(27)
! JVRP(68) = dARP(43)/dV(27)
  JVRP(68) = V(24)
! JVRP(69) = dARP(44)/dV(24)
  JVRP(69) = V(30)
! JVRP(70) = dARP(44)/dV(30)
  JVRP(70) = V(24)
! JVRP(71) = dARP(45)/dV(24)
  JVRP(71) = 1
! JVRP(72) = dARP(46)/dV(31)
  JVRP(72) = V(32)
! JVRP(73) = dARP(46)/dV(32)
  JVRP(73) = V(31)
! JVRP(74) = dARP(47)/dV(26)
  JVRP(74) = V(32)
! JVRP(75) = dARP(47)/dV(32)
  JVRP(75) = V(26)
! JVRP(76) = dARP(48)/dV(3)
  JVRP(76) = 1
! JVRP(77) = dARP(49)/dV(32)
  JVRP(77) = 2*V(32)
! JVRP(78) = dARP(50)/dV(28)
  JVRP(78) = V(32)
! JVRP(79) = dARP(50)/dV(32)
  JVRP(79) = V(28)
! JVRP(80) = dARP(51)/dV(27)
  JVRP(80) = 1
! JVRP(81) = dARP(52)/dV(20)
  JVRP(81) = V(27)
! JVRP(82) = dARP(52)/dV(27)
  JVRP(82) = V(20)
! JVRP(83) = dARP(53)/dV(13)
  JVRP(83) = 1
! JVRP(84) = dARP(54)/dV(13)
  JVRP(84) = 1
! JVRP(85) = dARP(55)/dV(13)
  JVRP(85) = V(26)
! JVRP(86) = dARP(55)/dV(26)
  JVRP(86) = V(13)
! JVRP(87) = dARP(56)/dV(23)
  JVRP(87) = V(29)
! JVRP(88) = dARP(56)/dV(29)
  JVRP(88) = V(23)
! JVRP(89) = dARP(57)/dV(23)
  JVRP(89) = V(27)
! JVRP(90) = dARP(57)/dV(27)
  JVRP(90) = V(23)
! JVRP(91) = dARP(58)/dV(23)
  JVRP(91) = V(25)
! JVRP(92) = dARP(58)/dV(25)
  JVRP(92) = V(23)
! JVRP(93) = dARP(59)/dV(23)
  JVRP(93) = V(30)
! JVRP(94) = dARP(59)/dV(30)
  JVRP(94) = V(23)
! JVRP(95) = dARP(60)/dV(17)
  JVRP(95) = V(29)
! JVRP(96) = dARP(60)/dV(29)
  JVRP(96) = V(17)
! JVRP(97) = dARP(61)/dV(17)
  JVRP(97) = V(27)
! JVRP(98) = dARP(61)/dV(27)
  JVRP(98) = V(17)
! JVRP(99) = dARP(62)/dV(17)
  JVRP(99) = V(25)
! JVRP(100) = dARP(62)/dV(25)
  JVRP(100) = V(17)
! JVRP(101) = dARP(63)/dV(5)
  JVRP(101) = V(27)
! JVRP(102) = dARP(63)/dV(27)
  JVRP(102) = V(5)
! JVRP(103) = dARP(64)/dV(11)
  JVRP(103) = V(31)
! JVRP(104) = dARP(64)/dV(31)
  JVRP(104) = V(11)
! JVRP(105) = dARP(65)/dV(11)
  JVRP(105) = 1
! JVRP(106) = dARP(66)/dV(14)
  JVRP(106) = V(27)
! JVRP(107) = dARP(66)/dV(27)
  JVRP(107) = V(14)
! JVRP(108) = dARP(67)/dV(14)
  JVRP(108) = V(30)
! JVRP(109) = dARP(67)/dV(30)
  JVRP(109) = V(14)
! JVRP(110) = dARP(68)/dV(4)
  JVRP(110) = V(26)
! JVRP(111) = dARP(68)/dV(26)
  JVRP(111) = V(4)
! JVRP(112) = dARP(69)/dV(7)
  JVRP(112) = V(27)
! JVRP(113) = dARP(69)/dV(27)
  JVRP(113) = V(7)
! JVRP(114) = dARP(70)/dV(19)
  JVRP(114) = V(27)
! JVRP(115) = dARP(70)/dV(27)
  JVRP(115) = V(19)
! JVRP(116) = dARP(71)/dV(19)
  JVRP(116) = 1
! JVRP(117) = dARP(72)/dV(19)
  JVRP(117) = V(25)
! JVRP(118) = dARP(72)/dV(25)
  JVRP(118) = V(19)
! JVRP(119) = dARP(73)/dV(15)
  JVRP(119) = V(27)
! JVRP(120) = dARP(73)/dV(27)
  JVRP(120) = V(15)
! JVRP(121) = dARP(74)/dV(15)
  JVRP(121) = 1
! JVRP(122) = dARP(75)/dV(22)
  JVRP(122) = V(29)
! JVRP(123) = dARP(75)/dV(29)
  JVRP(123) = V(22)
! JVRP(124) = dARP(76)/dV(22)
  JVRP(124) = V(27)
! JVRP(125) = dARP(76)/dV(27)
  JVRP(125) = V(22)
! JVRP(126) = dARP(77)/dV(22)
  JVRP(126) = V(25)
! JVRP(127) = dARP(77)/dV(25)
  JVRP(127) = V(22)
! JVRP(128) = dARP(78)/dV(22)
  JVRP(128) = V(30)
! JVRP(129) = dARP(78)/dV(30)
  JVRP(129) = V(22)
! JVRP(130) = dARP(79)/dV(18)
  JVRP(130) = V(31)
! JVRP(131) = dARP(79)/dV(31)
  JVRP(131) = V(18)
! JVRP(132) = dARP(80)/dV(18)
  JVRP(132) = 2*V(18)
! JVRP(133) = dARP(81)/dV(8)
  JVRP(133) = V(31)
! JVRP(134) = dARP(81)/dV(31)
  JVRP(134) = V(8)
      
END SUBROUTINE JacReactantProd

! End of JacReactantProd function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! Begin Derivative w.r.t. Rate Coefficients

! ------------------------------------------------------------------------------
! Subroutine for the derivative of Fun with respect to rate coefficients
! -----------------------------------------------------------------------------

      SUBROUTINE  dFun_dRcoeff( V, F, NCOEFF, JCOEFF, DFDR )
       
      USE cbm4_Parameters
      USE cbm4_StoichiomSP
      IMPLICIT NONE 

! V - Concentrations of variable/radical/fixed species            
      REAL(kind=dp) V(NVAR), F(NFIX)
! NCOEFF - the number of rate coefficients with respect to which we differentiate
      INTEGER NCOEFF       
! JCOEFF - a vector of integers containing the indices of reactions (rate
!          coefficients) with respect to which we differentiate
      INTEGER JCOEFF(NCOEFF)       
! DFDR  - a matrix containg derivative values; specifically, 
!         column j contains d Fun(1:NVAR) / d RCT( JCOEFF(j) )
!         for each 1 <= j <= NCOEFF
!         This matrix is stored in a column-wise linearized format
      REAL(kind=dp) DFDR(NVAR*NCOEFF)

! Local vector with reactant products
      REAL(kind=dp) A_RPROD(NREACT)
      REAL(kind=dp) aj
      INTEGER i,j,k
      
! Compute the reactant products of all reactions     
      CALL ReactantProd ( V, F, A_RPROD )

! Compute the derivatives by multiplying column JCOEFF(j) of the stoichiometric matrix with A_RPROD       
      DO j=1,NCOEFF
!                  Initialize the j-th column of derivative matrix to zero       
         DO i=1,NVAR
           DFDR(i+NVAR*(j-1)) = 0.0_dp 
         END DO
!                  Column JCOEFF(j) in the stoichiometric matrix times the
!                  reactant product  of the JCOEFF(j)-th reaction      
!                  give the j-th column of the derivative matrix   
         aj = A_RPROD(JCOEFF(j))
         DO k=CCOL_STOICM(JCOEFF(j)),CCOL_STOICM(JCOEFF(j)+1)-1
           DFDR(IROW_STOICM(k)+NVAR*(j-1)) = STOICM(k)*aj
         END DO
      END DO
      
      END SUBROUTINE  dFun_dRcoeff

! End Derivative w.r.t. Rate Coefficients


! Begin Jacobian Derivative w.r.t. Rate Coefficients

! ------------------------------------------------------------------------------
! Subroutine for the derivative of Jac with respect to rate coefficients
! Times a user vector
! -----------------------------------------------------------------------------

      SUBROUTINE  dJac_dRcoeff( V, F, U, NCOEFF, JCOEFF, DJDR )
       
      USE cbm4_Parameters
      USE cbm4_StoichiomSP
      IMPLICIT NONE 

! V - Concentrations of variable/fixed species            
      REAL(kind=dp) V(NVAR), F(NFIX)
! U - User-supplied Vector           
      REAL(kind=dp) U(NVAR)
! NCOEFF - the number of rate coefficients with respect to which we differentiate
      INTEGER NCOEFF       
! JCOEFF - a vector of integers containing the indices of reactions (rate
!          coefficients) with respect to which we differentiate
      INTEGER JCOEFF(NCOEFF)       
! DFDR  - a matrix containg derivative values; specifically, 
!         column j contains d Jac(1:NVAR) / d RCT( JCOEFF(j) ) * U
!                     for each 1 <= j <= NCOEFF
!         This matrix is stored in a column-wise linearized format
      REAL(kind=dp) DJDR(NVAR*NCOEFF)

! Local vector for Jacobian of reactant products
      REAL(kind=dp) JV_RPROD(NJVRP)
      REAL(kind=dp) aj
      INTEGER i,j,k
      
! Compute the Jacobian of all reactant products   
      CALL JacReactantProd( V, F, JV_RPROD )

! Compute the derivatives by multiplying column JCOEFF(j) of the stoichiometric matrix with A_PROD       
      DO j=1,NCOEFF
!                  Initialize the j-th column of derivative matrix to zero       
         DO i=1,NVAR
           DJDR(i+NVAR*(j-1)) = 0.0_dp
         END DO
!                  Column JCOEFF(j) in the stoichiometric matrix times the
!                  ( Gradient of reactant product of the JCOEFF(j)-th reaction X user vector )    
!                  give the j-th column of the derivative matrix   
!
!          Row JCOEFF(j) of JV_RPROD times the user vector
         aj = 0.0_dp
         DO k=CROW_JVRP(JCOEFF(j)),CROW_JVRP(JCOEFF(j)+1)-1
             aj = aj + JV_RPROD(k)*U(ICOL_JVRP(k))
         END DO
!          Column JCOEFF(j) of Stoichiom. matrix times aj         
         DO k=CCOL_STOICM(JCOEFF(j)),CCOL_STOICM(JCOEFF(j)+1)-1
           DJDR(IROW_STOICM(k)+NVAR*(j-1)) = STOICM(k)*aj
         END DO
      END DO
      
      END SUBROUTINE  dJac_dRcoeff

! End Jacobian Derivative w.r.t. Rate Coefficients


END MODULE cbm4_Stoichiom

