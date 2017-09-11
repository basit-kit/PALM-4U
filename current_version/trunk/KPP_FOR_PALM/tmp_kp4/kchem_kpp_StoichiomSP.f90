! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Stoichiometric Data Structures File
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
! File                 : kchem_kpp_StoichiomSP.f90
! Time                 : Thu Sep  7 11:39:19 2017
! Working directory    : /pd/home/khan-b/github/palm-4u/current_version/trunk/KPP_FOR_PALM/tmp_kp4
! Equation file        : kchem_kpp.kpp
! Output root filename : kchem_kpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE kchem_kpp_StoichiomSP

  USE kchem_kpp_Precision
  PUBLIC
  SAVE


! Row-compressed sparse data for the Jacobian of reaction products JVRP

  INTEGER, PARAMETER, DIMENSION(13) :: CROW_JVRP = (/ &
       1,  2,  3,  5,  7,  9, 10, 12, 14, 16, 18, 20, &
      21 /)

  INTEGER, PARAMETER, DIMENSION(20) :: ICOL_JVRP = (/ &
      13,  3,  6, 12,  4, 11,  9, 11,  9,  7, 12, 10, &
      12,  8, 12, 11, 13,  8, 13,  5 /)

  INTEGER, PARAMETER, DIMENSION(20) :: IROW_JVRP = (/ &
       1,  2,  3,  3,  4,  4,  5,  5,  6,  7,  7,  8, &
       8,  9,  9, 10, 10, 11, 11, 12 /)



!  Stoichiometric Matrix in Compressed Column Sparse Format


  INTEGER, PARAMETER, DIMENSION(13) :: CCOL_STOICM = (/ &
       1,  4,  6,  9, 12, 15, 19, 23, 28, 32, 35, 38, &
      41 /)

  INTEGER, PARAMETER, DIMENSION(40) :: IROW_STOICM = (/ &
       3, 12, 13,  3,  6,  6, 12, 13,  4, 10, 11,  8, &
       9, 11,  2,  7,  9, 10,  7, 11, 12, 13,  7,  9, &
      10, 12, 13,  8, 10, 12, 13,  1, 11, 13,  5,  8, &
      13,  5,  8, 13 /)

  INTEGER, PARAMETER, DIMENSION(40) :: ICOL_STOICM = (/ &
       1,  1,  1,  2,  2,  3,  3,  3,  4,  4,  4,  5, &
       5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8, &
       8,  8,  8,  9,  9,  9,  9, 10, 10, 10, 11, 11, &
      11, 12, 12, 12 /)

  REAL(kind=dp), PARAMETER, DIMENSION(40) :: STOICM = (/ &
       1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp /)


END MODULE kchem_kpp_StoichiomSP

