! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Jacobian Data Structures File
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
! File                 : smog_JacobianSP.f90
! Time                 : Mon Mar 13 15:58:02 2017
! Working directory    : /pd/home/khan-b/palm/current_version/trunk/KPP_FOR_PALM/kpp/kchem_kpp
! Equation file        : smog.kpp
! Output root filename : smog
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE smog_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(59) :: LU_IROW = (/ &
       1,  1,  1,  2,  2,  3,  3,  4,  4,  4,  5,  5, &
       5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7, &
       8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9, 10, &
      10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, &
      11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12 /)

  INTEGER, PARAMETER, DIMENSION(59) :: LU_ICOL = (/ &
       1, 11, 12,  2, 11,  3, 12,  4,  7, 11,  2,  5, &
      10, 11,  6,  8,  9, 10,  4,  7,  8, 10, 11, 12, &
       8,  9, 10, 12,  3,  7,  8,  9, 10, 11, 12,  5, &
       6,  7,  8,  9, 10, 11, 12,  4,  5,  6,  7,  8, &
       9, 10, 11, 12,  3,  6,  8,  9, 10, 11, 12 /)

  INTEGER, PARAMETER, DIMENSION(13) :: LU_CROW = (/ &
       1,  4,  6,  8, 11, 15, 19, 25, 29, 36, 44, 53, &
      60 /)

  INTEGER, PARAMETER, DIMENSION(13) :: LU_DIAG = (/ &
       1,  4,  6,  8, 12, 15, 20, 25, 32, 41, 51, 59, &
      60 /)


END MODULE smog_JacobianSP

