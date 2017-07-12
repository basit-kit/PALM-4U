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
! File                 : cbm4_StoichiomSP.f90
! Time                 : Mon Mar  6 12:48:45 2017
! Working directory    : /pd/home/khan-b/palm/current_version/trunk/KPP_FOR_PALM/kpp/cbm
! Equation file        : cbm4.kpp
! Output root filename : cbm4
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm4_StoichiomSP

  USE cbm4_Precision
  PUBLIC
  SAVE


! Row-compressed sparse data for the Jacobian of reaction products JVRP

  INTEGER, PARAMETER, DIMENSION(82) :: CROW_JVRP = (/ &
       1,  2,  3,  5,  7,  9, 11, 13, 14, 15, 16, 17, &
      19, 21, 22, 24, 26, 28, 29, 30, 31, 33, 35, 36, &
      38, 39, 41, 43, 45, 47, 48, 50, 51, 52, 53, 55, &
      57, 59, 60, 61, 63, 65, 67, 69, 71, 72, 74, 76, &
      77, 78, 80, 81, 83, 84, 85, 87, 89, 91, 93, 95, &
      97, 99,101,103,105,106,108,110,112,114,116,117, &
     119,121,122,124,126,128,130,132,133,135 /)

  INTEGER, PARAMETER, DIMENSION(134) :: ICOL_JVRP = (/ &
      26, 29, 25, 31, 26, 29, 26, 29, 29, 31, 25, 26, &
      25, 25,  1,  1, 25, 27, 25, 28, 30, 30, 31, 26, &
      30, 26, 30,  6,  6, 31, 26, 31, 27, 31,  9,  9, &
      27,  9, 26, 27, 12, 27, 28, 31, 26, 28, 10, 10, &
      27, 28, 28,  2,  2, 27, 16, 27, 21, 27, 21, 21, &
      21, 29, 21, 30, 24, 29, 24, 27, 24, 30, 24, 31, &
      32, 26, 32,  3, 32, 28, 32, 27, 20, 27, 13, 13, &
      13, 26, 23, 29, 23, 27, 23, 25, 23, 30, 17, 29, &
      17, 27, 17, 25,  5, 27, 11, 31, 11, 14, 27, 14, &
      30,  4, 26,  7, 27, 19, 27, 19, 19, 25, 15, 27, &
      15, 22, 29, 22, 27, 22, 25, 22, 30, 18, 31, 18, &
       8, 31 /)

  INTEGER, PARAMETER, DIMENSION(134) :: IROW_JVRP = (/ &
       1,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7, &
       8,  9, 10, 11, 12, 12, 13, 13, 14, 15, 15, 16, &
      16, 17, 17, 18, 19, 20, 21, 21, 22, 22, 23, 24, &
      24, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 31, &
      31, 32, 33, 34, 35, 35, 36, 36, 37, 37, 38, 39, &
      40, 40, 41, 41, 42, 42, 43, 43, 44, 44, 45, 46, &
      46, 47, 47, 48, 49, 50, 50, 51, 52, 52, 53, 54, &
      55, 55, 56, 56, 57, 57, 58, 58, 59, 59, 60, 60, &
      61, 61, 62, 62, 63, 63, 64, 64, 65, 66, 66, 67, &
      67, 68, 68, 69, 69, 70, 70, 71, 72, 72, 73, 73, &
      74, 75, 75, 76, 76, 77, 77, 78, 78, 79, 79, 80, &
      81, 81 /)



!  Stoichiometric Matrix in Compressed Column Sparse Format


  INTEGER, PARAMETER, DIMENSION(82) :: CCOL_STOICM = (/ &
       1,  4,  6,  9, 12, 15, 18, 21, 23, 25, 27, 29, &
      32, 35, 39, 42, 44, 47, 49, 52, 54, 57, 60, 63, &
      66, 69, 72, 75, 79, 82, 85, 88, 90, 92, 94, 97, &
     100,104,107,109,114,119,123,126,130,135,141,144, &
     147,151,156,160,167,173,175,177,187,194,203,211, &
     218,224,229,235,240,243,249,253,255,263,270,274, &
     284,288,292,301,311,321,324,327,328,330 /)

  INTEGER, PARAMETER, DIMENSION(329) :: IROW_STOICM = (/ &
      26, 29, 31, 25, 29, 25, 26, 31, 26, 29, 31, 26, &
      29, 30, 26, 29, 31, 25, 26, 30, 25, 29,  1, 25, &
       1, 29,  1, 27, 25, 27, 28, 25, 27, 28, 26, 29, &
      30, 31, 26, 30, 31, 30, 31,  6, 26, 30,  6, 12, &
       6, 26, 30, 26, 31,  9, 26, 31,  9, 27, 31,  9, &
      27, 31,  9, 26, 27,  9, 26, 31, 12, 26, 27, 12, &
      27, 30, 26, 27, 28, 31, 10, 26, 28, 10, 26, 28, &
      10, 26, 27,  2, 28,  2, 28,  2, 27,  2, 27, 28, &
      16, 27, 28, 16, 21, 27, 28, 16, 21, 28, 16, 21, &
      16, 21, 27, 28, 29, 12, 16, 21, 28, 30, 24, 27, &
      29, 32, 24, 27, 32, 12, 24, 30, 32, 16, 18, 21, &
      24, 28, 18, 21, 26, 28, 31, 32,  3, 26, 32,  3, &
      26, 32, 18, 21, 28, 32, 18, 21, 27, 28, 32, 18, &
      21, 27, 28,  8, 13, 18, 20, 24, 27, 28,  8, 13, &
      18, 20, 24, 28, 13, 28, 13, 26,  8, 16, 18, 20, &
      21, 23, 24, 27, 28, 29, 18, 20, 21, 23, 24, 27, &
      28, 16, 18, 20, 21, 23, 24, 25, 27, 28,  8, 18, &
      20, 21, 23, 24, 26, 30, 16, 17, 18, 21, 27, 28, &
      29, 17, 18, 21, 24, 27, 28, 16, 17, 21, 25, 28, &
       5, 11, 14, 18, 27, 28, 11, 19, 26, 28, 31, 11, &
      14, 28,  4, 14, 18, 19, 27, 28,  4, 12, 14, 30, &
       4, 26,  7, 11, 14, 15, 18, 20, 27, 28, 16, 18, &
      19, 21, 27, 28, 32, 16, 19, 28, 32, 15, 16, 18, &
      19, 21, 24, 25, 27, 28, 32, 15, 18, 27, 32, 15, &
      16, 28, 32, 16, 17, 18, 20, 22, 23, 24, 28, 29, &
       8, 15, 17, 18, 21, 22, 24, 27, 28, 32, 15, 16, &
      17, 20, 21, 22, 24, 25, 27, 28,  8, 22, 30, 18, &
      26, 31, 18,  8, 31 /)

  INTEGER, PARAMETER, DIMENSION(329) :: ICOL_STOICM = (/ &
       1,  1,  1,  2,  2,  3,  3,  3,  4,  4,  4,  5, &
       5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  9,  9, &
      10, 10, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, &
      14, 14, 15, 15, 15, 16, 16, 17, 17, 17, 18, 18, &
      19, 19, 19, 20, 20, 21, 21, 21, 22, 22, 22, 23, &
      23, 23, 24, 24, 24, 25, 25, 25, 26, 26, 26, 27, &
      27, 27, 28, 28, 28, 28, 29, 29, 29, 30, 30, 30, &
      31, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 35, &
      36, 36, 36, 37, 37, 37, 37, 38, 38, 38, 39, 39, &
      40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 42, 42, &
      42, 42, 43, 43, 43, 44, 44, 44, 44, 45, 45, 45, &
      45, 45, 46, 46, 46, 46, 46, 46, 47, 47, 47, 48, &
      48, 48, 49, 49, 49, 49, 50, 50, 50, 50, 50, 51, &
      51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 53, 53, &
      53, 53, 53, 53, 54, 54, 55, 55, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, &
      57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 59, 59, &
      59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, &
      60, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, &
      63, 63, 63, 63, 63, 63, 64, 64, 64, 64, 64, 65, &
      65, 65, 66, 66, 66, 66, 66, 66, 67, 67, 67, 67, &
      68, 68, 69, 69, 69, 69, 69, 69, 69, 69, 70, 70, &
      70, 70, 70, 70, 70, 71, 71, 71, 71, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 73, 73, 73, 73, 74, &
      74, 74, 74, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 78, 78, 78, 79, &
      79, 79, 80, 81, 81 /)

  REAL(kind=dp), PARAMETER, DIMENSION(150) :: STOICM_0 = (/ &
       -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  2.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  8.900000e-01_dp, &
       8.900000e-01_dp,  -1.000000e+00_dp,  1.100000e-01_dp,  2.000000e+00_dp,  -1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  2.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       1.000000e+00_dp,  2.000000e+00_dp,  -2.000000e+00_dp,  2.000000e+00_dp,  -1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       -2.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -2.000000e+00_dp,  1.000000e+00_dp, &
       -2.000000e+00_dp,  -1.000000e+00_dp,  2.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       2.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  2.000000e+00_dp,  1.000000e+00_dp, &
       1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       1.000000e+00_dp,  2.000000e+00_dp,  2.000000e+00_dp,  2.000000e+00_dp,  -2.000000e+00_dp /)
  REAL(kind=dp), PARAMETER, DIMENSION(150) :: STOICM_1 = (/ &
       7.900000e-01_dp,  7.900000e-01_dp,  7.900000e-01_dp,  -2.100000e-01_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.300000e-01_dp, &
       7.600000e-01_dp,  8.700000e-01_dp,  -1.110000e+00_dp,  1.100000e-01_dp,  -1.000000e+00_dp, &
       1.100000e-01_dp,  4.000000e-02_dp,  -9.800000e-01_dp,  9.600000e-01_dp,  -2.100000e+00_dp, &
       1.100000e+00_dp,  9.400000e-01_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       -1.000000e+00_dp,  2.000000e-02_dp,  3.000000e-01_dp,  2.800000e-01_dp,  2.200000e-01_dp, &
       2.000000e-01_dp,  -1.000000e+00_dp,  6.300000e-01_dp,  2.000000e-01_dp,  3.800000e-01_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  3.300000e-01_dp,  2.200000e-01_dp, &
       -1.000000e+00_dp,  7.400000e-01_dp,  -1.000000e+00_dp,  5.000000e-01_dp,  -1.000000e+00_dp, &
       1.000000e-01_dp,  4.400000e-01_dp,  9.000000e-02_dp,  9.100000e-01_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  7.000000e-01_dp,  1.000000e+00_dp,  3.000000e-01_dp, &
       1.700000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.560000e+00_dp, &
       2.200000e-01_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  4.200000e-01_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  1.200000e-01_dp,  -1.000000e+00_dp,  5.600000e-01_dp, &
       3.600000e-01_dp,  8.000000e-02_dp,  -1.000000e+00_dp,  4.400000e-01_dp,  -1.000000e+00_dp, &
       9.000000e-01_dp,  9.000000e-01_dp,  9.000000e-01_dp,  -1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  1.000000e+00_dp,  4.000000e-01_dp,  -1.000000e+00_dp,  6.000000e-01_dp, &
       3.000000e-01_dp,  -1.000000e+00_dp,  6.000000e-01_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp, &
       3.000000e-01_dp,  2.000000e-01_dp,  8.000000e-01_dp,  5.000000e-01_dp,  1.100000e+00_dp, &
       -1.000000e+00_dp,  7.000000e-01_dp,  2.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  2.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp,  2.000000e-01_dp,  6.900000e-01_dp, &
       3.000000e-02_dp,  -1.000000e+00_dp,  7.000000e-01_dp,  3.000000e-02_dp,  -1.000000e+00_dp, &
       8.000000e-02_dp,  7.600000e-01_dp,  6.200000e-01_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       1.000000e+00_dp,  5.000000e-01_dp,  4.500000e-01_dp,  5.000000e-01_dp,  9.000000e-01_dp, &
       -1.000000e+00_dp,  5.500000e-01_dp,  8.000000e-01_dp,  6.000000e-01_dp,  -1.000000e+00_dp /)
  REAL(kind=dp), PARAMETER, DIMENSION(29) :: STOICM_2 = (/ &
       1.300000e-01_dp,  4.000000e-01_dp,  1.000000e+00_dp,  1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  2.000000e-01_dp,  -1.000000e+00_dp,  6.700000e-01_dp,  2.000000e-01_dp, &
       2.000000e-01_dp,  6.000000e-02_dp,  5.500000e-01_dp,  1.000000e-01_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  4.000000e-01_dp,  -1.000000e+00_dp,  1.000000e-01_dp,  4.400000e-01_dp, &
       1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp,  1.000000e+00_dp, &
       -1.000000e+00_dp,  -2.000000e+00_dp,  -1.000000e+00_dp,  -1.000000e+00_dp /)
  REAL(kind=dp), PARAMETER, DIMENSION(329) :: STOICM = (/&
    STOICM_0, STOICM_1, STOICM_2 /)


END MODULE cbm4_StoichiomSP

