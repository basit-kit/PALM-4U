! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Utility Data Module File
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
! File                 : kchem_kpp_Monitor.f90
! Time                 : Mon Mar  6 14:06:05 2017
! Working directory    : /pd/home/khan-b/palm/current_version/trunk/KPP_FOR_PALM/kp4/bin/tmp_kp4
! Equation file        : kchem_kpp.kpp
! Output root filename : kchem_kpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE kchem_kpp_Monitor


  CHARACTER(LEN=15), PARAMETER, DIMENSION(7) :: SPC_NAMES = (/ &
     'O1D            ','O              ','O3             ', &
     'NO             ','NO2            ','M              ', &
     'O2             ' /)

  INTEGER, DIMENSION(1) :: LOOKAT
  INTEGER, DIMENSION(1) :: MONITOR
  CHARACTER(LEN=15), DIMENSION(1) :: SMASS
  CHARACTER(LEN=100), PARAMETER, DIMENSION(10) :: EQN_NAMES = (/ &
     '      O2 --> 2 O                                                                                    ', &
     '  O + O2 --> O3                                                                                     ', &
     '      O3 --> O + O2                                                                                 ', &
     '  O + O3 --> 2 O2                                                                                   ', &
     '      O3 --> O1D + O2                                                                               ', &
     ' O1D + M --> O + M                                                                                  ', &
     'O1D + O3 --> 2 O2                                                                                   ', &
     ' O3 + NO --> NO2 + O2                                                                               ', &
     ' O + NO2 --> NO + O2                                                                                ', &
     '     NO2 --> O + NO                                                                                 ' /)

! INLINED global variables

! End INLINED global variables


END MODULE kchem_kpp_Monitor
