 MODULE pmc_mpi_wrapper

!------------------------------------------------------------------------------!
! This file is part of PALM.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2017 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: pmc_mpi_wrapper_mod.f90 2425 2017-09-11 14:21:39Z basit $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1901 2016-05-04 15:39:38Z raasch
! Code clean up. The words server/client changed to parent/child. 
!
! 1900 2016-05-04 15:27:53Z raasch
! re-formatted to match PALM style
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
!
! 1808 2016-04-05 19:44:00Z raasch
! MPI module used by default on all machines
!
! 1779 2016-03-03 08:01:28Z raasch
! kind=dp replaced by wp
!
! 1764 2016-02-28 12:45:19Z raasch
! cpp-statement added (nesting can only be used in parallel mode),
! kind-parameters adjusted to PALM-kinds
!
! 1762 2016-02-25 12:31:13Z hellstea
! Initial revision by K. Ketelsen
!
! Description:
! ------------
!
! MPI Wrapper of Palm Model Coupler
!-------------------------------------------------------------------------------!

#if defined( __parallel )
    USE, INTRINSIC ::  ISO_C_BINDING

#if defined( __mpifh )
    INCLUDE "mpif.h"
#else
    USE MPI
#endif

    USE kinds
    USE pmc_handle_communicator,                                                &
        ONLY: m_model_comm, m_model_rank, m_to_parent_comm, m_to_child_comm

    IMPLICIT NONE

    PRIVATE
    SAVE

    INTERFACE pmc_send_to_parent
       MODULE PROCEDURE pmc_send_to_parent_integer
       MODULE PROCEDURE pmc_send_to_parent_integer_2
       MODULE PROCEDURE pmc_send_to_parent_real_r1
       MODULE PROCEDURE pmc_send_to_parent_real_r2
       MODULE PROCEDURE pmc_send_to_parent_real_r3
    END INTERFACE pmc_send_to_parent

    INTERFACE pmc_recv_from_parent
       MODULE PROCEDURE pmc_recv_from_parent_integer
       MODULE PROCEDURE pmc_recv_from_parent_real_r1
       MODULE PROCEDURE pmc_recv_from_parent_real_r2
       MODULE PROCEDURE pmc_recv_from_parent_real_r3
    END INTERFACE pmc_recv_from_parent

    INTERFACE pmc_send_to_child
       MODULE PROCEDURE pmc_send_to_child_integer
       MODULE PROCEDURE pmc_send_to_child_real_r1
       MODULE PROCEDURE pmc_send_to_child_real_r2
       MODULE PROCEDURE pmc_send_to_child_real_r3
    END INTERFACE pmc_send_to_child

    INTERFACE pmc_recv_from_child
       MODULE PROCEDURE pmc_recv_from_child_integer
       MODULE PROCEDURE pmc_recv_from_child_integer_2
       MODULE PROCEDURE pmc_recv_from_child_real_r1
       MODULE PROCEDURE pmc_recv_from_child_real_r2
       MODULE PROCEDURE pmc_recv_from_child_real_r3
    END INTERFACE pmc_recv_from_child

    INTERFACE pmc_bcast
       MODULE PROCEDURE pmc_bcast_integer
       MODULE PROCEDURE pmc_bcast_character
    END INTERFACE pmc_bcast

    INTERFACE pmc_inter_bcast
       MODULE PROCEDURE pmc_inter_bcast_integer_1
    END INTERFACE pmc_inter_bcast

    INTERFACE pmc_alloc_mem
       MODULE PROCEDURE pmc_alloc_mem_integer_1
       MODULE PROCEDURE pmc_alloc_mem_Real_1
    END INTERFACE pmc_alloc_mem

    INTERFACE pmc_time
       MODULE PROCEDURE pmc_time
    END INTERFACE pmc_time

    PUBLIC pmc_alloc_mem, pmc_bcast, pmc_inter_bcast, pmc_recv_from_child,      &
           pmc_recv_from_parent, pmc_send_to_child, pmc_send_to_parent,         &
           pmc_time

 CONTAINS


 SUBROUTINE pmc_send_to_parent_integer( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN) ::  buf          !<
    INTEGER, INTENT(IN)               ::  n            !<
    INTEGER, INTENT(IN)               ::  parent_rank  !<
    INTEGER, INTENT(IN)               ::  tag          !<
    INTEGER, INTENT(OUT)              ::  ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_INTEGER, parent_rank, tag, m_to_parent_comm,     &
                   ierr)

 END SUBROUTINE pmc_send_to_parent_integer



 SUBROUTINE pmc_recv_from_parent_integer( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(OUT) ::  buf          !<
    INTEGER, INTENT(IN)                ::  n            !<
    INTEGER, INTENT(IN)                ::  parent_rank  !<
    INTEGER, INTENT(IN)                ::  tag          !<
    INTEGER, INTENT(OUT)               ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_INTEGER, parent_rank, tag, m_to_parent_comm,     &
                   MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_parent_integer



 SUBROUTINE pmc_send_to_parent_integer_2( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: buf          !<
    INTEGER, INTENT(IN)                 :: n            !<
    INTEGER, INTENT(IN)                 :: parent_rank  !<
    INTEGER, INTENT(IN)                 :: tag          !<
    INTEGER, INTENT(OUT)                :: ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_INTEGER, parent_rank, tag, m_to_parent_comm,     &
                   ierr )

 END SUBROUTINE pmc_send_to_parent_integer_2



 SUBROUTINE pmc_send_to_parent_real_r1( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(IN) ::  buf          !<
    INTEGER, INTENT(IN)                ::  n            !<
    INTEGER, INTENT(IN)                ::  parent_rank  !<
    INTEGER, INTENT(IN)                ::  tag          !<
    INTEGER, INTENT(OUT)               ::  ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_REAL, parent_rank, tag, m_to_parent_comm, ierr )

 END SUBROUTINE pmc_send_to_parent_real_r1



 SUBROUTINE pmc_recv_from_parent_real_r1( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(OUT) ::  buf          !<
    INTEGER, INTENT(IN)                 ::  n            !<
    INTEGER, INTENT(IN)                 ::  parent_rank  !<
    INTEGER, INTENT(IN)                 ::  tag          !<
    INTEGER, INTENT(OUT)                ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_REAL, parent_rank, tag, m_to_parent_comm,        &
                   MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_parent_real_r1



 SUBROUTINE pmc_send_to_parent_real_r2( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::  buf          !<
    INTEGER, INTENT(IN)                  ::  n            !<
    INTEGER, INTENT(IN)                  ::  parent_rank  !<
    INTEGER, INTENT(IN)                  ::  tag          !<
    INTEGER, INTENT(OUT)                 ::  ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_REAL, parent_rank, tag, m_to_parent_comm, ierr )

 END SUBROUTINE pmc_send_to_parent_real_r2


 SUBROUTINE pmc_recv_from_parent_real_r2( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:), INTENT(OUT) ::  buf          !<
    INTEGER, INTENT(IN)                   ::  n            !<
    INTEGER, INTENT(IN)                   ::  parent_rank  !<
    INTEGER, INTENT(IN)                   ::  tag          !<
    INTEGER, INTENT(OUT)                  ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_REAL, parent_rank, tag, m_to_parent_comm,        &
                   MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_parent_real_r2



 SUBROUTINE pmc_send_to_parent_real_r3( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::  buf          !<
    INTEGER, INTENT(IN)                    ::  n            !<
    INTEGER, INTENT(IN)                    ::  parent_rank  !<
    INTEGER, INTENT(IN)                    ::  tag          !<
    INTEGER, INTENT(OUT)                   ::  ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_REAL, parent_rank, tag, m_to_parent_comm, ierr )

 END SUBROUTINE pmc_send_to_parent_real_r3



 SUBROUTINE pmc_recv_from_parent_real_r3( buf, n, parent_rank, tag, ierr )

    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:,:), INTENT(OUT) ::  buf          !<
    INTEGER, INTENT(IN)                     ::  n            !<
    INTEGER, INTENT(IN)                     ::  parent_rank  !<
    INTEGER, INTENT(IN)                     ::  tag          !<
    INTEGER, INTENT(OUT)                    ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_REAL, parent_rank, tag, m_to_parent_comm,        &
                   MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_parent_real_r3



 SUBROUTINE pmc_send_to_child_integer( child_id, buf, n, child_rank, tag,       &
                                       ierr )

    IMPLICIT NONE

    INTEGER, INTENT(IN)               ::  child_id     !<
    INTEGER, DIMENSION(:), INTENT(IN) ::  buf          !<
    INTEGER, INTENT(IN)               ::  n            !<
    INTEGER, INTENT(IN)               ::  child_rank   !<
    INTEGER, INTENT(IN)               ::  tag          !<
    INTEGER, INTENT(OUT)              ::  ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_INTEGER, child_rank, tag,                        &
                   m_to_child_comm(child_id), ierr )

 END SUBROUTINE pmc_send_to_child_integer



 SUBROUTINE pmc_recv_from_child_integer( child_id, buf, n, child_rank, tag,     &
                                         ierr )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                  ::  child_id     !<
    INTEGER, DIMENSION(:), INTENT(INOUT) ::  buf          !<
    INTEGER, INTENT(IN)                  ::  n            !<
    INTEGER, INTENT(IN)                  ::  child_rank   !<
    INTEGER, INTENT(IN)                  ::  tag          !<
    INTEGER, INTENT(OUT)                 ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_INTEGER, child_rank, tag,                        &
                   m_to_child_comm(child_id), MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_child_integer



 SUBROUTINE pmc_recv_from_child_integer_2( child_id, buf, n, child_rank,        &
                                           tag, ierr )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                  ::  child_id     !<
    INTEGER, DIMENSION(:,:), INTENT(OUT) ::  buf          !<
    INTEGER, INTENT(IN)                  ::  n            !<
    INTEGER, INTENT(IN)                  ::  child_rank   !<
    INTEGER, INTENT(IN)                  ::  tag          !<
    INTEGER, INTENT(OUT)                 ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_INTEGER, child_rank, tag,                        &
                   m_to_child_comm(child_id), MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_child_integer_2



 SUBROUTINE pmc_send_to_child_real_r1( child_id, buf, n, child_rank, tag,       &
                                       ierr )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                ::  child_id     !<
    REAL(wp), DIMENSION(:), INTENT(IN) ::  buf          !<
    INTEGER, INTENT(IN)                ::  n            !<
    INTEGER, INTENT(IN)                ::  child_rank   !<
    INTEGER, INTENT(IN)                ::  tag          !<
    INTEGER, INTENT(OUT)               ::  ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_REAL, child_rank, tag,                           &
                   m_to_child_comm(child_id), ierr )

 END SUBROUTINE pmc_send_to_child_real_r1



 SUBROUTINE pmc_recv_from_child_real_r1( child_id, buf, n, child_rank, tag,     &
                                         ierr )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                   ::  child_id     !<
    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  buf          !<
    INTEGER, INTENT(IN)                   ::  n            !<
    INTEGER, INTENT(IN)                   ::  child_rank   !<
    INTEGER, INTENT(IN)                   ::  tag          !<
    INTEGER, INTENT(OUT)                  ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_REAL, child_rank, tag,                           &
                   m_to_child_comm(child_id), MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_child_real_r1



 SUBROUTINE pmc_send_to_child_real_r2( child_id, buf, n, child_rank, tag,       &
                                       ierr )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                  ::  child_id     !<
    REAL(wp), DIMENSION(:,:), INTENT(IN) ::  buf          !<
    INTEGER, INTENT(IN)                  ::  n            !<
    INTEGER, INTENT(IN)                  ::  child_rank   !<
    INTEGER, INTENT(IN)                  ::  tag          !<
    INTEGER, INTENT(OUT)                 ::  ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_REAL, child_rank, tag,                           &
                   m_to_child_comm(child_id), ierr )

 END SUBROUTINE pmc_send_to_child_real_r2



 SUBROUTINE pmc_recv_from_child_real_r2( child_id, buf, n, child_rank, tag,     &
                                         ierr )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                   ::  child_id     !<
    REAL(wp), DIMENSION(:,:), INTENT(OUT) ::  buf          !<
    INTEGER, INTENT(IN)                   ::  n            !<
    INTEGER, INTENT(IN)                   ::  child_rank   !<
    INTEGER, INTENT(IN)                   ::  tag          !<
    INTEGER, INTENT(OUT)                  ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_REAL, child_rank, tag,                           &
                   m_to_child_comm(child_id), MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_child_real_r2



 SUBROUTINE pmc_send_to_child_real_r3( child_id, buf, n, child_rank, tag,       &
                                       ierr)

    IMPLICIT NONE

    INTEGER, INTENT(IN)                    ::  child_id     !<
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::  buf          !<
    INTEGER, INTENT(IN)                    ::  n            !<
    INTEGER, INTENT(IN)                    ::  child_rank   !<
    INTEGER, INTENT(IN)                    ::  tag          !<
    INTEGER, INTENT(OUT)                   ::  ierr         !<

    ierr = 0
    CALL MPI_SEND( buf, n, MPI_REAL, child_rank, tag,                           &
                   m_to_child_comm(child_id), ierr )

 END SUBROUTINE pmc_send_to_child_real_r3



 SUBROUTINE pmc_recv_from_child_real_r3( child_id, buf, n, child_rank, tag,     &
                                         ierr )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                     ::  child_id     !<
    REAL(wp), DIMENSION(:,:,:), INTENT(OUT) ::  buf          !<
    INTEGER, INTENT(IN)                     ::  n            !<
    INTEGER, INTENT(IN)                     ::  child_rank   !<
    INTEGER, INTENT(IN)                     ::  tag          !<
    INTEGER, INTENT(OUT)                    ::  ierr         !<

    ierr = 0
    CALL MPI_RECV( buf, n, MPI_REAL, child_rank, tag,                           & 
                   m_to_child_comm(child_id), MPI_STATUS_IGNORE, ierr )

 END SUBROUTINE pmc_recv_from_child_real_r3



 SUBROUTINE pmc_bcast_integer( buf, root_pe, comm, ierr )

    IMPLICIT NONE

    INTEGER, INTENT(INOUT)         ::  buf      !<
    INTEGER, INTENT(IN)            ::  root_pe  !<
    INTEGER, INTENT(IN), OPTIONAL  ::  comm     !<
    INTEGER, INTENT(OUT), OPTIONAL ::  ierr     !<

    INTEGER ::  mycomm  !<
    INTEGER ::  myerr   !<


    IF ( PRESENT( comm ) )  THEN
       mycomm = comm
    ELSE
       mycomm = m_model_comm
    ENDIF

    CALL MPI_BCAST( buf, 1, MPI_INTEGER, root_pe, mycomm, myerr )

    IF ( PRESENT( ierr ) )  THEN
       ierr = myerr
    ENDIF

 END SUBROUTINE pmc_bcast_integer



 SUBROUTINE pmc_bcast_character( buf, root_pe, comm, ierr )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(INOUT) ::  buf      !<
    INTEGER, INTENT(IN)             ::  root_pe  !<
    INTEGER, INTENT(IN), OPTIONAL   ::  comm     !<
    INTEGER, INTENT(OUT), OPTIONAL  ::  ierr     !<

    INTEGER ::  mycomm  !<
    INTEGER ::  myerr   !<

    IF ( PRESENT( comm ) )  THEN
       mycomm = comm
    ELSE
       mycomm = m_model_comm
    ENDIF

    CALL MPI_BCAST( buf, LEN(buf), MPI_CHARACTER, root_pe, mycomm, myerr )

    IF ( PRESENT( ierr ) )  THEN
       ierr = myerr
    ENDIF

 END SUBROUTINE pmc_bcast_character



 SUBROUTINE pmc_inter_bcast_integer_1( buf, child_id, ierr )

    IMPLICIT NONE

    INTEGER, INTENT(INOUT),DIMENSION(:) ::  buf        !<
    INTEGER, INTENT(IN),optional        ::  child_id   !<
    INTEGER, INTENT(OUT),optional       ::  ierr       !<

    INTEGER ::  mycomm   !<
    INTEGER ::  myerr    !<
    INTEGER ::  root_pe  !<

!
!-- PE 0 on parent broadcast to all child PEs
    IF ( PRESENT( child_id ) )  THEN

       mycomm = m_to_child_comm(child_id)

       IF ( m_model_rank == 0 )  THEN
          root_pe = MPI_ROOT
       ELSE
          root_pe = MPI_PROC_NULL
       ENDIF

    ELSE
       mycomm  = m_to_parent_comm
       root_pe = 0
    ENDIF

    CALL MPI_BCAST( buf, SIZE( buf ), MPI_INTEGER, root_pe, mycomm, myerr )

    IF ( PRESENT( ierr ) )  THEN
       ierr = myerr
    ENDIF

 END SUBROUTINE pmc_inter_bcast_integer_1



 SUBROUTINE pmc_alloc_mem_integer_1( iarray, idim1 )
!
!-- Allocate memory with MPI_ALLOC_MEM using intermediate C-pointer

    IMPLICIT NONE

    INTEGER, DIMENSION(:), POINTER, INTENT(INOUT) ::  iarray  !<
    INTEGER, INTENT(IN)                           ::  idim1   !<

    INTEGER, DIMENSION(1)          ::  ashape   !<
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !<
    INTEGER                        ::  ierr     !<

    TYPE(C_PTR)                    ::  p_myind  !<

    winsize = idim1 * C_SIZEOF( ierr )

    CALL MPI_ALLOC_MEM( winsize, MPI_INFO_NULL, p_myind, ierr )
    ashape(1) = idim1
    CALL C_F_POINTER( p_myind, iarray, ashape )

 END SUBROUTINE pmc_alloc_mem_integer_1



 SUBROUTINE pmc_alloc_mem_real_1( array, idim1, base_ptr )

    IMPLICIT NONE

    INTEGER(idp), INTENT(IN)                            ::  idim1     !<
    REAL(KIND=wp), DIMENSION(:), POINTER, INTENT(INOUT) ::  array     !<
    TYPE(C_PTR), INTENT(OUT), OPTIONAL                  ::  base_ptr  !<

    INTEGER, DIMENSION(1)          :: ashape   !<
    INTEGER(KIND=MPI_ADDRESS_KIND) :: winsize  !<
    INTEGER                        :: ierr     !<

    TYPE(C_PTR)                    :: p_myind  !<

    winsize = idim1 * wp

    CALL MPI_ALLOC_MEM( winsize , MPI_INFO_NULL, p_myind, ierr )
    ashape(1) = idim1
    CALL C_F_POINTER( p_myind, array, ashape )

    IF ( PRESENT( base_ptr ) )  THEN
       base_ptr = p_myind
    ENDIF

 END SUBROUTINE pmc_alloc_mem_Real_1



 FUNCTION pmc_time()

    REAL(kind=wp) :: pmc_time  !<

    pmc_time = MPI_WTIME()

  END FUNCTION pmc_time

#endif
 END MODULE pmc_mpi_wrapper
