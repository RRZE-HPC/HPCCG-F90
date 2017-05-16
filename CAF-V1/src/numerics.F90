! =======================================================================================
!
!      Filename:  numerics.F90
!      Description:  Implements the main cg solver iteration loop
!      Date: Jun 24, 2014
!      Author:  Jan Eitzinger (je), jan.eitzinger@fau.de
!      Copyright (c) 2017, Jan Eitzinger
!      All rights reserved.
!
!      Redistribution and use in source and binary forms, with or without
!      modification, are permitted provided that the following conditions are met:
!
!      * Redistributions of source code must retain the above copyright notice, this
!        list of conditions and the following disclaimer.
!
!      * Redistributions in binary form must reproduce the above copyright notice,
!        this list of conditions and the following disclaimer in the documentation
!        and/or other materials provided with the distribution.
!
!      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!      ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!      WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!      DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!      FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!      DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!      SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!      CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!      OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!=======================================================================================

#include <timer.h>

module numerics
contains
    subroutine cg(A, b, x, max_iter, tolerance, niters, normr)
        use constants
        use SparseMatrix
        use vector_operations
        use timer
#ifdef _MPI
        use mpi
#endif

        implicit none
        type(SparseMatrixT), intent(inout) :: A
        real(kind=dp), allocatable, intent(inout) :: b(:)
        real(kind=dp), allocatable, intent(inout) :: x(:)
        integer, intent(in) :: max_iter
        real(kind=dp), intent(in) :: tolerance
        integer, intent(out) :: niters
        real(kind=dp), intent(out) :: normr
        integer :: nrow
        integer :: ncol
        integer :: myrank
        real(kind=dp) :: alpha, beta
        real(kind=dp) :: rtrans
        real(kind=dp) :: oldrtrans
        real(kind=dp), allocatable :: r(:)
        real(kind=dp), allocatable :: p(:)[:]
        real(kind=dp), allocatable :: Ap(:)
        integer ::  k, i, j
        integer :: print_freq
        type(interval) :: timerdata, timertotal, timerspmv
        real(kind=dp) ::fnops_sparsemv, total_time_spmv = 0.0d0

        print_freq = max_iter/10
        if (print_freq>50) then
            print_freq=50
        end if
        if (print_freq<1) then
            print_freq=1
        end if

        myrank = this_image() - 1
        nrow = A%local_nrow
        ncol = A%local_ncol
        allocate(r(nrow))
        allocate(p(ncol)[*])
        allocate(Ap(nrow))

        normr = 0.0
        rtrans = 0.0
        oldrtrans = 0.0

        TIMER_START(1, timerdata,"Prepare")
        call vcopy(nrow, x, p)
        call spmv(A, p, Ap)
        call waxpby(nrow, 1.0d0, b, -1.0d0, Ap, r)
        call ddot1(nrow, r, rtrans)
        normr = sqrt(rtrans);
        TIMER_STOP(1, timerdata,"Prepare")
        TIMER_PRINT(1, timerdata,"Prepare")
        if (myrank == 0) then
        print *,"Initial Residual = ",normr
        end if

        TIMER_START(0, timertotal,"Total")

        do k = 1, max_iter
            if (normr < tolerance) then
                print *,"Converged Solution"
                exit
            end if
            TIMER_START(1, timerdata,"Iteration")
            if (k == 1) then
                call vcopy(nrow, r, p)
            else
                oldrtrans = rtrans
                call ddot1 (nrow, r, rtrans)
                beta = rtrans/oldrtrans
                call wxbw (nrow,  r, beta, p)
            end if

            normr = sqrt(rtrans)
            if (myrank==0 .and. (MOD(k,print_freq) == 0 .or. k+1 == max_iter)) then
                print "(a,i5,tr2,a,e20.10)",'ITER',k,'Residual',normr
            end if
            TIMER_START(2, timerspmv,"spmv")
            call spmv(A, p, Ap)
            TIMER_STOP(2, timerspmv,"spmv")
!            total_time_spmv = total_time_spmv + timer_print(timerspmv)
            alpha = 0.0
            call ddot(nrow, p, Ap, alpha)
            alpha = rtrans/alpha
            call wwpby(nrow, -alpha, Ap, r)
            call wwpby(nrow,  alpha, p, x)
            niters = k
            TIMER_STOP(1, timerdata,"Iteration")
            TIMER_PRINT(1, timerdata,"Iteration")
        end do

        TIMER_STOP(0, timertotal,"Total")
        TIMER_PRINT(0, timertotal,"Total")

        if (TIMINGLEV > 2) then
            fnops_sparsemv = niters*2.0d0*A%total_nnz
            print *,"MFLOPS SPARSEMV", &
                fnops_sparsemv/total_time_spmv/1.0d6
        end if

        deallocate(r)
        deallocate(p)
        deallocate(Ap)

    end subroutine cg

    subroutine compute_residual(n, v1, v2, residual)

        use constants
        use mpi

        implicit none
        integer, intent(in) :: n
        real(kind=dp), allocatable, intent(inout) :: v1(:)
        real(kind=dp), allocatable, intent(inout) :: v2(:)
        real(kind=dp), intent(out) :: residual
#ifdef _MPI
        real(kind=dp) :: global_residual
#endif
        real(kind=dp) :: local_residual = 0.0d0
        real(kind=dp) :: diff
        integer :: i, ierror

        do i = 1, n
            diff = abs(v1(i)-v2(i))
            if ( diff .gt. local_residual) then
                local_residual = diff
            endif
        end do

#ifdef _MPI
        global_residual = 0

        call MPI_ALLREDUCE(local_residual,global_residual,1,  &
            MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
        residual = global_residual
#else
!        call CO_MAX(local_residual)
        residual = local_residual
#endif

    end subroutine compute_residual

end module numerics

