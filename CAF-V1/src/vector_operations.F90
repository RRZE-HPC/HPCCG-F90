! =======================================================================================
!
!      Filename:  vector_operations.F90
!      Description:  Vector operation necessary for CG solver
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

module vector_operations
    use constants
#ifdef _MPI
    use mpi
#endif
#ifndef HAVE_COOP
    use caf_collectives
#endif

contains
    subroutine ddot(n, x, y, r)
        implicit none
        integer, intent(in) :: n
        real(kind=dp),  intent(in) :: x(:)
        real(kind=dp),  intent(in) :: y(:)
        real(kind=dp),  intent(inout) :: r

        real(kind=dp) :: local_result
        integer :: i
#ifdef _MPI
        real(kind=dp) :: global_result
        integer :: ierror
#endif
        local_result = 0.0

        !$OMP PARALLEL DO REDUCTION (+:local_result)
        do i = 1, n
            local_result = local_result + x(i) * y(i)
        end do

#ifdef _MPI
        global_result = 0.0
        call MPI_ALLREDUCE(local_result,global_result,1,  &
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
        r = global_result
#else
#ifdef HAVE_COOP
        call CO_SUM(local_result)
#else
        call CAF_CO_SUM(local_result)
#endif
        r = local_result
#endif
    end subroutine ddot

    subroutine ddot1(n, x, r)
        implicit none
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(inout) :: r

        real(kind=dp) :: local_result
        integer :: i
#ifdef _MPI
        integer :: ierror
        real(kind=dp) :: global_result
#endif
        local_result = 0.0

        !$OMP PARALLEL DO REDUCTION (+:local_result)
        do i = 1, n
            local_result = local_result + x(i) * x(i)
        end do

#ifdef _MPI
        global_result = 0.0
        call MPI_ALLREDUCE(local_result,global_result,1,  &
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
        r = global_result
#else
#ifdef HAVE_COOP
        call CO_SUM(local_result)
#else
        call CAF_CO_SUM(local_result)
#endif
        r = local_result
#endif
    end subroutine ddot1

    subroutine waxpby(n, alpha, x, beta, y, w)
        implicit none
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: alpha
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(in) :: beta
        real(kind=dp), intent(in) :: y(:)
        real(kind=dp), intent(inout) :: w(:)

        integer :: i

        !$OMP PARALLEL DO
        do i = 1, n
            w(i) = alpha * x(i) + beta * y(i)
        end do

    end subroutine waxpby

    subroutine wwpby(n, beta, y, w)
        implicit none
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: y(:)
        real(kind=dp), intent(in) :: beta
        real(kind=dp), intent(inout) :: w(:)

        integer :: i

        !$OMP PARALLEL DO
        do i = 1, n
            w(i) = w(i) + beta * y(i)
        end do

    end subroutine wwpby

    subroutine wxbw(n,  x, beta,  w)
        implicit none
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(in) :: beta
        real(kind=dp), intent(inout) :: w(:)

        integer :: i

        !$OMP PARALLEL DO
        do i = 1, n
            w(i) =  x(i) + beta * w(i)
        end do

    end subroutine wxbw

    subroutine vcopy(n, x, y)

        implicit none
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(inout) :: y(:)

        integer :: i

        !$OMP PARALLEL DO
        do i = 1, n
            y(i) =  x(i)
        end do

    end subroutine vcopy

end module vector_operations

