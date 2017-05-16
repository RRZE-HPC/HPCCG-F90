! =======================================================================================
!
!      Filename:  SparseMatrix.F90
!      Description:  Fortran90 module implementing CRS data structure and
!                    subroutines. Also contains spmv operation.
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

module SparseMatrix
    use constants
#ifdef _MPI
    use mpi
#endif
#ifndef HAVE_COOP
    use caf_collectives
#endif

    implicit none

    type, public :: SparseMatrixT
        integer :: start_row
        integer :: stop_row
        integer :: total_nrow
        integer(kind=4) :: total_nnz
        integer :: local_nrow
        integer :: local_ncol
        integer :: local_nnz
        integer, allocatable ::  colindex(:)
        integer, allocatable ::  rowindex(:)
        real(kind=dp), allocatable :: values(:)
        integer :: global_max_ncol
    end type SparseMatrixT

contains

    subroutine generate_matrix(nx, ny, nz, A, x, b, xexact)
        implicit none
        integer, intent(in) :: nx, ny, nz
        type(SparseMatrixT), intent(out) :: A
        real(kind=dp), allocatable, intent(out) :: x(:)
        real(kind=dp), allocatable, intent(out) :: b(:)
        real(kind=dp), allocatable, intent(out) :: xexact(:)

        real(kind=dp) :: local_result
        integer :: ix, iy, iz
        integer :: sx, sy, sz
        integer :: mysize, myrank
        integer :: local_nrow
        integer :: local_nnz
        integer :: total_nrow
        integer(kind=4) :: total_nnz
        integer :: start_row
        integer :: stop_row
        integer :: currow, curcol
        integer :: curlocalrow
        integer :: nnzrow
        integer :: cursor
        integer :: rowcursor

        mysize = num_images()
        myrank = this_image() - 1
        cursor = 1
        rowcursor = 1
        local_nrow = nx*ny*nz
        local_nnz = 27 * local_nrow
        total_nrow = local_nrow * mysize
        total_nnz = 27 * total_nrow
        start_row = (local_nrow * myrank) + 1
        stop_row = start_row + local_nrow - 1

        allocate(A%values(local_nnz))
        allocate(A%rowindex(local_nrow+1))
        allocate(A%colindex(local_nnz))
        allocate(x(local_nrow))
        allocate(b(local_nrow))
        allocate(xexact(local_nrow))

        do iz = 1, nz
            do iy = 1, ny
                do ix = 1, nx
                    curlocalrow = (iz-1)*nx*ny + (iy-1)*nx + (ix-1)
                    currow = start_row + curlocalrow
                    nnzrow = 0
                    A%rowindex(rowcursor) = cursor

                    do sz = -1, 1
                        do sy = -1, 1
                            do sx = -1, 1
                                curcol = currow + sz*nx*ny + sy*nx + sx;

                                if ( (ix+sx >= 1) .and. (ix+sx <= nx) .and. &
                                    (iy+sy >= 1) .and. (iy+sy <= ny) .and. &
                                    (curcol>= 1  .and. curcol<=total_nrow)) then

                                if (curcol == currow) then
                                    A%values(cursor) = 27.0
                                else
                                    A%values(cursor) = -1.0
                                end if

                                A%colindex(cursor) = curcol
                                cursor = cursor + 1
                                nnzrow = nnzrow + 1
                            end if
                        end do !sx
                    end do !sy
                end do !sz

                rowcursor = rowcursor + 1
                x(curlocalrow+1) = 0.0
                b(curlocalrow+1) = 27.0 - (nnzrow-1)
                xexact(curlocalrow+1) = 1.0
            end do !ix
        end do !iy
    end do !iz

    A%rowindex(rowcursor) = cursor
    A%start_row = start_row
    A%stop_row = stop_row
    A%total_nrow = total_nrow
    A%total_nnz = total_nnz
    A%local_nrow = local_nrow
    A%local_ncol = local_nrow
    A%local_nnz = local_nnz

    A%global_max_ncol = A%local_ncol
#ifdef HAVE_COOP
    call CO_MAX(A%global_max_ncol)
#else
    call CAF_CO_MAX(A%global_max_ncol)
#endif
end subroutine generate_matrix

subroutine dump_matrix(nx, ny, nz, A, x, b)
    implicit none
    integer, intent(in) :: nx, ny, nz
    type(SparseMatrixT), intent(in) :: A
    real(kind=dp), allocatable, intent(in) :: x(:)
    real(kind=dp), allocatable, intent(in) :: b(:)

    integer :: i, k
    integer :: numEntries
    real(kind=dp) :: test

    print *,'Matrix has ',A%total_nrow, &
        ' total rows and ',A%local_ncol,' columns.'
    print *,'Matrix has ',A%local_nrow,' local rows.'

    do i = 1, A%local_nrow
        numEntries = (A%rowindex(i+1)-1) - A%rowindex(i)+1
        print *,'Row id',i,'Start',A%rowindex(i), &
            'Stop',A%rowindex(i+1)-1,'Nonzeroes',numEntries

        do k = A%rowindex(i), A%rowindex(i+1)-1
            print "(i,a,f4.1)",A%colindex(k),' : ',A%values(k)
        end do !k
    end do !i
end subroutine dump_matrix

subroutine spmv(A, x, y)
    implicit none
    type(SparseMatrixT), intent(in) :: A
    real(kind=dp), intent(in) :: x(:)[*]
    real(kind=dp), intent(inout) :: y(:)

    integer :: i, k
    integer :: id, image
    real(kind=dp) :: rowsum

    sync all

!$OMP PARALLEL DO
    do i = 1, A%local_nrow
        rowsum = 0.0D+00

        do k = A%rowindex(i), A%rowindex(i+1)-1
            image = ((A%colindex(k)-1)/A%local_nrow) + 1
            id = MODULO(A%colindex(k)-1,A%local_nrow) + 1
            rowsum = rowsum + A%values(k) * x(id)[image]
        end do !k

        y(i) = rowsum
    end do !i

    sync all
end subroutine spmv

end module SparseMatrix

