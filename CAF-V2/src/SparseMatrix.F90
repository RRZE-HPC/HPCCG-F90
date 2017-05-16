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
        integer(kind=lk) :: total_nnz
        integer :: local_nrow
        integer :: local_ncol
        integer :: local_nnz
        integer, allocatable ::  colindex(:)
        integer, allocatable ::  rowindex(:)
        real(kind=dp), allocatable :: values(:)
        integer :: num_external, num_blocks
        integer :: global_max_ncol
        integer, allocatable :: remote_index(:)
!        integer, allocatable :: remote_image(:)
        integer, dimension(max_blocks+1) :: image, inmin, inmax, start
    end type SparseMatrixT

    real(kind=dp), allocatable :: xCA(:)[:]
    real(kind=dp) :: xl(max_external)
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
        integer(kind=lk) :: total_nnz
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
        local_nnz = 27_lk * local_nrow
        total_nrow = local_nrow * mysize
        total_nnz = 27_lk * total_nrow
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
            print "(i0,a,f4.1)",A%colindex(k),' : ',A%values(k)
        end do !k
    end do !i
end subroutine dump_matrix

subroutine spmv(A, x, y)
    implicit none
    type(SparseMatrixT), intent(in) :: A
    real(kind=dp), allocatable, intent(in) :: x(:)
    real(kind=dp), allocatable, intent(inout) :: y(:)

    integer :: i, k
    real(kind=dp) :: rowsum

!$OMP PARALLEL DO
    do i = 1, A%local_nrow
        rowsum = 0.0D+00

        do k = A%rowindex(i), A%rowindex(i+1)-1
            rowsum = rowsum + A%values(k) * x(A%colindex(k))
        end do !k

        y(i) = rowsum
    end do !i
end subroutine spmv

subroutine make_local_matrix(A)
    implicit none
    type(SparseMatrixT), intent(inout) :: A

    integer :: mysize, myrank
    integer :: i, j, img, ind
    integer :: cur_index
    integer :: counter
    integer :: num_external
    integer, allocatable :: externals(:)
    integer, allocatable :: external_index(:)

    ! Convert the index values for the rows on this processor
    ! to a local index space.
    ! We need to:
    ! * Determine if each index reaches to a local value or external value
    ! * If local, subtract start_row from index value to get local index
    ! * If external, find out if it is already accounted for
    !   - If so, then do nothing
    !   - otherwise
    !     - add it to the list of external indices
    !     - find out which processor owns the value
    !     - Set up communication for sparse MV operation

    mysize = num_images()
    myrank = this_image() - 1

    ! for the moment instead of a map use a simple global lookup table
    allocate(externals(A%total_nrow))
    allocate(external_index(max_external))
    A%image(:) = 0
    A%inmin(:) = HUGE(1)
    A%inmax(:) = 0

    externals = -1
    num_external = 1
    A%num_blocks = 1

    ! Scan the indices and transform to local
    do i = 1, A%local_nrow
        do j =  A%rowindex(i), A%rowindex(i+1)-1
            cur_index = A%colindex(j)

            ! shift local rows to the start
            if (A%start_row <= cur_index .and. cur_index <= A%stop_row) then
                A%colindex(j) = A%colindex(j) - A%start_row + 1
            else
                if (externals(cur_index) == -1) then
                    externals(cur_index) = num_external


            img = ((cur_index-1)/A%local_nrow) + 1 
            if (A%image(A%num_blocks) == 0) then
               A%image(A%num_blocks) = img
               A%start(A%num_blocks) = num_external
            end if
            if (img /= A%image(A%num_blocks)) then
               A%num_blocks = A%num_blocks + 1
               if (A%num_blocks > max_blocks) stop 'INCREASE max_blocks'
               A%image(A%num_blocks) = img
               A%start(A%num_blocks) = num_external
            end if
            
            ind =  MODULO(cur_index-1,A%local_nrow) + 1
            A%inmin(A%num_blocks) = min(ind, A%inmin(A%num_blocks))
            A%inmax(A%num_blocks) = max(ind, A%inmax(A%num_blocks))
     

                    num_external = num_external + 1

                    if (num_external <= max_external) then
                        external_index(num_external-1) = cur_index
                        ! Mark index as external by negating it
                        A%colindex(j) = -A%colindex(j)
                    else
                        print *,"Must increase max_external"
                        stop
                    end if
                else
                    ! Mark index as external by negating it
                    A%colindex(j) = -A%colindex(j)
                end if
            end if
        end do !j
    end do !i

    A%start(A%num_blocks+1) = num_external
    A%num_external = num_external-1
    if (num_images() == 1) A%num_blocks = 0

!    allocate(A%remote_image(A%num_external))
    allocate(A%remote_index(A%num_external))

    !Go through list of externals and find the processor that owns each
    !Keep it simple for the moment and use the static partitioning in HPCCG
    !Also setup the mapping between external id to remote id here.

    do i = 1, A%num_external
        cur_index = external_index(i)
!        A%remote_image(i) = ((cur_index-1)/A%local_nrow) + 1 ! one based index
        A%remote_index(i) =  MODULO(cur_index-1,A%local_nrow) + 1
    end do !i

    !In contrast to the MPI implementation the entries will not be sorted blockwise
    !according to owning processor but keep their ordering.

    ! map all external ids to the new local index
    do i = 1, A%local_nrow
        do j =  A%rowindex(i), A%rowindex(i+1)-1
            if ( A%colindex(j)<0 ) then
                cur_index = -A%colindex(j)
                A%colindex(j) = A%local_nrow + externals(cur_index)
            end if
        end do !j
    end do !i

    !! Finish up
    A%local_ncol = A%local_nrow + A%num_external
    A%global_max_ncol = A%local_nrow
#ifdef HAVE_COOP
    call CO_MAX(A%global_max_ncol)
#else
    call CAF_CO_MAX(A%global_max_ncol)
#endif
    deallocate(externals)
    allocate(xCA(A%global_max_ncol)[*])
end subroutine make_local_matrix

subroutine exchange_externals(A,x)
    implicit none
    type(SparseMatrixT), intent(in) :: A
    real(kind=dp), intent(inout) :: x(:)

    integer :: i, k

!    write(*,*) 'DEBUG: ',size(xca), A%global_max_ncol, size(x)
    xCA(1:A%global_max_ncol) = x(1:A%global_max_ncol)
    sync all

    ! Update local halo by pulling from remote images
!    if (this_image() == 3) then
!      write(*,*) '#ext: ',A%num_external
!      write(*,*) '#ext_remoteind: ',A%remote_index
!      write(*,*) '#ext_remoteimg: ',A%remote_image
!    end if
!    do i = 1, A%num_external
!        x(A%local_nrow+i) = xCA(A%remote_index(i))[A%remote_image(i)]
!    end do !i
    if (this_image() == 0) then
      write(*,*) '#ext: ',A%num_blocks, A%num_external
      write(*,*) '#ext_min: ',A%inmin(1:A%num_blocks)
      write(*,*) '#ext_max: ',A%inmax(1:A%num_blocks)
      write(*,*) '#ext_img: ',A%image(1:A%num_blocks)
      write(*,*) '#ext_start: ',A%start(1:A%num_blocks+1)
    end if
     do k = 1, A%num_blocks
        xl(1:A%inmax(k)-A%inmin(k)+1) = xCA(A%inmin(k):A%inmax(k))[A%image(k)]
        do i = A%start(k), A%start(k+1)-1
!           if (A%remote_index(i) < A%inmin(k) .or. &
!               A%remote_index(i) > A%inmax(k)) STOP 'out of bounds'
           x(A%local_nrow+i) = xl(A%remote_index(i)-A%inmin(k)+1)
        end do
    end do 

    sync all
end subroutine exchange_externals

end module SparseMatrix

