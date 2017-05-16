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

    implicit none
    integer, public :: ierror
    integer, public :: single  ! indicate if MPI program runs with 1 process

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
#ifdef _MPI
        integer :: num_send_neighbors
        integer :: total_to_be_sent
        integer, allocatable :: elements_to_send(:)
        integer, allocatable :: neighbors(:)
        integer, allocatable :: recv_length(:)
        integer, allocatable :: send_length(:)
        real(kind=dp), allocatable :: send_buffer(:)
#endif
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
#ifdef _MPI
        real(kind=dp) :: global_result

        call MPI_COMM_SIZE(MPI_COMM_WORLD,mysize,ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)

        if (mysize == 1) then
            single = 1
        else
            single = 0
        endif
#else
        mysize = 1
        myrank = 0
#endif
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

!TODO: At the moment this code does not perform first touch initialization.
! For scaling across ccNUMA domain boundaries use hybrid variant with 1 process
! per domain.
!#ifdef _OMP
!!$OMP PARALLEL DO
!        do iz = 1, A%local_nrow
!            ix = (iz-1) * 27
!            do iy = 1, 27
!                A%values(ix + iy) = 0.0
!                A%colindex(ix + iy) = 0.0
!            end do !iy
!            x(iz) = 0.0
!            b(iz) = 0.0
!        end do !iz
!#endif

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
            print "(i,a,f4.1)",A%colindex(k),' : ',A%values(k)
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

#ifdef _MPI
subroutine make_local_matrix(A)
    implicit none
    type(SparseMatrixT), intent(inout) :: A

    integer :: mysize, myrank
    integer :: i, j
    integer :: cur_index
    integer :: counter
    integer, allocatable :: tmp_buffer(:)
    integer, allocatable :: tmp_neighbors(:)
    integer :: num_external
    integer, allocatable :: externals(:)
    integer, allocatable :: external_processor(:)
    integer, allocatable :: external_index(:)
    integer, allocatable :: external_local_index(:)
    integer, allocatable :: new_external_processor(:)
    integer, allocatable :: global_index_offsets(:)
    integer :: num_recv_neighbors  !number of ranks we receive from
    integer :: num_send_neighbors  !number of ranks we need to send to
    integer :: length
    integer :: found
    integer :: MPI_MY_TAG
    integer, allocatable :: recv_list(:) !list with all ranks we reveive from
    integer, allocatable :: send_list(:) !list with all ranks we send to
    integer, allocatable :: request(:)
    integer :: mystatus(MPI_STATUS_SIZE)
    integer, allocatable :: new_external(:)
    integer, allocatable :: lengths(:)
    integer :: start, newlength

    if (single == 1) then
        print *, "RUN with 1 PROCESS"
        return
    end if

    ! We need to convert the index values for the rows on this processor
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

    call MPI_COMM_SIZE(MPI_COMM_WORLD,mysize,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)

    ! for the moment instead of a map use a simple global lookup table
    allocate(externals(A%total_nrow))
    allocate(external_index(max_external))

    externals = -1
    num_external = 1

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

    num_external = num_external-1
    allocate(tmp_buffer(mysize))
    allocate(global_index_offsets(mysize))

    !Go through list of externals to find out which processors must be accessed.
    tmp_buffer = 0
    tmp_buffer(myrank+1) = A%start_row

    call MPI_ALLREDUCE(tmp_buffer, global_index_offsets, mysize, MPI_INTEGER, &
        MPI_SUM, MPI_COMM_WORLD,ierror)

    allocate(external_processor(num_external))

    !Go through list of externals and find the processor that owns each
    do i = 1, num_external
        cur_index = external_index(i)

        do j = mysize, 1, -1
            if (global_index_offsets(j) <= cur_index) then
                external_processor(i) = j
                exit
            end if
        end do !j
    end do !i

    !Go through the external elements. For each newly encountered external
    !point assign it the next index in the local sequence. Then look for other
    !external elements who are updated by the same node and assign them the next
    !set of index numbers in the local sequence (ie. elements updated by the same node
    !have consecutive indices).

    allocate(external_local_index(max_external))
    external_local_index = -1
    counter = A%local_nrow + 1

    do i = 1, num_external
        if (external_local_index(i) == -1) then
            external_local_index(i) = counter
            counter = counter + 1

            do j = i+1, num_external
                if (external_processor(j) == external_processor(i)) then
                    external_local_index(j) = counter
                    counter = counter + 1
                end if
            end do !j
        end if
    end do !i

    ! map all external ids to the new local index
    do i = 1, A%local_nrow
        do j =  A%rowindex(i), A%rowindex(i+1)-1
            if ( A%colindex(j)<0 ) then
                cur_index = -A%colindex(j)
                A%colindex(j) = external_local_index(externals(cur_index))
            end if
        end do !j
    end do !i

    allocate(new_external_processor(num_external))
    new_external_processor = 0

    ! setup map from external id to partition
    do i = 1, num_external
        new_external_processor(external_local_index(i) - A%local_nrow) = &
            external_processor(i)
    end do !i

    ! Count the number of neighbors from which we receive information to update
    ! our external elements. Additionally, fill the array tmp_neighbors in the
    ! following way:
    !      tmp_neighbors[i] = 0   ==>  No external elements are updated by
    !                                  processor i
    !      tmp_neighbors[i] = x   ==>  (x-1)/size elements are updated from
    !                                  processor i
    allocate(tmp_neighbors(mysize))
    tmp_neighbors = 0

    num_recv_neighbors = 0
    length             = 1

    do i = 1, num_external
        if (tmp_neighbors(new_external_processor(i)) == 0) then
            num_recv_neighbors = num_recv_neighbors + 1
            tmp_neighbors(new_external_processor(i)) = 1
        end if
        tmp_neighbors(new_external_processor(i)) = &
            tmp_neighbors(new_external_processor(i)) + mysize
    end do

    ! sum over all processors all the tmp_neighbors arrays
    call MPI_ALLREDUCE(tmp_neighbors, tmp_buffer, mysize, MPI_INTEGER, MPI_SUM, &
        MPI_COMM_WORLD, ierror)

    ! decode the combined 'tmp_neighbors' (stored in tmp_buffer)
    ! array from all the processors
    num_send_neighbors = MOD(tmp_buffer(myrank+1),mysize)

    ! decode 'tmp_buffer[rank] to deduce total number of elements we must send
    A%total_to_be_sent = (tmp_buffer(myrank+1) - num_send_neighbors) / mysize

    ! Check to see if we have enough workspace allocated.  This could be
    ! dynamically modified, but let's keep it simple for now...
    if (num_send_neighbors > max_num_messages) then
        print *,"Must increase max_num_messages in constants.F90"
        print *, "Must be at least ",num_send_neighbors
        stop
    end if

    if (A%total_to_be_sent > max_external ) then
        print *,"Must increase max_external in constants.F90"
        print *, "Must be at least ",A%total_to_be_sent
        stop
    end if

    deallocate(tmp_neighbors)

#ifdef DEBUG
    print *,"Processor ",myrank," of ",mysize,": &
        Number of send neighbors = ",num_send_neighbors
    print *,"Processor ",myrank," of ",mysize,": &
        Number of receive neighbors = ",num_recv_neighbors
    print *,"Processor ",myrank," of ",mysize,": &
        Total number of elements to send = ",A%total_to_be_sent
#endif

    call MPI_Barrier(MPI_COMM_WORLD, ierror)

    !Make a list of the neighbors that will send information to update our
    !external elements (in the order that we will receive this information)
    allocate(recv_list(max_external))

    j = 1
    recv_list(j) = new_external_processor(1) - 1
    j = j + 1

    do i = 2, num_external
        if (new_external_processor(i-1) /= new_external_processor(i)) then
            recv_list(j) = new_external_processor(i) - 1
            j = j + 1
        end if
    end do !i

    ! Ensure that all the neighbors we expect to receive from also send to us
    ! Send a 0 length message to each of our recv neighbors
    allocate(send_list(num_send_neighbors))
    send_list = 0

    MPI_MY_TAG = 99
    allocate(request(max_num_messages))

    do i = 1, num_send_neighbors
        call MPI_IRECV(tmp_buffer(i), 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_MY_TAG, &
            MPI_COMM_WORLD, request(i),ierror)
    end do !i

    do i = 1, num_recv_neighbors
        call MPI_SEND(tmp_buffer(i), 1, MPI_INTEGER, recv_list(i), MPI_MY_TAG, &
            MPI_COMM_WORLD, ierror)
    end do !i

    ! Receive message from each send neighbor to construct 'send_list'.
    do i = 1, num_send_neighbors
        call MPI_WAIT(request(i), mystatus, ierror)
        send_list(i) = mystatus(MPI_SOURCE)
    end do !i

    ! Compare the two lists. In most cases they should be the same.
    ! However, if they are not then add new entries to the recv list
    ! that are in the send list (but not already in the recv list).
    !WHY!! This ensures that the sendlist is equal to the sendlist
    !But why is this required? -> Just One neighbour list??
    do j = 1, num_send_neighbors
        found = 0;
        do i = 1, num_recv_neighbors
            if (recv_list(i) == send_list(j)) then
                found = 1
            end if
        end do !i

        if (found == 0) then
            recv_list(num_recv_neighbors) = send_list(j)
            num_recv_neighbors = num_recv_neighbors + 1
        end if
    end do !j

    deallocate(send_list)
    num_send_neighbors = num_recv_neighbors

    if (num_send_neighbors > max_num_messages) then
        print *,"Must increase max_external in constants.F90"
        stop
    end if

    ! Start filling SparseMatrix struct
    ! Create 'new_external' which explicitly put the external elements in the
    ! order given by 'external_local_index'
    allocate(new_external(num_external))

    do i = 1, num_external
        new_external(external_local_index(i) - A%local_nrow) = external_index(i)
    end do !i

    deallocate(external_local_index)
    deallocate(external_index)
    ! Send each processor the global index list of the external elements in the
    ! order that I will want to receive them when updating my external elements
    allocate(lengths(num_recv_neighbors))
    MPI_MY_TAG = MPI_MY_TAG + 1

    ! First post receives
    do i = 1, num_recv_neighbors
        call MPI_IRECV(lengths(i), 1, MPI_INTEGER, recv_list(i) , MPI_MY_TAG,  &
            MPI_COMM_WORLD, request(i), ierror)
    end do !i

    allocate(A%neighbors(max_num_neighbors))
    allocate(A%recv_length(max_num_neighbors))
    allocate(A%send_length(max_num_neighbors))

    j = 1
    do i = 1, num_recv_neighbors
        start  = j
        newlength = 0

        ! go through list of external elements until updating processor changes
        do while ((j <= num_external) .and. &
                (new_external_processor(j) == (recv_list(i)+1)))
            newlength = newlength + 1
            j = j + 1
            if (j == num_external) exit !is this necessary ??
        end do

        A%recv_length(i) = newlength + 1
        A%neighbors(i)  = recv_list(i)
        length = j - start + 1
        call MPI_SEND(length, 1, MPI_INTEGER, recv_list(i), MPI_MY_TAG, &
            MPI_COMM_WORLD,ierror)
    end do !i

    ! Complete the receives of the number of externals
    do i = 1, num_recv_neighbors
        call MPI_WAIT(request(i), mystatus, ierror)
        A%send_length(i) = lengths(i)
    end do !i

    deallocate(lengths)

    allocate(A%elements_to_send(A%total_to_be_sent))
    A%elements_to_send = 0

    ! Build "elements_to_send" list.  These are the x elements I own
    ! that need to be sent to other processors.
    MPI_MY_TAG = MPI_MY_TAG + 1
    j = 1

    do i = 1, num_recv_neighbors
        call MPI_IRECV(A%elements_to_send(j), A%send_length(i), MPI_INTEGER, &
            A%neighbors(i),MPI_MY_TAG, MPI_COMM_WORLD, request(i), ierror)
        j = j + A%send_length(i)
    end do !i

    j = 1
    do i = 1, num_recv_neighbors
        start  = j
        newlength = 0

        ! Go through list of external elements
        ! until updating processor changes. This is redundant, but
        ! saves us from recording this information
        do while((j <= num_external) .and. &
                (new_external_processor(j) == (recv_list(i)+1)))
            newlength = newlength + 1
            j = j + 1
            if (j == num_external) exit
        end do

        call MPI_SEND(new_external(start), j-start+1, MPI_INTEGER, recv_list(i), &
            MPI_MY_TAG, MPI_COMM_WORLD, ierror)
    end do

    ! receive from each neighbor the global index list of external elements
    do i = 1, num_recv_neighbors
        call MPI_WAIT(request(i), mystatus, ierror)
    end do !i

    ! replace global indices by local indices
    do i = 1, A%total_to_be_sent
        A%elements_to_send(i) = A%elements_to_send(i) - A%start_row + 1
    end do

    !! Finish up
    A%num_send_neighbors = num_send_neighbors
    A%local_ncol = A%local_nrow + num_external

    ! Used in exchange_externals
    allocate(A%send_buffer(A%total_to_be_sent))

    deallocate(externals)
    deallocate(tmp_buffer)
    deallocate(global_index_offsets)
    deallocate(recv_list)
    deallocate(external_processor)
    deallocate(new_external)
    deallocate(new_external_processor)
    deallocate(request)
end subroutine make_local_matrix


subroutine exchange_externals(A,x)
    implicit none
    type(SparseMatrixT), intent(inout) :: A
    real(kind=dp), allocatable, intent(inout) :: x(:)

    integer :: mysize, myrank
    integer :: i
    integer :: mycount
    integer :: offset
    integer :: MPI_MY_TAG
    integer, allocatable :: request(:)
    integer :: mystatus(MPI_STATUS_SIZE)

    if (single == 1) then
        return
    end if

    call MPI_COMM_SIZE(MPI_COMM_WORLD,mysize,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)

    offset = 0
    MPI_MY_TAG = 99

    allocate(request(A%num_send_neighbors))

    do i = 1, A%num_send_neighbors
        call MPI_IRECV(x(A%local_nrow+1+offset), A%recv_length(i), &
            MPI_DOUBLE_PRECISION, A%neighbors(i), &
            MPI_MY_TAG, MPI_COMM_WORLD, request(i),ierror)
        offset = offset + A%recv_length(i)
    end do !i

    do i = 1, A%total_to_be_sent
        A%send_buffer(i) = x(A%elements_to_send(i))
    end do !i

    offset = 0

    do i = 1, A%num_send_neighbors
        call MPI_SEND(A%send_buffer(1+offset), A%send_length(i), &
            MPI_DOUBLE_PRECISION, A%neighbors(i), &
            MPI_MY_TAG, MPI_COMM_WORLD, ierror)
        offset = offset + A%send_length(i)
    end do !i

    do i = 1, A%num_send_neighbors
        call MPI_WAIT(request(i), mystatus, ierror)
    end do

    deallocate(request)
end subroutine exchange_externals
#endif

end module SparseMatrix

