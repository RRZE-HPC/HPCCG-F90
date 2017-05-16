! =======================================================================================
!
!      Filename:  main.F90
!      Description:  Main program for HPCCG
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

subroutine printInfo(A)
    use SparseMatrix
#ifdef _OMP
    use omp_lib
#endif

    implicit none
    type(SparseMatrixT), intent(in) :: A

    integer :: mysize, myrank
    integer :: num_threads
#ifdef _MPI
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
    if (myrank == 0) then
#else
    myrank = 0
    mysize = 1
#endif
#ifdef _OMP
!$OMP PARALLEL
    num_threads = omp_get_num_threads()
!$OMP END PARALLEL
#endif

    print *,"HPCCG Benchmark"
    print *,"Using double precision"
    print *,"Total size (N x N):",A%total_nrow
    print *,"Total nonzeroes:",A%total_nnz
#ifdef _MPI
    print *,"MPI enabled using",mysize,"processes."
#else
    print *,"MPI disabled."
#endif
#ifdef _OMP
    print *,"OpenMP enabled using",num_threads,"threads."
#else
    print *,"OpenMP disabled."
#endif
#ifdef _MPI
    endif
#endif
end subroutine printInfo

program HPCCGmain
    use constants
    use SparseMatrix
    use numerics
    use timer
#ifdef _MPI
    use mpi
#endif

    implicit none
    character(len=32) :: arg
    integer :: argc
    type(SparseMatrixT) :: A
    real(kind=dp), allocatable :: x(:)
    real(kind=dp), allocatable :: b(:)
    real(kind=dp), allocatable :: xexact(:)
    real(kind=dp) :: norm, d
    integer :: nx, ny, nz
    integer :: niters, max_iter
    real(kind=dp) :: normr, tolerance
    real(kind=dp) :: residual

#ifdef _MPI
    call MPI_INIT(ierror)
#endif

    argc = command_argument_count()

    if (argc == 3) then
        call get_command_argument(1, arg)
        read (arg, '(i)') nx
        call get_command_argument(2, arg)
        read (arg, '(i)') ny
        call get_command_argument(3, arg)
        read (arg, '(i)') nz
    else
        print *, 'USAGE: ./HPCCG <nx> <ny> <nz>'
        stop
    end if

    call timer_init()
    call generate_matrix(nx, ny, nz, A, x, b, xexact)
    call printInfo(A)

#ifdef _MPI
    call make_local_matrix(A)
#endif

    niters = 0
    normr = 0.0
    max_iter = 150
    tolerance = 0.0

    call cg( A, b, x, max_iter, tolerance, niters, normr)

#ifdef _MPI
    call MPI_FINALIZE(ierror)
#endif
end program HPCCGmain

