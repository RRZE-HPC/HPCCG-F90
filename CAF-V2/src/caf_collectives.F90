! =============================================================================
!
!      Filename:  caf_collectives.F90
!      Description:  CAF fallback implementation for collective
!                    operations required in HPCCG.
!      Date: Aug 20, 2014
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
! =============================================================================

module caf_collectives
    use constants
contains
    subroutine CAF_CO_SUM(val)
        real(kind=dp), intent(inout) :: val

        integer :: i
        real(kind=dp), allocatable :: arr(:)[:]
        allocate(arr(num_images())[*])

        ! Reduction
        arr(this_image())[1] = val
        sync all

        ! Broadcast
        if (this_image() == 1) then
            val = SUM(arr(:))

            do i = 1, num_images()
            arr(1)[i] = val
            end do
        end if
        sync all

        val = arr(1)
        deallocate(arr)
    end subroutine CAF_CO_SUM

    subroutine CAF_CO_MAX(val)
        integer, intent(inout) :: val

        integer :: i
        integer, allocatable :: arr(:)[:]
        allocate(arr(num_images())[*])

        ! Reduction
        arr(this_image())[1] = val
        sync all

        ! Broadcast
        if (this_image() == 1) then
            val = MAXVAL(arr(:)[1])

            do i = 1, num_images()
            arr(1)[i] = val
            end do
        end if
        sync all

        val = arr(1)
        deallocate(arr)
    end subroutine CAF_CO_MAX

end module caf_collectives

