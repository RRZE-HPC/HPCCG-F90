! =======================================================================================
!
!      Filename:  timer.F90
!      Description:  Implements timing and likwid profiling regions
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

module timer
#ifdef PERFCTR
    use likwid
#endif
    implicit none

    type interval
        double precision startt
        double precision stopt
    end type interval

    logical ::  master
    integer :: ierror

    private :: master, ierror
contains

    subroutine timer_init(  )
        if (this_image() == 1) then
            master = .true.
        else
            master = .false.
        endif
#ifdef PERFCTR
        call likwid_markerInit()
#endif
    end subroutine timer_init

    subroutine timer_start( timing, region_name )
        implicit none
        type( interval ), intent( inout ) :: timing
        character (len=*), intent( in ) :: region_name
#ifdef PERFCTR
        call likwid_markerStartRegion(region_name)
#else
        sync all
        if(master) then
        call timingc(timing%startt)
        end if
#endif
    end subroutine timer_start

    subroutine timer_stop( timing, region_name )
        implicit none
        type( interval ), intent( inout ) :: timing
        character (len=*), intent( in ) :: region_name
#ifdef PERFCTR
        call likwid_markerstopRegion(region_name)
#else
        sync all
        if(master) then
        call timingc(timing%stopt)
        endif
#endif
    end subroutine timer_stop

    subroutine timer_output( timing , region_name )
        implicit none
        type( interval ), intent( in ) :: timing
        character (len=*), intent( in ) :: region_name
        real :: timing_result
#ifdef PERFCTR
        timing_result = 0
#else
        if(master) then
            timing_result = timing%stopT - timing%startT
            print *,"TIMER: ",region_name, timing_result , "sec"
        endif
#endif
    end subroutine timer_output

    function timer_print( timing ) result (timing_result)
        implicit none
        type( interval ), intent( in ) :: timing
        real :: timing_result
#ifdef PERFCTR
        timing_result = 0
#else
        timing_result = timing%stopT - timing%startT
#endif
    end function timer_print

    subroutine timer_finalize()
#ifdef PERFCTR
        call likwid_markerClose()
#endif
    end subroutine timer_finalize
end module timer

