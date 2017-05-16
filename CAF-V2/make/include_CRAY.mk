# =======================================================================================
#
#      Filename:  include_CRAY.mk
#      Description:  MPI build configuration for Cray tool chain
#      Date: Jun 24, 2014
#      Author:  Jan Eitzinger (je), jan.eitzinger@fau.de
#      Copyright (c) 2017, Jan Eitzinger
#      All rights reserved.
#
#      Redistribution and use in source and binary forms, with or without
#      modification, are permitted provided that the following conditions are met:
#
#      * Redistributions of source code must retain the above copyright notice, this
#        list of conditions and the following disclaimer.
#
#      * Redistributions in binary form must reproduce the above copyright notice,
#        this list of conditions and the following disclaimer in the documentation
#        and/or other materials provided with the distribution.
#
#      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#      ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#      WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#      DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#      FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#      DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#      SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#      CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#      OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# =======================================================================================

FC  = ftn
MPIFC  = ftn
CC  = cc
AR  = ar

OMPFLAGS = -h omp
OPT =  -O3
DEBUGFLAGS =
ARCHFLAGS =
FCFLAGS   =  $(OPT) $(ARCHFLAGS) $(DEBUGFLAGS) -h caf -em  -J ./$(COMPILER)
MPIFLAGS   =
CFLAGS   =  -O2 -h c99
ASFLAGS  = -gdwarf-2
CPPFLAGS =
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     =

