# =======================================================================================
#
#      Filename:  Makefile
#      Description:  Central Makefile
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

include config.mk

export COMPILER = $(strip $(TOOLCHAIN))

#CONFIGURE BUILD SYSTEM
TARGET = HPCCG-$(COMPILER)
BUILD_DIR  = ./$(COMPILER)
SRC_DIR    = ./src
MAKE_DIR   = ./make
Q         ?= @

#DO NOT EDIT BELOW
include config.mk
include $(MAKE_DIR)/include_$(COMPILER).mk
INCLUDES  += -I$(BUILD_DIR)/ -I./includes

VPATH     = $(SRC_DIR)
SRC       = $(wildcard $(SRC_DIR)/*.F90)
OBJ       = $(patsubst $(SRC_DIR)/%.F90, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.F90))
OBJ       += $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c))
CPPFLAGS := $(INCLUDES) $(DEFINES) $(DEFINES)

ifeq ($(strip $(ENABLE_MPI)),YES)
FC = $(MPIFC)
FCFLAGS += $(MPIFLAGS)
DEFINES += -D_MPI
endif

ifeq ($(strip $(ENABLE_OMP)),YES)
FCFLAGS := $(OMPFLAGS) $(FCFLAGS)
DEFINES += -D_OMP
endif

ifeq ($(strip $(HAVE_COOP)),YES)
DEFINES += -DHAVE_COOP
endif

DEFINES += -DTIMINGLEV=$(strip $(TIMER_LEVEL))

$(TARGET): $(BUILD_DIR) $(OBJ)
	$(info ===>  LINKING  $(TARGET))
	$(Q)$(FC)  $(CPPFLAGS) $(FCFLAGS) -o $@ $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o:  %.F90
	$(info ===>  COMPILE  $@)
	$(Q)$(FC) -c  $(CPPFLAGS)  $(FCFLAGS) $< -o $@

$(BUILD_DIR)/%.o:  %.c
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c  $(CPPFLAGS)  $(CFLAGS) $< -o $@

.PHONY: clean distclean tags info asm

clean:
	$(info ===>  CLEAN)
	@rm -rf $(BUILD_DIR)
	@rm -f tags

distclean: clean
	$(info ===>  DIST CLEAN)
	@rm -f $(TARGET)

info:
	$(info $(FCFLAGS))
	$(Q)$(FC) $(VERSION)

asm:  $(BUILD_DIR) $(ASM)

tags:
	$(info ===>  GENERATE TAGS)
	$(Q)ctags -R

$(BUILD_DIR):
	@mkdir $(BUILD_DIR)

-include $(OBJ:.o=.d)

# Add explicit dependencies for Fortran90 modules
include $(MAKE_DIR)/dep.mk
