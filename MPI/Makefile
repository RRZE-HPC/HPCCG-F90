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

SRC_DIR    = ./src
MAKE_DIR   = ./make

#DO NOT EDIT BELOW

TARGET = HPCCG-$(COMPILER)

include $(MAKE_DIR)/include_$(COMPILER).mk

BUILD_DIR  = ./$(COMPILER)
Q         ?= @
VPATH     = $(SRC_DIR)
SRC       = $(wildcard $(SRC_DIR)/*.F90)
OBJ       = $(patsubst $(SRC_DIR)/%.F90, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.F90))
OBJ       += $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c))
INCLUDES  += -I$(BUILD_DIR)/ -I./includes

LIBS :=  $(LIBS)
CPPFLAGS += $(INCLUDES) $(DEFINES)

ifeq ($(strip $(ENABLE_MPI)),YES)
FC = $(MPIFC)
FCFLAGS += $(MPIFLAGS)
DEFINES += -D_MPI
endif

ifeq ($(strip $(ENABLE_OMP)),YES)
FCFLAGS := $(OMPFLAGS) $(FCFLAGS)
DEFINES += -D_OMP
endif

DEFINES += -DTIMINGLEV=$(strip $(TIMER_LEVEL))

all:  $(TARGET)

$(TARGET): $(BUILD_DIR) $(OBJ)
	@echo "===>  LINK  $@"
	$(Q)$(FC)  $(CPPFLAGS) $(FCFLAGS) -o $@ $(OBJ) $(LIBS)

$(BUILD_DIR):
	@mkdir $(BUILD_DIR)

$(BUILD_DIR)/%.o:  %.F90
	@echo "===>  COMPILE  $@"
	$(Q)$(FC) -c  $(CPPFLAGS)  $(FCFLAGS) $< -o $@

$(BUILD_DIR)/%.o:  %.c
	@echo "===>  COMPILE  $@"
	$(Q)$(CC) -c  $(CPPFLAGS)  $(CFLAGS) $< -o $@

.PHONY: clean distclean

clean:
	@echo "===>  CLEAN"
	@rm -rf $(BUILD_DIR)

distclean: clean
	@echo "===>  DIST CLEAN"
	@rm -f $(TARGET)


include  $(MAKE_DIR)/dep.mk

