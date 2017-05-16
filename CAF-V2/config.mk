# =======================================================================================
#
#      Filename: config.mk
#
#      Description:  Configuration for build system
#
#      Date: Jun 24, 2014
#
#      Author:  Jan Treibig (jt), jan.treibig@gmail.com
#      Project:  AREVA CAF / MPI study
#
# =======================================================================================

# Toolchain configuration to use. Currently supported: IFORT, CRAY
TOOLCHAIN = IFORT

# Configuration option to control MPI distributed memory parallelism
# Can be YES or NO
ENABLE_MPI = NO

# Configuration option to control OPENMP shared memory parallelism
# Can be YES or NO
ENABLE_OMP = NO

# Configure if CAF implementation already supports collective intrinsics
HAVE_COOP = NO

# Timer Level. Higher means fine grain timing output
TIMER_LEVEL = 1

