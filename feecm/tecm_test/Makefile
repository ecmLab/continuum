###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#  Copyright 2025, ToBeDecided, All Rights Reserved
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################
## Use the MOOSE submodule if it exists and MOOSE_DIR is not set
# MOOSE_SUBMODULE    := $(CURDIR)/moose
# ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
#   MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
# else
#   MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
# endif

MOOSE_DIR := /Users/howardtu/Documents/modeling/localPackage/moose

#ifeq ($(wildcard $(MOOSE_DIR)/framework/Makefile),)
# $(error MOOSE directory not found at $(MOOSE_DIR))
#endif

## framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# To use certain physics included with MOOSE, set variables below to
# yes as needed.  Or set ALL_MODULES to yes to turn on everything (overrides
# other set variables).

ALL_MODULES                 := no

CHEMICAL_REACTIONS          := yes
CONTACT                     := no
ELECTROMAGNETICS            := no
EXTERNAL_PETSC_SOLVER       := no
FLUID_PROPERTIES            := no
FSI                         := no
FUNCTIONAL_EXPANSION_TOOLS  := no
GEOCHEMISTRY                := no
HEAT_TRANSFER               := yes
LEVEL_SET                   := no
MISC                        := no
NAVIER_STOKES               := yes
OPTIMIZATION                := no
PERIDYNAMICS                := no
PHASE_FIELD                 := no
POROUS_FLOW                 := no
RAY_TRACING                 := no
REACTOR                     := no
RDG                         := no
RICHARDS                    := no
SOLID_MECHANICS             := yes
STOCHASTIC_TOOLS            := no
THERMAL_HYDRAULICS          := yes
XFEM                        := no

include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := tecm_test
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
