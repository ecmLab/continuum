###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
MOOSE_DIR ?= /home/josva/sim_local/software/moose
###############################################################################
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/../../
endif
# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
ADDITIONAL_CPPFLAGS += -Wall -Wextra
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk
################################## MODULES ####################################
ALL_MODULES               := no
CHEMICAL_REACTIONS        := no
# kept on as a safety net; CZM itself lives in solid_mechanics
CONTACT                   := yes
FLUID_PROPERTIES          := no
HEAT_TRANSFER             := no
LINEAR_ELASTICITY         := no
MISC                      := no
NAVIER_STOKES             := no
PHASE_FIELD               := no
POROUS_FLOW               := no
# REQUIRED: provides BiLinearMixedModeTraction, CZMInterfaceKernel,
# the CohesiveZoneMaster action, plus all of the bulk plasticity /
# eigenstrain machinery used previously.
SOLID_MECHANICS           := yes
XFEM                      := no
include $(MOOSE_DIR)/modules/modules.mk
###############################################################################
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := contact_czm
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk
###############################################################################