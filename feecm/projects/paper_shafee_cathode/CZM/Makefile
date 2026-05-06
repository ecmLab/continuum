###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
MOOSE_DIR ?= /home/josva/sim_local/software/moose
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
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
# CONTACT was the strict requirement of the original [Contact] block.
# The CZM version no longer needs it (CohesiveZoneMaster + the CZM
# materials live in solid_mechanics, BreakMeshByBlockGenerator lives in
# the framework).  We keep it ON because it is harmless, it does not
# break anything, and it lets the same executable still run the
# original AG_Contact_Loss.i side-by-side.  Set to 'no' if you want a
# leaner build.
CONTACT                   := yes
FLUID_PROPERTIES          := no
HEAT_TRANSFER             := no
LINEAR_ELASTICITY         := no
MISC                      := no
NAVIER_STOKES             := no
PHASE_FIELD               := no
POROUS_FLOW               := no
# solid_mechanics is required for:
#   - IsotropicPlasticityStressUpdate
#   - ComputeMultipleInelasticStress
#   - ComputeThermalExpansionEigenstrain
#   - ComputeIsotropicElasticityTensor
#   - The TensorMechanics Master action
#   - The CohesiveZoneMaster action and the CZM traction-separation
#     materials it needs (e.g. SalehaniIraniNonCoupledTraction,
#     PureElasticTractionSeparation, BiLinearMixedModeTraction, ...)
SOLID_MECHANICS           := yes
XFEM                      := no
include $(MOOSE_DIR)/modules/modules.mk
###############################################################################
# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := contact_czm
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk
###############################################################################