#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=12
#SBATCH --array=1-63
#SBATCH --job-name=MOOSE_ParamSweep
#SBATCH --output=slurm_output/slurm-%A_%a.out
#SBATCH --error=slurm_output/slurm-%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cunezben@psc.edu
#SBATCH -A eve230003p

# type 'man sbatch' for more information and options
# this job will run 63 array tasks (parameter combinations)
# Sweep 1: 9 Young's moduli × 4 stack pressures = 36 combinations  
# Sweep 2: 9 Young's moduli × 3 yield strengths = 27 combinations
# Total: 63 parameter combinations
# each task uses 12 cores on RM-shared partition for 24 Hours
# this job would potentially charge 288 RM SUs per task

# echo commands to stdout
set -x

# move to working directory
# this job assumes:
# - all input data is stored in this directory
# - all output should be stored in this directory
# - please note that eve230003p is the groupname/allocation
# - cunezben is the PSC username
cd /ocean/projects/eve230003p/cunezben/projects/moose/Joseph/NMC_LPSCl_FEA/BC_Free_Top_Right/

# Create output directories
mkdir -p slurm_output
mkdir -p rst
mkdir -p results_sweep1
mkdir -p results_sweep2

# Load required modules for Bridges-2
module purge
module load gcc/10.2.0
module load openmpi/4.0.5-gcc10.2.0
module load anaconda3/2022.10

# Set MPI environment variables
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export F77=mpif77

# Define parameter arrays
# Young's moduli from 200 to 1000 MPa in steps of 50 MPa (17 values)
ymods=(200 300 400 500 600 700 800 900 1000)

# Stack pressures for first sweep (4 values)
stack_pressures_1=(0 2 4 6)

# Yield strengths for second sweep (3 values)  
yield_strengths_2=(6 96 196)

# Calculate which parameter combination this array task should run
task_id=$SLURM_ARRAY_TASK_ID

if [ $task_id -le 36 ]; then
    # First sweep: ymod vs stack pressure (ystr=196)
    # Tasks 1-36: 9 ymods × 4 stack_pressures
    ymod_index=$(( ($task_id - 1) / 4 ))
    sp_index=$(( ($task_id - 1) % 4 ))
    
    ymod=${ymods[$ymod_index]}
    sptop=${stack_pressures_1[$sp_index]}
    ystr=196
    
    echo "Sweep 1 - Task $task_id: Young's Modulus=$ymod MPa, Stack Pressure=$sptop MPa, Yield Strength=$ystr MPa"
    echo "Start time: $(date)"
    
    # Create results directory if it doesn't exist
    results_dir="results_sweep1"
    mkdir -p $results_dir
    mkdir -p $results_dir/rst  # Create rst subdirectory for MOOSE output
    
    # Check if input files exist
    if [ ! -f "AG_Contact_Loss.i" ]; then
        echo "ERROR: Input file AG_Contact_Loss.i not found!"
        exit 1
    fi
    
    if [ ! -f "input_mesh_Shafee.msh" ]; then
        echo "ERROR: Mesh file input_mesh_Shafee.msh not found!"
        exit 1
    fi
    
    # run pre-compiled program which is already in your project space
mpirun -np 32 ./my_solid_app-opt \
        -i AG_Contact_Loss.i \
        ymod=$ymod \
        ystr=$ystr \
        sptop=$sptop \
        > $results_dir/log_task${SLURM_ARRAY_TASK_ID}_ymod${ymod}_sp${sptop}_ystr${ystr}.out 2>&1
        
else
    # Second sweep: ymod vs yield strength (sptop=6)
    # Tasks 37-63: 9 ymods × 3 yield_strengths  
    adjusted_task=$(( $task_id - 36 ))
    ymod_index=$(( ($adjusted_task - 1) / 3 ))
    ystr_index=$(( ($adjusted_task - 1) % 3 ))
    
    ymod=${ymods[$ymod_index]}
    ystr=${yield_strengths_2[$ystr_index]}
    sptop=6
    
    echo "Sweep 2 - Task $task_id: Young's Modulus=$ymod MPa, Stack Pressure=$sptop MPa, Yield Strength=$ystr MPa"
    
    # Create results directory if it doesn't exist
    results_dir="results_sweep2"
    mkdir -p $results_dir
    mkdir -p $results_dir/rst  # Create rst subdirectory for MOOSE output
    
    # Check if input files exist
    if [ ! -f "AG_Contact_Loss.i" ]; then
        echo "ERROR: Input file AG_Contact_Loss.i not found!"
        exit 1
    fi
    
    if [ ! -f "input_mesh_Shafee.msh" ]; then
        echo "ERROR: Mesh file input_mesh_Shafee.msh not found!"
        exit 1
    fi
    
    # run pre-compiled program which is already in your project space
mpirun -np 32 ./my_solid_app-opt \
        -i AG_Contact_Loss.i \
        ymod=$ymod \
        ystr=$ystr \
        sptop=$sptop \
        > $results_dir/log_task${SLURM_ARRAY_TASK_ID}_ymod${ymod}_sp${sptop}_ystr${ystr}.out 2>&1
fi

# Check exit status
if [ $? -eq 0 ]; then
    echo "Task $SLURM_ARRAY_TASK_ID completed successfully"
    echo "End time: $(date)"
else
    echo "Task $SLURM_ARRAY_TASK_ID failed with exit code $?"
    echo "End time: $(date)"
fi
