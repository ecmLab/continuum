#!/bin/bash
# ---- k = 1e-2, Rct = 1 ohm-cm^2 ---- 
if [ -f run.dat ]
then 
    rm -f run.dat
    touch run.dat
fi
for c in 3e-3 5e-3 10e-3 15e-3 20e-3
do 
    i=`awk -v c=$c 'BEGIN {ans=c*1000; printf"%.1f\n",ans}'`; 
    for rct in 100 10 100 200 500 750 1000
    do 
        k=`awk -v r=$rct 'BEGIN {ans=(1.0/r)*1e-2; printf"%.1e\n",ans}'`
        # echo mpirun -np 4 ~/github/electro_chemo_mech2/electro_chemo_mech2-opt -i electric_diffusion_full_cell_mesh.i Constraints/anode_constraint/k=$k BCs/current/$c Outputs/csv/file_base=csv/rct_${rct}_i_${i} Outputs/out/file_base=rst/rct_${rct}_i_${i} 
        # mpirun -np 4 ~/github/electro_chemo_mech2/electro_chemo_mech2-opt -i electric_diffusion_full_cell_mesh.i Constraints/anode_constraint/k=$k BCs/current/$c Outputs/csv/file_base=csv/rct_${rct}_i_${i} Outputs/out/file_base=rst/rct_${rct}_i_${i}
        # echo mpirun -np 4 ~/github/electro_chemo_mech2/electro_chemo_mech2-opt -i electric_diffusion_full_cell_mesh.i k_anode=${k} current_density=$c Outputs/csv/file_base=csv/rct_${rct}_i_${i} Outputs/out/file_base=rst/rct_${rct}_i_${i} >> run.dat
        # mpirun -np 4 ~/github/electro_chemo_mech2/electro_chemo_mech2-opt -i electric_diffusion_full_cell_mesh.i k_anode=${k} current_density=$c Outputs/csv/file_base=csv/rct_${rct}_i_${i} Outputs/out/file_base=rst/rct_${rct}_i_${i} 
        cat full_cell_parametric.i | sed "s/cddx/$c/g" | sed "s/kax/$k/g" | sed "s/testxxx/rct_${rct}_i_${i}/g" > full_cell_parametric_rct_${rct}_i_${i}.i
        mpirun -np 8 ~/github/electro_chemo_mech2/electro_chemo_mech2-opt -i full_cell_parametric_rct_${rct}_i_${i}.i
    done
done
