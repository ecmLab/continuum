#!/bin/bash
visco='1e-3 1e-4'
visco=1e4
current='0.1 0.5 1.0 2.0 3.0'
cap=3.2
porosity='0.4 0.5 0.6 0.7 0.8'
# current='0.1 1.0'
# porosity='0.4'
file_base='unsaturated'
viscosity="Modules/FluidProperties/the_simple_fluid/viscosity=${visco}"
for c in $current
do 
    time=`echo "scale=1;$cap/$c*3600" | bc -l`
    time2=`echo "scale=1;$time*2.0" | bc -l`
    # echo $time
    flux_function="'if (t <= $time,-66.41344e-9*$c*10, 66.41344e-9*$c*10)'"
    ff="BCs/constant_injection_flux/flux_function=$flux_function"
    end_time="Executioner/end_time=$time2"
    dtmax1=`echo "scale=1; 300.0-(100.0*($c-0.1))"| bc -l`
    dtmax2="'$time'"
    sync="Outputs/sync_times/$dtmax2"
    # dtstart=`"echo scale=1; 10.0*300.0-(100.0*($c-0.1))"| bc -l`
    echo "Dt start = $dtstart"
    dtmax="Executioner/dtmax=$dtmax1"
    for p in $porosity
    do 
        poro_state="Materials/porosity_Agc/porosity=$p"
        echo $poro_state
        filebase="unsaturated_current_${c}_porosity_${p}_visco_${visco}"
        fb2="'$filebase'"
        fb="Outputs/file_base=$fb2"
        run_cmd="mpirun -np 4 ~/github/electro_chemo_mech2/electro_chemo_mech2-opt -i test_unsaturated_real_units.i $ff $poro_state $dtmax $end_time $dtmax2 $fb $viscosity"
        echo $run_cmd
        eval $run_cmd
    done
done
