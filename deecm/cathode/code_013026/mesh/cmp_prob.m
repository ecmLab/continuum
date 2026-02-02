%% This script is designed for generating models for distribution of both LPS and NCM particles
%%% Function: cmp_prob: Compute the particle probability of all diameters of two material
%%% By: Howard Tu
function ptc = cmp_prob(ptc,nSmp, flg)

% The sampling intervals of current particle
ptc.smp        = linspace(ptc.min,ptc.max,nSmp+1);   

% The overal probability in each sampling bar of current particle
for i = 1 : nSmp
    ptc.prb(i) = integral(@(size)particleDistr(size,ptc.ave,ptc.sgm,flg),ptc.smp(i),ptc.smp(i+1));  
end

% The sampling diameters of current particles
ptc.dia        = ptc.smp(1:end-1) + (ptc.max-ptc.min)/(2*nSmp); 
