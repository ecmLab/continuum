%% This script is designed for generating models for distribution of both LPS and NCM particles
%%% Function: cmp_prob: Compute the particle probability of all diameters of two material
%%% By: Howard Tu
function ptc = cmp_nmbr(ptc,vol)
                                     
ptc.nTot      = round(6/pi*vol/(ptc.prb * (ptc.dia.^3)'));                  % Approximate total number of particles
ptc.nDia      = round(ptc.nTot*ptc.prb);                                    % Particle number of each diameter
ptc.nTot      = sum(ptc.nDia);                                              % Real total number of particles
ptc.dVol      = pi/6 * ptc.dia.^3 .* ptc.nDia;                              % The volume of each particle size
ptc.rVol      = sum(ptc.dVol);                                              % Real volume of total particles
ptc.dVol      = ptc.dVol/ptc.rVol;

% For particle size distribution plot purpose only
ptc.pNd       = ptc.nTot*ptc.prb;
ptc.pVol      = pi/6 * ptc.dia.^3 .* ptc.pNd;
ptc.pVol      = ptc.pVol/sum(ptc.pVol);