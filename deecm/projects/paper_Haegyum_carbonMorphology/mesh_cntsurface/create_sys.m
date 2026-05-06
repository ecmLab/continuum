%% create_sys.m
%%% Define the model system. lz_in and load_in let mdl_run.m tune packing.
%%% By: Howard Tu, 06/10/2023 (mod. DRX/C terminology + tunable packing)
function [sys,drx,carbon] = create_sys(drx, carbon, ii, lz_in, load_in)

if nargin < 4 || isempty(lz_in),   lz_in   = 30;   end
if nargin < 5 || isempty(load_in), load_in = 10.0; end

%% 1. System parameters
sys         = struct;
sys.masRt   = (6.5-ii*0.5)/(88.5+ii*0.5);                  % Carbon:DRX mass ratio (~5 wt% binder assumed)
sys.load    = load_in;                                     % Areal loading [mg/cm^2]
sys.lx      = 30;
sys.ly      = 30;
sys.lz      = lz_in;
sys.prosty  = 1 - sys.load/(sys.lz*(1+sys.masRt))*(1/drx.den + sys.masRt/carbon.den);
sys.vol     = sys.lx * sys.ly * sys.lz;
sys.volRt   = sys.masRt * drx.den/carbon.den;
carbon.vol  = (1 - sys.prosty) * sys.vol * sys.volRt/(1 + sys.volRt);
drx.vol     = (1 - sys.prosty) * sys.vol * 1/(1 + sys.volRt);

if sys.prosty <= 0
    error('create_sys: requested loading/box size yields non-positive porosity (%.3f).', sys.prosty);
end

%% 2. Particle counts and volumes
drx = cmp_nmbr(drx, drx.vol);

sphere_volume    = (pi/6) * carbon.dia^3;
cluster_volume   = carbon.cluster.cluster_volume;
carbon.nClusters = max(1, round(carbon.vol / cluster_volume));
carbon.nTot      = carbon.nClusters * carbon.cluster.sphere_count;
carbon.rVol      = carbon.nTot * sphere_volume;
carbon.nDia      = carbon.nTot;

%% 3. Other
sys.mass = drx.rVol*drx.den + carbon.rVol*carbon.den;
end