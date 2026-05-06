%% particle_info.m
%%% Function: define geometry and material properties of all particles
%%% By: Howard Tu, 06/10/2023 (mod. for DRX/C terminology, L/D parameter)
function [drx,carbon,gst] = particle_info(LD_ratio)

if nargin < 1 || isempty(LD_ratio)
    LD_ratio = 10;          % Default CNT slenderness ratio
end

%% 1. DRX particles
drx         = struct;
drx.nD      = 40;
drx.den     = 4.0;          % g/cm^3 == pg/μm^3
drx.typ     = 1;
drx.ave     = 1;
drx.min     = 0.5;
drx.max     = 2;
drx.sgm     = 0.2;

%% 2. Conductive carbon (CNT) particles
carbon            = struct;
carbon.den        = 2.2;
carbon.typ        = 2;
carbon.dia        = 0.1;
carbon.ave        = carbon.dia;
carbon.LD_ratio   = LD_ratio;
carbon.length     = LD_ratio * carbon.dia;   % CNT rod length (μm)

% Default intra-cluster spacing (overridable in mdl_run.m). Fraction of
% carbon.dia for center-to-center distance between adjacent beads in a rod.
carbon.intra_cc   = 0.85;

carbon.cluster    = buildCNTCluster(carbon.dia, carbon.length, carbon.intra_cc);

%% 3. Ghost particles
gst         = struct;
gst.den     = 4.5;
gst.typ     = 3;

%% 4. Theoretical particle probabilities
drx = cmp_prob(drx, drx.nD, 1);

end

function cluster = buildCNTCluster(base_dia, length_um, intra_cc)
    if length_um <= 0
        length_um = base_dia;
    end
    spacing = base_dia * intra_cc;
    nSeg = max(2, ceil(length_um / spacing) + 1);
    z_positions = linspace(-length_um/2, length_um/2, nSeg)';
    offsets = [zeros(nSeg,2), z_positions];
    cluster.name           = 'CNT';
    cluster.offsets        = offsets;
    cluster.sphere_count   = size(offsets,1);
    cluster.allow_rotation = true;
    cluster.sphere_volume  = (pi/6) * base_dia^3;
    cluster.cluster_volume = cluster.sphere_count * cluster.sphere_volume;
end