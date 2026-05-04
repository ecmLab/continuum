%% mdl_run.m
%  Generates DRX cathode + CNT models. CNT rods wrap around DRX surfaces
%  in random configurations (geodesic random walks).
%  If the user-specified L/D cannot place enough CNTs to hit the target
%  carbon weight ratio, L/D is decreased by 1 and additional shorter rods
%  are inserted, repeating until the target is met or L/D = min_LD.
%%% By: Howard Tu, 03/11/2026
clc; clear;
% "Unit = micro" unit system: Mass: pg, Distance: um, Time: us, Density: pg/um^3
nWt = 11;

%% ============================================================
%  USER-CONFIGURABLE PARAMETERS
%% ============================================================
LD_ratio       = 100;
load_mg_cm2    = 10.0;
box_lz_um      = 30;
olap_cc        = 0.85;
olap_dc        = 0.85;
intra_cc       = 0.85;
wrap_perturb   = 0.30;
min_LD         = 5;     % floor for adaptive L/D reduction

for iwt = 1:nWt

[drx,carbon,gst] = particle_info(LD_ratio);

carbon.olap_cc      = olap_cc;
carbon.olap_dc      = olap_dc;
carbon.intra_cc     = intra_cc;
carbon.wrap_perturb = wrap_perturb;

carbon.cluster = rebuildCNTCluster(carbon.dia, carbon.length, carbon.intra_cc);

[sys,drx,carbon] = create_sys(drx, carbon, iwt, box_lz_um, load_mg_cm2);
sys.cord = zeros(0,7);

fid  = strcat('massratio/mr',num2str(iwt),'/');
stmp = ['mkdir ' fid];
eval(stmp);

fprintf('CNT initial L/D = %g  ->  L = %.3f μm, %d beads/rod, spacing = %.4f μm\n', ...
    carbon.LD_ratio, carbon.length, carbon.cluster.sphere_count, carbon.dia*carbon.intra_cc);
fprintf('Total particles: DRX = %d, target carbon spheres = %d (%d clusters)\n', ...
    drx.nTot, carbon.nTot, carbon.nClusters);
fprintf('Initial target porosity: %6.4f\n', 1-(drx.rVol+carbon.rVol)/sys.vol);

% 1. Insert DRX particles
sys.cord = insrtPtc(sys, drx, drx.typ, drx.den, 1);

% 2. Iteratively wrap CNT rods, dropping L/D by 1 until target carbon
%    volume is reached (or min_LD hit).
target_carbon_volume = carbon.vol;
sphere_volume        = (pi/6) * carbon.dia^3;
placed_carbon_volume = 0;
total_inserted       = 0;
total_failed         = 0;
current_LD           = LD_ratio;

while placed_carbon_volume < target_carbon_volume && current_LD >= min_LD
    carbon.LD_ratio = current_LD;
    carbon.length   = current_LD * carbon.dia;
    carbon.cluster  = rebuildCNTCluster(carbon.dia, carbon.length, carbon.intra_cc);

    remaining        = target_carbon_volume - placed_carbon_volume;
    nNeeded          = max(1, round(remaining / carbon.cluster.cluster_volume));
    carbon.nClusters = nNeeded;

    fprintf('CNT batch: L/D=%d, L=%.3f μm, %d beads/rod, target %d rods\n', ...
        current_LD, carbon.length, carbon.cluster.sphere_count, nNeeded);

    [sys, batch_stats] = insertCarbonOnDRXSurface(sys, carbon);

    placed_in_batch      = batch_stats.inserted_clusters * carbon.cluster.sphere_count * sphere_volume;
    placed_carbon_volume = placed_carbon_volume + placed_in_batch;
    total_inserted       = total_inserted + batch_stats.inserted_clusters;
    total_failed         = total_failed + batch_stats.failed_clusters;

    fprintf('  Placed %d/%d clusters. Carbon volume %.4f / %.4f μm^3 (%.2f%%)\n', ...
        batch_stats.inserted_clusters, nNeeded, ...
        placed_carbon_volume, target_carbon_volume, ...
        100*placed_carbon_volume/max(target_carbon_volume,eps));

    if placed_carbon_volume < target_carbon_volume
        if current_LD <= min_LD
            fprintf('Warning: target carbon volume not met at minimum L/D=%d. Stopping.\n', min_LD);
            break;
        end
        % If a batch placed nothing, decrement still proceeds via this path.
        current_LD = current_LD - 5;
    end
end

carbon_stats = struct('name','CNT_surface', ...
                      'inserted_clusters',     total_inserted, ...
                      'failed_clusters',       total_failed, ...
                      'final_LD',              current_LD, ...
                      'placed_carbon_volume',  placed_carbon_volume, ...
                      'target_carbon_volume',  target_carbon_volume);

% 3. Ghost particle layer (top of box)
gst.max  = min([carbon.dia, min(drx.dia)]);
gst.nxy  = ceil(3*(sqrt(2)-1) * sys.lx/gst.max);
gst.dia  = sys.lx/gst.nxy;
gst.cord = zeros(gst.nxy^2,7);
gst.cord(:,1) = [size(sys.cord,1)+1 : size(sys.cord,1)+gst.nxy^2];
gst.cord(:,2) = gst.typ;
gst.cord(:,3) = gst.dia;
gst.cord(:,4) = gst.den;
gst.cord(:,7) = sys.lz + gst.dia/2;
for ig = 1 : gst.nxy
    for jg = 1 : gst.nxy
        gst.cord((ig-1)*gst.nxy+jg,5) = gst.dia/2 + gst.dia*(jg-1);
        gst.cord((ig-1)*gst.nxy+jg,6) = gst.dia/2 + gst.dia*(ig-1);
    end
end

% 4. Output LAMMPS data file
stmp   = strcat(fid,'mdl.data');
fileID = fopen(stmp,'w');
fprintf(fileID,'LAMMPS data file\n\n');
fprintf(fileID,'%6d %6s\n',size(sys.cord,1)+size(gst.cord,1),'atoms');
fprintf(fileID,'\n   3 atom types\n\n');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.lx, 'xlo', 'xhi');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.ly, 'ylo', 'yhi');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.lz+gst.dia, 'zlo', 'zhi');
fprintf(fileID,'\nAtoms # sphere\n\n');
fprintf(fileID,'%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n',sys.cord');
fprintf(fileID,'%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n',gst.cord');
fclose(fileID);

fprintf('Final: clusters placed = %d, failed = %d, final L/D = %d\n', ...
    carbon_stats.inserted_clusters, carbon_stats.failed_clusters, carbon_stats.final_LD);
fprintf('Carbon volume realized: %.4f / %.4f μm^3 (%.2f%%)\n\n', ...
    placed_carbon_volume, target_carbon_volume, ...
    100*placed_carbon_volume/max(target_carbon_volume,eps));
end


%% ===== CNT cluster rebuilder (consistent with particle_info) =====
function cluster = rebuildCNTCluster(base_dia, length_um, intra_cc)
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


%% ===== Surface-attached CNT placement =====

function [sys, stats] = insertCarbonOnDRXSurface(sys, carbon)
    stats = struct('name', 'CNT_surface', 'inserted_clusters', 0, 'failed_clusters', 0);
    if carbon.nClusters <= 0
        return;
    end
    drx_particles = sys.cord(sys.cord(:,2) == 1, :);
    if isempty(drx_particles)
        error('No DRX particles found. Cannot attach CNTs to surfaces.');
    end
    fprintf('Wrapping %d CNT rods (L/D=%g) around DRX surfaces...\n', ...
        carbon.nClusters, carbon.LD_ratio);
    for c = 1:carbon.nClusters
        [placed, newAtoms] = placeCNTTangentially(sys, carbon, drx_particles);
        if placed
            sys.cord = [sys.cord; newAtoms]; %#ok<AGROW>
            stats.inserted_clusters = stats.inserted_clusters + 1;
        else
            stats.failed_clusters = stats.failed_clusters + 1;
        end
        if mod(c, 100) == 0
            fprintf('  progress %d/%d\n', c, carbon.nClusters);
        end
    end
end

function [placed, newAtoms] = placeCNTTangentially(sys, carbon, drx_particles)
    placed       = false;
    newAtoms     = [];
    max_attempts = 5000;
    nDRX    = size(drx_particles, 1);
    c_r     = carbon.dia / 2;
    nSeg    = carbon.cluster.sphere_count;
    spacing = carbon.dia * carbon.intra_cc;

    for attempt = 1:max_attempts
        idx_drx    = randi(nDRX);
        drx_center = drx_particles(idx_drx, 5:7);
        drx_r      = drx_particles(idx_drx, 3) / 2;

        R_arc   = drx_r + c_r;
        d_theta = spacing / R_arc;

        positions = randomGeodesicWalk(drx_center, R_arc, nSeg, d_theta, carbon.wrap_perturb);

        if any(positions(:,1) < c_r) || any(positions(:,1) > sys.lx - c_r) || ...
           any(positions(:,2) < c_r) || any(positions(:,2) > sys.ly - c_r) || ...
           any(positions(:,3) < c_r) || any(positions(:,3) > sys.lz - c_r)
            continue;
        end

        if hasSelfCollision(positions, spacing)
            continue;
        end

        if canPlace(sys, positions, carbon.dia, carbon.olap_cc, carbon.olap_dc)
            placed   = true;
            newAtoms = buildAtomBlock(sys, positions, carbon);
            break;
        end
    end
end

function positions = randomGeodesicWalk(center, R, nSeg, d_theta, perturb)
    u = randn(1,3); u = u / norm(u);
    t = randn(1,3); t = t - dot(t,u)*u; t = t / norm(t);
    positions = zeros(nSeg, 3);
    for i = 1:nSeg
        positions(i,:) = center + R*u;
        if i == nSeg
            break;
        end
        u_new = u*cos(d_theta) + t*sin(d_theta);
        t_new = -u*sin(d_theta) + t*cos(d_theta);
        u = u_new / norm(u_new);
        t = t_new - dot(t_new, u)*u;
        t = t / norm(t);
        if perturb > 0
            phi = (rand-0.5) * 2 * perturb;
            t = t*cos(phi) + cross(u,t)*sin(phi);
            t = t - dot(t,u)*u;
            t = t / norm(t);
        end
    end
end

function tf = hasSelfCollision(positions, spacing)
    n = size(positions, 1);
    min_sq = (spacing * 0.95)^2;
    tf = false;
    for i = 1:n-2
        d2 = sum((positions(i+2:end,:) - positions(i,:)).^2, 2);
        if any(d2 < min_sq)
            tf = true;
            return;
        end
    end
end

function tf = canPlace(sys, positions, dia, olap_cc, olap_dc)
    tf = true;
    if isempty(sys.cord)
        return;
    end
    for i = 1:size(positions,1)
        pos   = positions(i,:);
        mnBnd = pos - dia;
        mxBnd = pos + dia;
        nearby = sys.cord(sys.cord(:,5) >= mnBnd(1) & sys.cord(:,5) <= mxBnd(1) & ...
                          sys.cord(:,6) >= mnBnd(2) & sys.cord(:,6) <= mxBnd(2) & ...
                          sys.cord(:,7) >= mnBnd(3) & sys.cord(:,7) <= mxBnd(3), :);
        if isempty(nearby)
            continue;
        end
        dist_sq     = (pos(1) - nearby(:,5)).^2 + (pos(2) - nearby(:,6)).^2 + (pos(3) - nearby(:,7)).^2;
        min_dist_sq = ((dia + nearby(:,3))/2).^2;

        tol = ones(size(nearby,1), 1);
        tol(nearby(:,2) == 2) = olap_cc;
        tol(nearby(:,2) == 1) = olap_dc;

        if any(dist_sq < min_dist_sq .* tol)
            tf = false;
            return;
        end
    end
end

function block = buildAtomBlock(sys, positions, carbon)
    n         = size(positions, 1);
    block     = zeros(n, 7);
    start_idx = size(sys.cord, 1);
    for i = 1:n
        block(i,:) = [start_idx + i, carbon.typ, carbon.dia, carbon.den, positions(i,:)];
    end
end