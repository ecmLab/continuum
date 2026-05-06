%% This script generates models of the DRX cathode + Carbon for Haegyum's project
%  CNT morphology variant: rods lie tangentially on DRX particle surfaces
%%% By: Howard Tu, 03/11/2026
clc; clear;
% "Unit = micro" unit system: Mass: pg, Distance: um, Time: us, Density: pg/um^3
nWt = 1;

for iwt = 1:nWt

%% Define variables
[drx,carbon,gst] = particle_info();
[sys,drx,carbon]  = create_sys(drx, carbon, iwt);
sys.cord = zeros(0,7);   % 0-row, 7-column so column indexing works in insrtPtc

% Generate output folder
fid  = strcat('massratio/mr',num2str(iwt),'/');
stmp = ['mkdir ' fid];
eval(stmp);

fprintf('The total particles are DRX = %6.4f , Carbon spheres = %6.4f\n',drx.nTot, carbon.nTot);
fprintf('Carbon clusters planned: %d (spheres/cluster=%d)\n', carbon.nClusters, carbon.cluster.sphere_count);
fprintf('The initial porosity is : %6.4f\n',1-(drx.rVol+carbon.rVol)/sys.vol);

% 1. Insert all DRX particles first (surface attachment requires DRX to exist)
sys.cord = insrtPtc(sys, drx, drx.typ, drx.den, 1);

% 2. Attach CNT rods to DRX surfaces (center of rod touches surface, tangentially)
[sys, carbon_stats] = insertCarbonOnDRXSurface(sys, carbon);

% 3. Generate ghost particle layer at top of box (for pressure application)
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

% 4. Output all particles in LAMMPS data file format
stmp   = strcat(fid,'mdl.data');
fileID = fopen(stmp,'w');
fprintf(fileID,'LAMMPS data file\n');
fprintf(fileID,'\n');
fprintf(fileID,'%6d %6s\n',size(sys.cord,1)+size(gst.cord,1),'atoms');
fprintf(fileID,'\n');
fprintf(fileID,'   3 atom types\n');
fprintf(fileID,'\n');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.lx, 'xlo', 'xhi');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.ly, 'ylo', 'yhi');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.lz+gst.dia, 'zlo', 'zhi');
fprintf(fileID,'\n');
fprintf(fileID,'Atoms # sphere\n');
fprintf(fileID,'\n');
fprintf(fileID,'%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n',sys.cord');
fprintf(fileID,'%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n',gst.cord');
fclose(fileID);

fprintf('Carbon cluster placement: inserted %d, failed %d\n', ...
    carbon_stats.inserted_clusters, carbon_stats.failed_clusters);
end


%% ===== Helper functions for surface-attached CNT clusters =====

function [sys, stats] = insertCarbonOnDRXSurface(sys, carbon)
% Attach CNT rods to randomly chosen DRX particle surfaces.
% The center bead of each rod is placed at the DRX surface, oriented tangentially.
    stats = struct('name', 'CNT_surface', 'inserted_clusters', 0, 'failed_clusters', 0);
    if carbon.nClusters <= 0
        return;
    end
    % Extract DRX particles (type 1) as attachment targets
    drx_particles = sys.cord(sys.cord(:,2) == 1, :);
    nDRX = size(drx_particles, 1);
    if nDRX == 0
        warning('No DRX particles found. Cannot attach CNTs to surfaces.');
        return;
    end
    fprintf('Attaching %d CNT clusters to DRX surfaces (tangential)...\n', carbon.nClusters);
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
% Place a CNT rod bent to follow the DRX particle surface (great circle arc).
%
% All N beads are placed at distance R_arc = drx_r + c_r from the DRX center,
% i.e. every bead touches the DRX surface simultaneously.
%
% Arc geometry: bead_k = drx_center + R_arc*(cos(theta_k)*n + sin(theta_k)*v_t)
%   theta_k = k * d_theta,  d_theta = spacing / R_arc
%   k is centered at 0, so the arc is symmetric about the attachment normal n.
%
% Overlap proof: |bead_k - drx_center| = R_arc * sqrt(cos^2 + sin^2) = R_arc
%   = drx_r + c_r = min_dist  for all k  =>  no overlap with attached DRX.
    placed   = false;
    newAtoms = [];
    max_attempts = 5000;
    nDRX    = size(drx_particles, 1);
    c_r     = carbon.dia / 2;
    nSeg    = carbon.cluster.sphere_count;  % number of beads
    spacing = carbon.dia * 0.6;            % arc-length between consecutive bead centers

    for attempt = 1:max_attempts
        % 1. Pick a random DRX particle as attachment host
        idx_drx    = randi(nDRX);
        drx_center = drx_particles(idx_drx, 5:7);
        drx_r      = drx_particles(idx_drx, 3) / 2;

        % 2. Arc radius: all beads sit on the DRX surface
        R_arc = drx_r + c_r;

        % 3. Angular step so arc-length between adjacent beads equals spacing
        d_theta = spacing / R_arc;

        % 4. Angular positions for all beads, symmetric about theta=0
        k      = (0:nSeg-1) - (nSeg-1)/2;  % [1 x nSeg], centered at 0
        thetas = k * d_theta;               % [1 x nSeg]

        % 5. Random surface normal (attachment point) and tangential arc direction
        n   = randn(1, 3);  n   = n   / norm(n);
        v_t = randomTangentialDir(n);

        % 6. Bead positions on the arc (all on the DRX surface)
        %    positions [nSeg x 3] = drx_center + R_arc*(cos*n + sin*v_t)
        positions = drx_center + R_arc * (cos(thetas') * n + sin(thetas') * v_t);

        % 7. Bounds check: all beads inside the simulation box
        if any(positions(:,1) < c_r) || any(positions(:,1) > sys.lx - c_r) || ...
           any(positions(:,2) < c_r) || any(positions(:,2) > sys.ly - c_r) || ...
           any(positions(:,3) < c_r) || any(positions(:,3) > sys.lz - c_r)
            continue;
        end

        % 8. Overlap check against all existing particles
        %    Every bead is at exactly R_arc = min_dist from attached DRX -> passes
        if canPlace(sys, positions, carbon.dia)
            placed   = true;
            newAtoms = buildAtomBlock(sys, positions, carbon);
            break;
        end
    end
end

function v_t = randomTangentialDir(n)
% Generate a uniformly random unit vector in the plane perpendicular to n.
    v   = randn(1, 3);
    v_t = v - dot(v, n) * n;   % remove normal component
    v_t = v_t / norm(v_t);
end

function tf = canPlace(sys, positions, dia)
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
        if any(dist_sq < min_dist_sq * 0.999)
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
