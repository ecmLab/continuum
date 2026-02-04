%% This script is designed for generating models of the DRX cathode + Carbon used for Haegyum's project
%  Different models under different particle distributions of DRX to carbon
%%% By: Howard Tu, 01/30/2026
clc; clear;
% "Unit = micro" unit system is used to be consistent with the LAMMPS input file
% Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3
nWt    = 11;  

for iwt = 1:nWt

%% Define variables
[drx,carbon,gst] = particle_info();
% 2. Generate system, including box and system level variables
[sys,drx,carbon] = create_sys(drx, carbon,iwt);
sys.cord = [];

% 0. Generate folder of different files
    fid = strcat('massratio/mr',num2str(iwt),'/');
    stmp = ['mkdir ' fid];
    eval(stmp);

fprintf('The total particles are DRX = %6.4f , Carbon spheres = %6.4f\n',drx.nTot, carbon.nTot);
fprintf('Carbon clusters planned: %d (spheres/cluster=%d)\n', carbon.nClusters, carbon.cluster.sphere_count);
fprintf('The initial porosity is : %6.4f\n',1-(drx.rVol+carbon.rVol)/sys.vol);

% 3. Insert all Particles into the box, stored in sys.cord matrix:
  % 1st colume: id of particle;
  % 2th column: particle type: Ag or aC;
  % 3th column: diameter of particle;
  % 4th column: density of particle
  % 5-7th, coordinate of each particle
  % Note: it would be more efficient if bigger particles are inserted first, then smaller particles
  % Initialize the cordinate matrix with one biggest Ag particle, Particle's boundary should not outside the wall
[sys, carbon_stats] = insertCarbonClusters(sys, carbon);
sys.cord  = insrtPtc(sys, drx, drx.typ, drx.den, 1);     % Insert all DRX particles
gst.max   = min([carbon.dia, min(drx.dia)]);

% 4. Generate the coordinates of all ghost particles
gst.nxy  = ceil(3*(sqrt(2)-1) * sys.lx/gst.max);  % The number of ghost particles in X and Y direction for current particle size
gst.dia  = sys.lx/gst.nxy;                      % The diameter of ghost particles for current particle size
gst.cord = zeros(gst.nxy^2,7);                  % initialize the ghost particles
gst.cord(:,1) = [size(sys.cord,1)+1 : size(sys.cord,1)+gst.nxy^2];   % Ghost particles index after all Ag and C particles
gst.cord(:,2) = gst.typ;
gst.cord(:,3) = gst.dia;
gst.cord(:,4) = gst.den;
gst.cord(:,7) = sys.lz + gst.dia/2;              % All ghost particles located at the top of the system with a radius distance
for ig = 1 : gst.nxy
    for jg = 1 : gst.nxy
        gst.cord((ig-1)*gst.nxy+jg,5) = gst.dia/2 + gst.dia*(jg-1);
        gst.cord((ig-1)*gst.nxy+jg,6) = gst.dia/2 + gst.dia*(ig-1);
    end
end

% 5. Output all particles in LAMMPS data file format
stmp = strcat(fid,'mdl.data');
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

if ~isempty(carbon_stats)
    fprintf('Carbon cluster placement: inserted %d, failed %d\n', ...
        carbon_stats.inserted_clusters, carbon_stats.failed_clusters);
end
end

%% ===== Helper functions for carbon clusters =====
function [sys, stats] = insertCarbonClusters(sys, carbon)
    stats = struct('name', carbon.cluster.name, 'inserted_clusters', 0, ...
                   'failed_clusters', 0);
    if carbon.nClusters <= 0
        return;
    end
    fprintf('Inserting %d %s clusters...\n', carbon.nClusters, carbon.cluster.name);
    for c = 1:carbon.nClusters
        [placed, newAtoms] = placeCarbonCluster(sys, carbon, carbon.cluster);
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

function [placed, newAtoms] = placeCarbonCluster(sys, carbon, cluster)
    placed = false;
    newAtoms = [];
    max_attempts = 5000;
    offsets = cluster.offsets;
    for attempt = 1:max_attempts
        if cluster.allow_rotation
            rot = randomRotationMatrix();
            offs = (rot * offsets')';
        else
            offs = offsets;
        end
        [range, valid] = computeCenterRanges(sys, offs);
        if ~valid
            return;
        end
        center = [randInRange(range(1,:)), randInRange(range(2,:)), randInRange(range(3,:))];
        positions = offs + center;
        if canPlace(sys, positions, carbon.dia)
            placed = true;
            newAtoms = buildAtomBlock(sys, positions, carbon);
            break;
        end
    end
end

function [range, valid] = computeCenterRanges(sys, offsets)
    min_off = min(offsets, [], 1);
    max_off = max(offsets, [], 1);
    range = [0 - min_off(1), sys.lx - max_off(1);
             0 - min_off(2), sys.ly - max_off(2);
             0 - min_off(3), sys.lz - max_off(3)];
    valid = all(range(:,2) > range(:,1));
end

function val = randInRange(bounds)
    val = bounds(1) + rand * (bounds(2) - bounds(1));
end

function tf = canPlace(sys, positions, dia)
    tf = true;
    if isempty(sys.cord)
        return;
    end
    for i = 1:size(positions,1)
        pos = positions(i,:);
        mnBnd = pos - dia;
        mxBnd = pos + dia;
        nearby = sys.cord(sys.cord(:,5) >= mnBnd(1) & sys.cord(:,5) <= mxBnd(1) & ...
                          sys.cord(:,6) >= mnBnd(2) & sys.cord(:,6) <= mxBnd(2) & ...
                          sys.cord(:,7) >= mnBnd(3) & sys.cord(:,7) <= mxBnd(3), :);
        if isempty(nearby)
            continue;
        end
        dist_sq = (pos(1) - nearby(:,5)).^2 + (pos(2) - nearby(:,6)).^2 + (pos(3) - nearby(:,7)).^2;
        min_dist_sq = ((dia + nearby(:,3))/2).^2;
        if any(dist_sq < min_dist_sq * 0.999)
            tf = false;
            return;
        end
    end
end

function block = buildAtomBlock(sys, positions, carbon)
    n = size(positions,1);
    block = zeros(n,7);
    start_idx = size(sys.cord,1);
    for i = 1:n
        block(i,:) = [start_idx + i, carbon.typ, carbon.dia, carbon.den, positions(i,:)];
    end
end

function R = randomRotationMatrix()
    u1 = rand;
    u2 = rand * 2*pi;
    u3 = rand * 2*pi;
    q = [sqrt(1-u1)*sin(u2), sqrt(1-u1)*cos(u2), sqrt(u1)*sin(u3), sqrt(u1)*cos(u3)];
    qw = q(4); qx = q(1); qy = q(2); qz = q(3);
    R = [1 - 2*(qy^2 + qz^2),     2*(qx*qy - qz*qw),     2*(qx*qz + qy*qw);
         2*(qx*qy + qz*qw), 1 - 2*(qx^2 + qz^2),         2*(qy*qz - qx*qw);
         2*(qx*qz - qy*qw),     2*(qy*qz + qx*qw), 1 - 2*(qx^2 + qy^2)];
end
