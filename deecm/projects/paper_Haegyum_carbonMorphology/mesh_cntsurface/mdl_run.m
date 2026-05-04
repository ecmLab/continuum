%% mdl_run.m
%  Array-job version. One Slurm array task = one mass ratio.
%  Each task gets the iwt index from $SLURM_ARRAY_TASK_ID and runs
%  exactly one case, then exits. No parpool, no parfor.
%%% By: Howard Tu, 03/11/2026   |   HPC array-job version
clc; clear;
% Unit: pg, um, us, pg/um^3

%% ============================================================
%  USER-CONFIGURABLE PARAMETERS
%% ============================================================
LD_ratio       = 10;
load_mg_cm2    = 10.0;
box_lz_um      = 30;
olap_cc        = 0.85;
olap_dc        = 0.85;
intra_cc       = 0.85;
wrap_perturb   = 0.30;
min_LD         = 5;

%% ============================================================
%  SELECT WHICH CASE TO RUN
%% ============================================================
arr = getenv('SLURM_ARRAY_TASK_ID');
if isempty(arr)
    % Fallback: allow manual runs e.g. `IWT=3 matlab -batch mdl_run`
    arr = getenv('IWT');
end
if isempty(arr)
    error(['mdl_run:noTaskId', ...
           'No SLURM_ARRAY_TASK_ID (or IWT) set. ', ...
           'Submit with: sbatch --array=1-11 run_matlab.sh']);
end
iwt = str2double(arr);
if isnan(iwt) || iwt < 1 || iwt ~= round(iwt)
    error('mdl_run:badTaskId','Invalid task id "%s" (need positive integer).', arr);
end

%% ============================================================
%  ENSURE OUTPUT DIR EXISTS  (idempotent — safe across array tasks)
%% ============================================================
out_dir = fullfile('massratio', sprintf('mr%d', iwt));
if ~exist('massratio','dir'), mkdir('massratio'); end
if ~exist(out_dir,'dir'),     mkdir(out_dir);     end

%% ============================================================
%  RUN
%% ============================================================
t0 = tic;
fprintf('====== Starting case mr%d on host %s ======\n', iwt, getHostname());
runOneCase(iwt, LD_ratio, load_mg_cm2, box_lz_um, ...
           olap_cc, olap_dc, intra_cc, wrap_perturb, min_LD);
fprintf('====== Case mr%d completed in %.2f minutes ======\n', iwt, toc(t0)/60);


%% ===== Helpers ==============================================

function h = getHostname()
    [s, h] = system('hostname -s');
    if s ~= 0, h = 'unknown'; else, h = strtrim(h); end
end

%% ===== Per-case driver (one mass ratio) =====
function runOneCase(iwt, LD_ratio, load_mg_cm2, box_lz_um, ...
                    olap_cc, olap_dc, intra_cc, wrap_perturb, min_LD)

    rng(iwt, 'twister');   % reproducible per case

    [drx, carbon, gst] = particle_info(LD_ratio);

    carbon.olap_cc      = olap_cc;
    carbon.olap_dc      = olap_dc;
    carbon.intra_cc     = intra_cc;
    carbon.wrap_perturb = wrap_perturb;
    carbon.cluster      = rebuildCNTCluster(carbon.dia, carbon.length, carbon.intra_cc);

    [sys, drx, carbon] = create_sys(drx, carbon, iwt, box_lz_um, load_mg_cm2);
    sys.cord = zeros(0, 7);

    fid = fullfile('massratio', sprintf('mr%d', iwt));

    fprintf('[mr%d] CNT initial L/D=%g  L=%.3f um, %d beads/rod, spacing=%.4f um\n', ...
        iwt, carbon.LD_ratio, carbon.length, carbon.cluster.sphere_count, ...
        carbon.dia*carbon.intra_cc);
    fprintf('[mr%d] DRX=%d, target carbon spheres=%d (%d clusters)\n', ...
        iwt, drx.nTot, carbon.nTot, carbon.nClusters);
    fprintf('[mr%d] Initial target porosity: %6.4f\n', ...
        iwt, 1-(drx.rVol+carbon.rVol)/sys.vol);

    % 1. Insert DRX
    sys.cord = insrtPtc(sys, drx, drx.typ, drx.den, 1);

    % 2. Iteratively wrap CNT rods
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

        fprintf('[mr%d] CNT batch: L/D=%d, L=%.3f um, %d beads/rod, target %d rods\n', ...
            iwt, current_LD, carbon.length, carbon.cluster.sphere_count, nNeeded);

        [sys, batch_stats] = insertCarbonOnDRXSurface(sys, carbon, iwt);

        placed_in_batch      = batch_stats.inserted_clusters * carbon.cluster.sphere_count * sphere_volume;
        placed_carbon_volume = placed_carbon_volume + placed_in_batch;
        total_inserted       = total_inserted + batch_stats.inserted_clusters;
        total_failed         = total_failed + batch_stats.failed_clusters;

        fprintf('[mr%d]   Placed %d/%d clusters. Carbon vol %.4f / %.4f um^3 (%.2f%%)\n', ...
            iwt, batch_stats.inserted_clusters, nNeeded, ...
            placed_carbon_volume, target_carbon_volume, ...
            100*placed_carbon_volume/max(target_carbon_volume,eps));

        if placed_carbon_volume < target_carbon_volume
            if current_LD <= min_LD
                fprintf('[mr%d] Warning: target carbon vol not met at min L/D=%d. Stopping.\n', ...
                    iwt, min_LD);
                break;
            end
            current_LD = current_LD - 5;
        end
    end

    % 3. Ghost layer
    gst.max  = min([carbon.dia, min(drx.dia)]);
    gst.nxy  = ceil(3*(sqrt(2)-1) * sys.lx/gst.max);
    gst.dia  = sys.lx/gst.nxy;
    gst.cord = zeros(gst.nxy^2, 7);
    gst.cord(:,1) = (size(sys.cord,1)+1 : size(sys.cord,1)+gst.nxy^2)';
    gst.cord(:,2) = gst.typ;
    gst.cord(:,3) = gst.dia;
    gst.cord(:,4) = gst.den;
    gst.cord(:,7) = sys.lz + gst.dia/2;
    for ig = 1:gst.nxy
        for jg = 1:gst.nxy
            gst.cord((ig-1)*gst.nxy+jg, 5) = gst.dia/2 + gst.dia*(jg-1);
            gst.cord((ig-1)*gst.nxy+jg, 6) = gst.dia/2 + gst.dia*(ig-1);
        end
    end

    % 4. Output LAMMPS data file
    out_path = fullfile(fid, 'mdl.data');
    fileID = fopen(out_path, 'w');
    if fileID < 0
        error('Could not open %s for writing.', out_path);
    end
    fprintf(fileID, 'LAMMPS data file\n\n');
    fprintf(fileID, '%6d %6s\n', size(sys.cord,1)+size(gst.cord,1), 'atoms');
    fprintf(fileID, '\n   3 atom types\n\n');
    fprintf(fileID, '%23.19f %23.19f %3s %3s\n', 0.0, sys.lx, 'xlo', 'xhi');
    fprintf(fileID, '%23.19f %23.19f %3s %3s\n', 0.0, sys.ly, 'ylo', 'yhi');
    fprintf(fileID, '%23.19f %23.19f %3s %3s\n', 0.0, sys.lz+gst.dia, 'zlo', 'zhi');
    fprintf(fileID, '\nAtoms # sphere\n\n');
    fprintf(fileID, '%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n', sys.cord');
    fprintf(fileID, '%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n', gst.cord');
    fclose(fileID);

    fprintf('[mr%d] Final: clusters placed=%d, failed=%d, final L/D=%d\n', ...
        iwt, total_inserted, total_failed, current_LD);
    fprintf('[mr%d] Carbon vol realized: %.4f / %.4f um^3 (%.2f%%)\n\n', ...
        iwt, placed_carbon_volume, target_carbon_volume, ...
        100*placed_carbon_volume/max(target_carbon_volume,eps));
end


%% ===== CNT cluster rebuilder (matches particle_info) =====
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


%% ===== Surface-attached CNT placement (with spatial hash) =====
function [sys, stats] = insertCarbonOnDRXSurface(sys, carbon, iwt)
    stats = struct('name','CNT_surface','inserted_clusters',0,'failed_clusters',0);
    if carbon.nClusters <= 0
        return;
    end
    drx_particles = sys.cord(sys.cord(:,2) == 1, :);
    if isempty(drx_particles)
        error('No DRX particles found. Cannot attach CNTs to surfaces.');
    end

    max_dia   = max(max(sys.cord(:,3)), carbon.dia);
    cell_size = max_dia;
    cluster_size   = carbon.cluster.sphere_count;
    extra_capacity = carbon.nClusters * cluster_size + 64;
    H = hashInit(sys.lx, sys.ly, sys.lz + cell_size, cell_size, ...
                 size(sys.cord,1) + extra_capacity);
    H = hashAddBatch(H, sys.cord(:,5:7), sys.cord(:,3), sys.cord(:,2));

    fprintf('[mr%d] Wrapping %d CNT rods (L/D=%g)...\n', ...
        iwt, carbon.nClusters, carbon.LD_ratio);

    new_atoms = zeros(carbon.nClusters * cluster_size, 7);
    new_count = 0;
    sys_n0    = size(sys.cord, 1);

    for c = 1:carbon.nClusters
        [placed, positions] = placeCNTTangentially(sys, carbon, drx_particles, H);
        if placed
            n   = size(positions, 1);
            ids = (sys_n0 + new_count + 1 : sys_n0 + new_count + n).';
            block = [ids, ...
                     repmat(carbon.typ, n, 1), ...
                     repmat(carbon.dia, n, 1), ...
                     repmat(carbon.den, n, 1), ...
                     positions];
            new_atoms(new_count+1 : new_count+n, :) = block;
            H = hashAddBatch(H, positions, ...
                             repmat(carbon.dia, n, 1), ...
                             repmat(carbon.typ, n, 1));
            new_count = new_count + n;
            stats.inserted_clusters = stats.inserted_clusters + 1;
        else
            stats.failed_clusters = stats.failed_clusters + 1;
        end
        if mod(c, 100) == 0
            fprintf('[mr%d]   progress %d/%d\n', iwt, c, carbon.nClusters);
        end
    end

    if new_count > 0
        sys.cord = [sys.cord; new_atoms(1:new_count, :)];
    end
end

function [placed, positions] = placeCNTTangentially(sys, carbon, drx_particles, H)
    placed       = false;
    positions    = [];
    max_attempts = 5000;
    nDRX    = size(drx_particles, 1);
    c_r     = carbon.dia / 2;
    nSeg    = carbon.cluster.sphere_count;
    spacing = carbon.dia * carbon.intra_cc;

    for attempt = 1:max_attempts %#ok<NASGU>
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

        if canPlaceHash(H, positions, carbon.dia, carbon.olap_cc, carbon.olap_dc)
            placed = true;
            return;
        end
    end
    positions = [];
end

function positions = randomGeodesicWalk(center, R, nSeg, d_theta, perturb)
    u = randn(1,3); u = u / norm(u);
    t = randn(1,3); t = t - dot(t,u)*u; t = t / norm(t);
    positions = zeros(nSeg, 3);
    cd = cos(d_theta); sd = sin(d_theta);
    for i = 1:nSeg
        positions(i,:) = center + R*u;
        if i == nSeg, break; end
        u_new = u*cd + t*sd;
        t_new = -u*sd + t*cd;
        u = u_new / norm(u_new);
        t = t_new - dot(t_new, u)*u;
        t = t / norm(t);
        if perturb > 0
            phi = (rand-0.5) * 2 * perturb;
            cp = cos(phi); sp = sin(phi);
            t = t*cp + cross(u,t)*sp;
            t = t - dot(t,u)*u;
            t = t / norm(t);
        end
    end
end

function tf = hasSelfCollision(positions, spacing)
    n = size(positions, 1);
    if n < 3
        tf = false;
        return;
    end
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


%% ===== Spatial hash =====
function H = hashInit(lx, ly, lz, cell_size, capacity)
    H.cs = cell_size;
    H.nx = max(1, ceil(lx / cell_size));
    H.ny = max(1, ceil(ly / cell_size));
    H.nz = max(1, ceil(lz / cell_size));
    H.cells = cell(H.nx, H.ny, H.nz);
    H.pos = zeros(capacity, 3);
    H.dia = zeros(capacity, 1);
    H.typ = zeros(capacity, 1);
    H.n   = 0;
end

function H = hashAddBatch(H, positions, diameters, types)
    n = size(positions, 1);
    if n == 0, return; end
    cap = size(H.pos, 1);
    if H.n + n > cap
        new_cap = max(cap*2, H.n + n + 1024);
        H.pos(new_cap, 3) = 0;
        H.dia(new_cap, 1) = 0;
        H.typ(new_cap, 1) = 0;
    end
    rng_idx = H.n + (1:n).';
    H.pos(rng_idx, :) = positions;
    H.dia(rng_idx)    = diameters;
    H.typ(rng_idx)    = types;

    ix = max(1, min(H.nx, ceil(positions(:,1) / H.cs)));
    iy = max(1, min(H.ny, ceil(positions(:,2) / H.cs)));
    iz = max(1, min(H.nz, ceil(positions(:,3) / H.cs)));
    for k = 1:n
        H.cells{ix(k), iy(k), iz(k)}(end+1, 1) = rng_idx(k);
    end
    H.n = H.n + n;
end

function tf = canPlaceHash(H, positions, dia, olap_cc, olap_dc)
    tf = true;
    if H.n == 0, return; end
    cs = H.cs;

    for i = 1:size(positions, 1)
        pos = positions(i, :);
        ixmin = max(1,    ceil((pos(1) - dia) / cs));
        ixmax = min(H.nx, ceil((pos(1) + dia) / cs));
        iymin = max(1,    ceil((pos(2) - dia) / cs));
        iymax = min(H.ny, ceil((pos(2) + dia) / cs));
        izmin = max(1,    ceil((pos(3) - dia) / cs));
        izmax = min(H.nz, ceil((pos(3) + dia) / cs));
        if ixmin > ixmax || iymin > iymax || izmin > izmax
            continue;
        end

        nC = (ixmax-ixmin+1) * (iymax-iymin+1) * (izmax-izmin+1);
        bufs = cell(nC, 1);
        k = 0;
        for ix = ixmin:ixmax
            for iy = iymin:iymax
                for iz = izmin:izmax
                    k = k + 1;
                    bufs{k} = H.cells{ix, iy, iz};
                end
            end
        end
        idxs = vertcat(bufs{:});
        if isempty(idxs), continue; end

        nb_pos = H.pos(idxs, :);
        nb_dia = H.dia(idxs);
        nb_typ = H.typ(idxs);

        in_box = nb_pos(:,1) >= pos(1)-dia & nb_pos(:,1) <= pos(1)+dia & ...
                 nb_pos(:,2) >= pos(2)-dia & nb_pos(:,2) <= pos(2)+dia & ...
                 nb_pos(:,3) >= pos(3)-dia & nb_pos(:,3) <= pos(3)+dia;
        if ~any(in_box), continue; end
        nb_pos = nb_pos(in_box, :);
        nb_dia = nb_dia(in_box);
        nb_typ = nb_typ(in_box);

        dist_sq     = (pos(1)-nb_pos(:,1)).^2 + (pos(2)-nb_pos(:,2)).^2 + (pos(3)-nb_pos(:,3)).^2;
        min_dist_sq = ((dia + nb_dia)/2).^2;

        tol = ones(numel(nb_dia), 1);
        tol(nb_typ == 2) = olap_cc;
        tol(nb_typ == 1) = olap_dc;

        if any(dist_sq < min_dist_sq .* tol)
            tf = false;
            return;
        end
    end
end