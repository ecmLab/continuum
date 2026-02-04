%% This script is designed for generating models for distribution of both DRX and carbon particles
%%% Function: impt_mat: define geometry and material properties of all particles
%%% By: Howard Tu, 06/10/2023 (mod. for DRX/C terminology)
%% Comments from Haegyum:
% We can assume DRX particle size ranging from 500 nm to 4 um. 
% CNT diameter would be < 100 nm and length can range from 10 um to 200 um. 
% Super P could be a spherical shape with 50 - 200 nm. 
% Graphene thickness could be < 50 nm and flake size could range from 10 um x 10 um to 200 um x 200 um. 

function [drx,carbon,gst] = particle_info()
%% 1.Define properties of DRX particles in cathode structure
drx         = struct;
drx.nD      = 40;            % The number of sampling diameters of DRX particles
drx.den     = 4.0;           % The density of DRX material, in unit picogram/micron^3, equivalent to g/cm^3, based on the paper:
                             % https://www.sciencedirect.com/science/article/pii/S2542435121005316
drx.typ     = 1;             % Label DRX particles as Type 1 particles
drx.ave     = 1;             % Designed average diameter of all DRX particles, in μm
drx.min     = 0.5;           % Minimum DRX diameter
drx.max     = 2;             % Maximum DRX diameter
drx.sgm     = 0.2;           % Deviation of DRX particles


%% 2.Define properties of conductive carbon particles (graphene flakes)
carbon         = struct;
carbon.den     = 2.2;           % Carbon density, pg/μm^3 (equivalent to g/cm^3)
carbon.typ     = 2;             % Label carbon particles as Type 2 particles
carbon.dia     = 0.1;           % Building-block diameter of carbon particles, μm
carbon.ave     = carbon.dia;
carbon.cluster = buildGrapheneCluster(carbon.dia, 1.0, 1.0, 0.1); % 1μm x 1μm flake, 100 nm thick

%% 3.Define ghost particles at the top of the system to apply pressure
gst         = struct;
gst.den     = 10;            % The density of ghost material
gst.typ     = 3;

%% 4. Compute theoretical particle probabilities of all particles
drx = cmp_prob(drx,drx.nD,1);      % For DRX particles 

end

function cluster = buildGrapheneCluster(base_dia, lx, ly, thickness)
    spacing = base_dia * 0.9;
    nx = max(1, ceil(lx / spacing));
    ny = max(1, ceil(ly / spacing));
    nz = max(1, ceil(thickness / base_dia));
    x_coords = linspace(-lx/2, lx/2, nx);
    y_coords = linspace(-ly/2, ly/2, ny);
    z_coords = linspace(-thickness/2, thickness/2, nz);
    offsets = zeros(nx*ny*nz,3);
    idx = 0;
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                idx = idx + 1;
                offsets(idx,:) = [x_coords(ix), y_coords(iy), z_coords(iz)];
            end
        end
    end
    cluster.name = 'Graphene';
    cluster.offsets = offsets;
    cluster.sphere_count = size(offsets,1);
    cluster.allow_rotation = true;
    cluster.sphere_volume = (pi/6) * base_dia^3;
    cluster.cluster_volume = cluster.sphere_count * cluster.sphere_volume;
end
