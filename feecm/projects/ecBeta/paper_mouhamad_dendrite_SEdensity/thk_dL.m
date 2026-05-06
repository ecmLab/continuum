%% This function is used to calculate the current thickness of Li and increment
function [dtLi, dhLi] = thk_dL(it,dt,wSE,yLi,dlt)
    VLi_F = 1.347*10^-4;   % MolarVolume_Li/F, in unit cm^3/(s*A)
    dtLi  = zeros(length(yLi),1);
    
% Step 1: Load CSV data from MOOSE calculation: structure of the csv data
% Column 1-2: current density along the x and y direction, in unit mA/cm^2
% Column 4-5: coordinate in the x and y direction, in unit um
    tmp         = csvread(['rst/t',num2str(it),'.csv'],1,0);
    [cord,indx] = sortrows(tmp(:,4:5),[1,2],{'ascend' 'ascend'}); % sort the coordinate with x ascend and y dscend
    crnt        = tmp(indx,1:2); % Sort the current density accordingly
    
% Step 2: Calculate the normal current density on pore surface
    cordY  = cord(cord(:,1)==wSE,2);   % The y-coodinate of node points, x is always wSE
    nrCnt  = crnt(cord(:,1)==wSE,1);   % The normal current density along the pore surface, in unit mA/cm^2
    
% Step 3: Calculate Li thickness incremental from normal current, and interpolate to sampling point
    dtLi0  = nrCnt*VLi_F*dt*10*3600;   % The Li thickness increasement from current density, in unit um
    iy     = yLi(yLi>cordY(1));
    dtLi1  = interp1(cordY,dtLi0,iy);
    dtLi(yLi>cordY(1)) = dtLi1;
    
% Step 4: Calculate Li length incremental and update Li thickness incremental
% Assuming Li within the distance delta is contributed to Li length growth
    idx    = find(iy<cordY(1)+dlt);        % Find the sampling points within a distance delta to the tip
    dhLi   = trapz(iy(idx),dtLi1(idx))/dtLi1(idx(end)) - dlt;  % The Li length incremental
    
    idd    = find(yLi>cordY(1)-dhLi & yLi<cordY(1)+dlt);
    dtLi(idd) = dtLi1(idx(end));
    