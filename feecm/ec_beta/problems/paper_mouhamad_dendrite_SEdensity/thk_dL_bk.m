%% This function is used to calculate the current thickness of Li and increment
function [tLiN, hLiN] = thk_dL(it,dt,wSE,ySE,dlt,pd,tLiP, hLiP)
    VLi_F = 1.347*10^-4;   % MolarVolume_Li/F, in unit cm^3/(s*A)
    dtLi  = zeros(length(ySE),1);
    
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
    dtLi0  = nrCnt*VLi_F*dt*10*3600;   % The Li thickness incremental from current density, in unit um
    idx_Li0= find(ySE>cordY(1));      % The sampling points of Li on pore surface
    yLi    = ySE(idx_Li0);              % The y coordinates of these sampling points
    dtLi(idx_Li0) = interp1(cordY,dtLi0,yLi); % The interpolated Li thickness incremental
    tLiN   = tLiP + dtLi;              % Compute the new Li thickness
    
% Step 4: Update the Li thickness and Li length if the Li thickness is larger than the pore radius pd/2
    iOvr   = find(tLiN>pd/2);      % Find the sampling points with Li thickness beyond pore radius
    idx_Li = find(tLiN>0);         % The sampling points of Li on pore surface at next timestep
    pRoot  = idx_Li(1);            % Find the location of the root point at next timestep
    j      = 0;
    if length(iOvr)>1                % If Li thickness is indeed beyond pore radius
        aOvr = trapz(ySE(iOvr(1):iOvr(end)),tLiN(iOvr(1):iOvr(end))-pd/2); % Area beyond limit
        tLiN(iOvr(1):iOvr(end)) = pd/2;  % Set the Li thickness beyond pore radius as the pore radius
        if aOvr > 0                 % If the overall area beyond pore radius is positive
            vr = 0; k=0; aResidual=0;
            while vr<aOvr/2         % Li grow along anode direction
                if k > length(ySE)-iOvr(end)-1
                    aResidual = aOvr/2 - vr;
                    break;
                end
                k  = k + 1;
                vr = pd/2*(ySE(iOvr(end)+k)-ySE(iOvr(end))) - trapz(ySE(iOvr(end):iOvr(end)+k),tLiN(iOvr(end):iOvr(end)+k));
            end
            
            vl=0; j=0;
            while vl<aOvr/2+aResidual         % Li grow along cathode direction
                if j > iOvr(1)-1
                    break;
                end
                j  = j + 1;
                vl = pd/2*(ySE(iOvr(1))-ySE(iOvr(1)-j)) - trapz(ySE(iOvr(1)-j:iOvr(1)),tLiN(iOvr(1)-j:iOvr(1)));
            end
            
            tLiN(iOvr(1)-j+1:iOvr(end)+k) = pd/2;  % Set the Li thickness beyond pore radius as the pore radius
        end
    end
    
% Step 5: Calculate Li length incremental and update Li length    
    if length(iOvr)>1 && j>iOvr(1)-pRoot % If Li grow beyond the root point, the Li length is determined by Li extrusion
        hLiN = ySE(end)-ySE(iOvr(1)-j);
    else       % If Li grow inside the root point, Li length growth is based on Li volume within the distance delta 
        idx    = find(ySE>cordY(1) & ySE<cordY(1)+dlt);  % Find the sampling points within a distance delta to the root
        if length(idx)<2                      % In extreme cases, less than two points is included, enlarge the range                                  
            dlt = dlt*1.5;
            idx = find(ySE>=yLi(1) & ySE<yLi(1)+dlt);        
        end
%         dhLi   = trapz(ySE(idx),tLiN(idx))/tLiN(idx(end)) - dlt;  % The Li length incremental
        dhLi   = trapz(ySE(idx),dtLi(idx))/dtLi(idx(end)) - dlt;  % The Li length incremental
        idd    = find(ySE>cordY(1)-dhLi & ySE<cordY(1)+dlt);
        tLiN(idd) = tLiN(idx(end));        % Update Li thicknes based on the assumption
        hLiN   = hLiP + dhLi;              % Compute the new Li length
    end
    
    