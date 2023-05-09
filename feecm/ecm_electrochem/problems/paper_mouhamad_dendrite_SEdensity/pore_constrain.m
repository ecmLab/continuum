%% This function is used to calculate the current thickness of Li and increment
function [tLiN, hLiN] = pore_constrain(ySE,pd,tLiP, hLiP)
tLiN  = tLiP;
hLiN  = hLiP;

% Step 1: Update the Li thickness and Li length if the Li thickness is larger than the pore radius pd/2
iOvr   = find(tLiP>pd/2);      % Find the sampling points with Li thickness beyond pore radius
idx_Li = find(tLiP>0);         % The sampling points of Li on pore surface at next timestep
pRoot  = idx_Li(1);            % Find the location of the root point at next timestep
j      = 0;
if length(iOvr)>1                % If Li thickness is indeed beyond pore radius
    aOvr = trapz(ySE(iOvr(1):iOvr(end)),tLiP(iOvr(1):iOvr(end))-pd/2); % Area beyond limit
    tLiN(iOvr(1):iOvr(end)) = pd/2;  % Set the Li thickness beyond pore radius as the pore radius
    if aOvr > 0                 % If the overall area beyond pore radius is positive
        vr = 0; k=0; aResidual=0;
        while vr<aOvr/2         % Li grow along anode direction
            if k > length(ySE)-iOvr(end)-1
                aResidual = aOvr/2 - vr;
                break;
            end
            k  = k + 1;
            vr = pd/2*(ySE(iOvr(end)+k)-ySE(iOvr(end))) - trapz(ySE(iOvr(end):iOvr(end)+k),tLiP(iOvr(end):iOvr(end)+k));
        end
        
        vl=0; j=0;
        while vl<aOvr-vr+aResidual         % Li grow along cathode direction
            if j > iOvr(1)-1
                break;
            end
            j  = j + 1;
            vl = pd/2*(ySE(iOvr(1))-ySE(iOvr(1)-j)) - trapz(ySE(iOvr(1)-j:iOvr(1)),tLiP(iOvr(1)-j:iOvr(1)));
        end
        
        tLiN(iOvr(1)-j+1:iOvr(end)+k) = pd/2;  % Set the Li thickness beyond pore radius as the pore radius
    end
    
    % Step 2: Calculate Li length incremental and update Li length
    if j>iOvr(1)-pRoot % If Li grow beyond the root point, the Li length is determined by Li extrusion
        hLiN = ySE(end)-ySE(iOvr(1)-j);
    end
    
end



