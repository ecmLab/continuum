%% Function: Insert all Particles into the box
% mtd = 1, direct search method
% mtd = 2, Linked-List Grid Search Method

function cord = insrtPtc(sys, ptc, type,den, mtd)

rng('shuffle');

lx      = sys.lx; 
ly      = sys.ly;
lz      = sys.lz;
cord    = sys.cord;
dia     = ptc.dia;
Np      = ptc.nDia;

idx      = size(cord,1);                    % The current number of particles in the box
ipt      = 0;
% Initially only one particle exists, so the biggest NCM particle number should - 1
if idx == 1                            
    Np(end) = Np(end) - 1;
end
% Start the Insertion Process
for i = length(dia):-1:1
    lDia  = dia(i);                        % The diameter of to be inserted particle
    msh   = (max(dia)+lDia)*[1,1,1];       % The mesh size of each particle
    for j = 1 : Np(i)
        flg  = 1;                          % The flag of keeping on searching
        while flg
            try_p = [(lx - lDia)*rand + lDia/2,(ly - lDia)*rand + lDia/2,(lz - lDia)*rand + lDia/2]; % Try a point
        % Direct Search Method
            if mtd == 1
                dis   = sqrt((try_p(1)-cord(1:idx,5)).^2 + ...
                             (try_p(2)-cord(1:idx,6)).^2 + ...
                             (try_p(3)-cord(1:idx,7)).^2);
                diff  = dis(1:idx) - (lDia + cord(1:idx,3))/2;
                if diff > 0.0
                    idx         = idx + 1;
                    ipt         = ipt + 1;
                    cord(idx,:) = [idx, type, lDia, den, try_p];
                    flg         = 0;
                end
         % Linked-List Grid Search Method
            elseif mtd == 2
                mnBnd  = try_p - msh/2;               % Lower limit of the mesh
                mxBnd  = try_p + msh/2;               % Upper limit of the mesh 
                if all(mnBnd >= 0) && all(mxBnd <= [lx, ly, lz])  % Only when the particle is inside the box
                    PX   =  cord(cord(1:idx,5) >= mnBnd(1) & cord(1:idx,5) <= mxBnd(1),3:7); % Screen X direction
                    PY   =  PX(PX(:,4) >= mnBnd(2) & PX(:,4) <= mxBnd(2),:);             % Screen Y direction
                    PZ   =  PY(PY(:,5) >= mnBnd(3) & PY(:,5) <= mxBnd(3),:);             % Screen Z direction
                    if isempty(PZ)
                        idx         = idx + 1;
                        cord(idx,:) = [idx, type, lDia, den, try_p];
                        flg         = 0;
                    else
                        dis   = sqrt((try_p(1)-PZ(:,3)).^2 + (try_p(2)-PZ(:,4)).^2 + ...
                                     (try_p(3)-PZ(:,5)).^2);
                        diff  = dis - (lDia + PZ(:,1))/2;
                        if diff >= 0
                            idx         = idx + 1;
                            cord(idx,:) = [idx, type, lDia, den, try_p];
                            flg         = 0;
                        end
                    end
                end
                clear PX PY PZ
            end
            
        end
        
        if not(mod(idx,200))
            fprintf('Progress: %d of %d\n',ipt, ptc.nTot);
        end
    end
end

