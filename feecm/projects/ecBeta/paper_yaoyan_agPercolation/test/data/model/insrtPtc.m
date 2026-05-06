%% Function: Insert all Particles into the box
% mtd = 1, direct search method
% mtd = 2, Linked-List Grid Search Method

function cord = insrtPtc(sys, ptc)

rng('shuffle');     % Avoid repeating the same random number

lx      = sys.lx; 
ly      = sys.ly;
lz      = sys.lz;
cord    = sys.cord;
rds     = ptc.rds;
Np      = ptc.nRds;

idx      = size(cord,1);                    % The current number of particles in the box
ipt      = 0;
% Initially only one particle exists, so the biggest NCM particle number should - 1
if idx == 1                            
    Np(end) = Np(end) - 1;
end
% Start the Insertion Process
for i = length(rds):-1:1
    lRds  = rds(i);                        % The diameter of to be inserted particle
    msh   = (max(rds)+lRds)*[1,1,1];       % The mesh size of each particle
    for j = 1 : Np(i)
        flg  = 1;                          % The flag of keeping on searching
        while flg
            try_p = [(lx - 4*lRds)*rand + 2*lRds,(ly - 4*lRds)*rand + 2*lRds, (lz - 4*lRds)*rand + 2*lRds]; % Try a point
        % Direct Search Method
                dis   = sqrt((try_p(1)-cord(1:idx,3)).^2 + ...
                             (try_p(2)-cord(1:idx,4)).^2 + ...
                             (try_p(3)-cord(1:idx,5)).^2);
                diff  = dis(1:idx) - (lRds + cord(1:idx,2));
                if diff >= -lRds
                    idx         = idx + 1;
                    ipt         = ipt + 1;
                    cord(idx,:) = [idx, lRds, try_p];
                    flg         = 0;
                end
        end
        
        if not(mod(idx,20))
            fprintf('Progress: %d of %d\n',ipt, sum(ptc.nRds));
        end
    end
end

