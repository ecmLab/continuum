%% Analysis code for box size convergence test

nCfg = 3;
iNcm = 1;
nmr  = [56, 80, 90];
mr   =  (95-nmr)./(nmr);
ld   = 10*10;
ncm_den = 4.85;
lps_den = 1.87;

for icfg = 1 : nCfg
    for imr = 1 : length(nmr)
% Load data
        fname = strcat('c',num2str(icfg),'_ncm',num2str(iNcm),'_lps4_mr',num2str(nmr(imr)),'ld');
        vname = strcat('c',num2str(icfg),'_ncm',num2str(iNcm),'_lps4_mr',num2str(nmr(imr)));
        tmp   = strcat(vname,'=load(''',fname,''');');
        eval(tmp);
        tmp   = strcat('val=',vname);
        eval(tmp);
        
% Compute porosity
        pr = 1 - ld ./ (val(:,4)*ncm_den) .* (1 + ncm_den/lps_den*mr(imr)/(1 + mr(imr)));
        tmp   = strcat(vname,'(:,5) = pr;');
        eval(tmp);
        
    end
end


