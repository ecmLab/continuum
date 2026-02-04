%% The Normal Distribution Function used for particle size distribution

function Prbly = particleDistr(size,aveD,sigD, flg)

%     vol       = 1/6*pi*size.^3;
%     aveD      = mean(vol);
%     sigD      = std2(vol);
    
    if flg == 0      % Gaussian distribution
        Prbly = 1/(sigD*sqrt(2*pi)) .* exp(-(size-aveD).^2./(2*sigD^2));
    elseif flg == 1      % Log-normal distribution
        aveE  = log(aveD/(sqrt(1 + sigD^2/aveD^2)));  % equivalent mean value
        sigE  = sqrt(log(1 + sigD^2/aveD^2));   % equivalent deviation value
        Prbly = 1./(size*sigE*sqrt(2*pi)) .* exp(-(log(size)-aveE).^2./(2*sigE^2));
    else
        disp('Wrong Input for Distribution');
    end
end