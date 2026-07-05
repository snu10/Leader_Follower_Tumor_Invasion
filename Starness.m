function [DistDiff, PeakDeg] = Starness(Sulcus_theta , Sulcus_r, Gyrus_theta, Gyrus_r)


DistDiff = mean(Gyrus_r) - mean(Sulcus_r);

if Sulcus_theta(1) < Gyrus_theta(1) % leftmost point is sulcus

    Sulcus_theta_N = [Sulcus_theta, (Sulcus_theta(1) + 2*pi) ];
    Sulcus_r_N = [Sulcus_r, Sulcus_r(1)];

else    % leftmost point is gyrus 

    Sulcus_theta_N = [(Sulcus_theta(end)-2*pi), Sulcus_theta ];
    Sulcus_r_N = [Sulcus_r(end),  Sulcus_r];

end

if Sulcus_theta_N(end) < Gyrus_theta(end)  % Rightmost is gyrus

        Sulcus_theta_N = [Sulcus_theta_N,  (Sulcus_theta(1)+2*pi)];
    Sulcus_r_N = [ Sulcus_r_N,  Sulcus_r_N(1)];

end


% Sulcus_theta_N
% Sulcus_r_N
% 
% Gyrus_theta
% Gyrus_r

n_peaks = length(Gyrus_theta);

PeakDeg = 0;


for z = 1:n_peaks

    SulcusAng = Sulcus_theta_N(z+1) - Sulcus_theta_N(z);

    SulcusAvrR = (   Sulcus_r_N(z) + Sulcus_r_N(z+1) )/2;

    BaseWidth = 2*SulcusAvrR*sin(SulcusAng/2);


    ArmHeight = Gyrus_r(z) - SulcusAvrR;

      % ArmHeight / BaseWidth

    PeakDeg = PeakDeg  + ArmHeight / BaseWidth;


end

PeakDeg = PeakDeg/n_peaks;

