SavePNG = 0;
VisStar = 1;

densityMat_tot_G = P_density_G + C_density_G + A_density_G + E_density_G + D_density_G;

densityMat_Liv_G = P_density_G + C_density_G + A_density_G + E_density_G ;

densityMat_tot_C0 = gather(densityMat_tot_G);
densityMat_tot_C = densityMat_tot_C0;

densityMat_Liv_C0 = gather(densityMat_Liv_G);
densityMat_Liv_C = densityMat_Liv_C0;


TumMass_tot = sum(densityMat_tot_C0, 'all')*(sc^2);
TumMass_liv = sum(densityMat_Liv_C0, 'all')*(sc^2);


densityMat_tot_C(densityMat_tot_C < 0.005) = 0;
densityMat_Liv_C(densityMat_Liv_C < 0.005) = 0;

tumorMask = densityMat_tot_C > 0;

CC = bwconncomp(tumorMask);

n_Cluster = CC.NumObjects;

if n_Cluster > 1

    Cluster_Size = zeros(1, n_Cluster);

    for nc0 = 1:n_Cluster

        Cluster_Size(nc0) = numel(CC.PixelIdxList{nc0});

    end

    [~, L_Clst] = max(Cluster_Size);

    ZeroMask = zeros(size(densityMat_tot_C0));
    ZeroMask(CC.PixelIdxList{L_Clst}) = 1;
    ZeroMask = (ZeroMask > 0);

    tumorMask = ZeroMask;

end
TumArea = sum(tumorMask(:))*(sc^2);


BW = imfill(tumorMask, 'holes');

boundaries = bwboundaries(BW);
% Assume the object of interest is the first one:
boundary = boundaries{1};  % boundary is an N-by-2 matrix: [row, col]

% 2. Compute the centroid of the shape using regionprops
stats = regionprops(BW, 'Centroid');
center =  size(densityMat_tot_C0)/2;  % [x, y] coordinate

% 3. Convert boundary coordinates to polar coordinates relative to the center
% Get x (column) and y (row) coordinates from the boundary:
y_coords = boundary(:,1);  % row indices
x_coords = boundary(:,2);  % column indices

% Compute differences relative to the center
dx = x_coords - center(1);
dy = y_coords - center(2);

% Convert to polar coordinates (theta in radians, r is the radial distance)
[theta0, r0] = cart2pol(dx, dy);

theta = theta0;
r = r0*sc;

overlap_pos  = find( (diff(theta0) <= 0) & (diff(theta0) > -pi));
theta(overlap_pos) = [];
r(overlap_pos) = [];

for jj = 1:20

    overlap_pos  = find( (diff(theta) <= 0) & (diff(theta) > -pi));
    theta(overlap_pos) = [];
    r(overlap_pos) = [];

end

[theta_sorted, idx] = sort(theta);
r_sorted = r(idx);

theta_EXT = [(theta_sorted - 2*pi); theta_sorted; (theta_sorted + 2*pi)];
r_EXT = [r_sorted; r_sorted; r_sorted];

[unique_x, ~, idx] = unique(theta_EXT);
counts = histcounts(idx, 0.5:1:max(idx)+0.5);

repeated_x = unique_x(counts > 1);

JM = [];

for kk = 1:length(repeated_x)

    xPos = find(theta_EXT == repeated_x(kk));

    y_min = min(r_EXT(xPos))*ones(size(xPos));   % Choose smallest value when repeated

    r_EXT(xPos) = y_min;

    JM = [JM, xPos(2:end)];

end

theta_EXT(JM) = [];
r_EXT(JM)=[];

if VisStar == 1

    % Plot the distance vs. theta (using radians or degrees as preferred)
    PeakFig = figure('position', [500 300 1000 600]);
    plot(theta_EXT, r_EXT,'k');
    xlabel('Theta (radians)');
    ylabel('Distance from center');
    title([ Fig_title ,  ': Distance vs. Theta']);

end

theta_uniform = linspace(min(theta_EXT), max(theta_EXT), 20000);  % or finer if needed
r_interp_lin = interp1(theta_EXT, r_EXT, theta_uniform, 'linear');

r_interp = smoothdata(r_interp_lin, 'gaussian', 60);  % adjust window size

if VisStar == 1

    hold on
    plot(theta_uniform, r_interp, 'r')

end

Max_r = max(r_EXT);
min_r = min(r_EXT);

minProminence = (Max_r - min_r)*0.1;  % adjust this value based on data traits

minProminence = min(minProminence, 1);

[peakVals, locs, widths1, prominences1] = ...
    findpeaks(r_interp, theta_uniform, 'MinPeakProminence', minProminence, 'MinPeakDistance', 0.05  );

[peakValsN, locsN, widths2, prominences2] = ...
    findpeaks(-r_interp, theta_uniform, 'MinPeakProminence', minProminence, 'MinPeakDistance', 0.05);


showIDX1 = find((locs <= pi) & (locs>= -pi));
showIDX2 = find((locsN <= pi) & (locsN>= -pi));

if VisStar == 1

    plot(locs(showIDX1), peakVals(showIDX1), 'ro', 'MarkerFaceColor','r');
    plot(locsN(showIDX2), -peakValsN(showIDX2), 'bo', 'MarkerFaceColor','b');
    ylim([0 Inf])
    xlim([-pi pi])

end

Gyrus_theta = locs(showIDX1);
Gyrus_r = peakVals(showIDX1);

Sulcus_theta = locsN(showIDX2);
Sulcus_r = -peakValsN(showIDX2);

if length(peakVals) >= 4

    [DistDiff, PeakScore] = Starness(Sulcus_theta , Sulcus_r, Gyrus_theta, Gyrus_r);

else

    PeakScore= 0;
    DistDiff = 0;

end

MassMap_liv = densityMat_Liv_C0 *(sc^2);  % Normalize: now a discrete probability distribution

avr_x = sum(MassMap_liv .* x_grid, 'all') / TumMass_liv;
avr_y = sum(MassMap_liv .* y_grid, 'all') / TumMass_liv;

var_x = sum(MassMap_liv .* (x_grid - avr_x).^2, 'all') ;
var_y = sum(MassMap_liv .* (y_grid - avr_y).^2, 'all') ;

Inv_score = var_x+var_y;

if  VisStar == 1

    ResFig = figure('position', [500 300 700 700], 'color','white');
    imagesc(densityMat_tot_C0)
    xlabel('X-grids');
    ylabel('Y-grids');
    title([ Fig_title ]);

    hold on

end

[r_num0, c_num0] = size(densityMat_tot_C);

cent_val = r_num0/2;

if VisStar == 1

    [Xs, Ys] = pol2cart(Sulcus_theta , Sulcus_r/sc);

    plot(Xs+cent_val,Ys+cent_val,  'bo', 'markerfacecolor', [0.3 0.3 1], 'markeredgecolor','none')


    [Xg, Yg] = pol2cart(Gyrus_theta , Gyrus_r/sc);

    plot(Xg+cent_val,Yg+cent_val,  'ro', 'markerfacecolor', 'r','markeredgecolor','none' )

    colorbar('Ticks', [0:0.5:5])
    clim([0 2])

    colormap pink
    axis square

end

if SavePNG == 1 && VisStar == 1

    v_info_all = version;

    v_info = v_info_all(1:2);

    pngName = [ParamTag, ' ', NT_name, ' ', CTX_name, ' ', AD_name, ' ', ECM_name, '.png' ];

    if strcmp(v_info, '25')

        exportgraphics(ResFig, pngName,'resolution', 600, 'Padding', 100)

    else
        exportgraphics(ResFig, pngName,'resolution', 600)

    end

end

ParamMat = [CA_cond, CTX_cond,  ECM_cond];

if ECM_cond < 10

    ParamMat = [CA_cond, CTX_cond,  0, ECM_cond];

end

ResMat = [TumMass_tot, Max_r, DistDiff, Inv_score, PeakScore];

RecordMat = [ParamMat, ResMat];

% P_density_C = gather(P_density_G) ;
A_density_C = gather(A_density_G) ;
% C_density_C = gather(C_density_G) ;
% E_density_C = gather(E_density_G) ;

Gyrus_data = [Gyrus_theta; Gyrus_r];
Sulcus_data = [Sulcus_theta; Sulcus_r];



