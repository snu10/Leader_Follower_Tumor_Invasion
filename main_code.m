%% Parameter Setting
DurRec = 0;

CA_cond = 2;     % 1: HA & LG, 2: MA & MG, 3: LA & HG, 4: From HA & LG to LA & HG (Vessel degradation), 5: From LA & HG to HA & LG (Tumor angiogenesis)
CTX_cond = 2;       % 1: CA-pos, 2: Neutral, 3: CA-neg
ECM_cond = 15;

T_nc1 = 25;         % Start time of CA field change
T_nc2 = 30;         % End time of CA field change

if CA_cond == 4
    CA_Change = 1;  % Condition of changing CA field
    NC_idx1 = 1;
    NC_idx2 = 3;
elseif CA_cond == 5
    CA_Change = 1;  % Condition of changing CA field
    NC_idx1 = 3;
    NC_idx2 = 1;
else
    CA_Change = 0;
end

ECM_den = ECM_cond * 0.1;
ECM_str = sprintf('%02d', ECM_cond);

ParamTag = [num2str(CA_cond), num2str(CTX_cond), num2str(ECM_str)];

%% Grid and Domain Setup
HW = 60;             % Halfwidth of the arena (Unit: millemeter)
sc = 0.1;            % Scale value (the smaller, the finer): the actual width of the grid
RepR = 0.5;          % Unit radius for suppressing transition to leader
SenseR = 0.3;        % Unit radius for ECM sensing

[x_grid, y_grid] = meshgrid((-HW:sc:HW), (-HW:sc:HW));

n_X = length(x_grid);
n_Y = length(y_grid);

BD_X = gpuArray.zeros(n_Y, n_X + 1);
BD_Y = gpuArray.zeros(n_Y + 1, n_X);

%% Time and Display Settings
dt = 0.005;          % Time step
T = 50;              % Total simulation time
SHWINT = 200;

VisMod = 1;
FileSav = 0;

TimeLapse = 0;       % Saving distribution in specific time points
LapsePoints = 10:10:100;

if TimeLapse == 1
    TL_idx = 1;
    TL_rec = cell(5, (1 + length(LapsePoints)));

    TL_rec{1,1} = ParamTag;
    TL_rec{2,1} = 'TotalMap';
    TL_rec{3,1} = 'LivMap';
    TL_rec{4,1} = 'Gyrus';
    TL_rec{5,1} = 'Sulcus';

    TL_idx = TL_idx + 1;
    MatName = [ParamTag, '_TimeLapse.mat'];
end

%% Initial Cell Distribution
Initial_total_size = 100;   % Unit: micro gram

d_data = sqrt((x_grid).^2 + y_grid.^2);  % Distance from the center

La = 6;
Lb = 6;

Y_boundary = Lb * sqrt(1 - (x_grid.^2) / (La^2));
InnerBin = (abs(y_grid) < Y_boundary);
P_num = imgaussfilt(double(InnerBin), 0.1/sc);

Tr_PC = 0.03;       % Transition rate from plain followers to leader-competent followers
Tr_CP = 0.06;       % Transition rate from leader-competent followers to plain followers

P_size0 = Initial_total_size / (1 + Tr_PC/Tr_CP);
C_size0 = Initial_total_size - P_size0;

P_num = P_size0 * P_num / (sum(P_num, 'all'));
C_num = C_size0 * P_num / (sum(P_num, 'all'));

P_density = P_num / (sc^2);   % Plain follower density
C_density = C_num / (sc^2);   % Leader-competent follower density

%% CA Field Setup
CA_field0 = sqrt(x_grid.^2 + y_grid.^2);

HALG_CA = WeakGrad(CA_field0);
LAHG_CA = StpGrad(CA_field0);
MAMG_CA = (HALG_CA + LAHG_CA) / 2;

[HALG_grad_x, HALG_grad_y] = BoundaryGrad(HALG_CA, sc, BD_X, BD_Y);
[MAMG_grad_x, MAMG_grad_y] = BoundaryGrad(MAMG_CA, sc, BD_X, BD_Y);
[LAHG_grad_x, LAHG_grad_y] = BoundaryGrad(LAHG_CA, sc, BD_X, BD_Y);

% --- CA field change setup ---
if CA_Change == 1
    CA_Field_cell = cell(3,1);
    CA_Grad_cell = cell(3,2);

    CA_Field_cell{1} = HALG_CA;
    CA_Field_cell{2} = MAMG_CA;
    CA_Field_cell{3} = LAHG_CA;

    CA_Grad_cell{1,1} = HALG_grad_x;   CA_Grad_cell{1,2} = HALG_grad_y;
    CA_Grad_cell{2,1} = MAMG_grad_x;   CA_Grad_cell{2,2} = MAMG_grad_y;
    CA_Grad_cell{3,1} = LAHG_grad_x;   CA_Grad_cell{3,2} = LAHG_grad_y;

    CA_Field_1 = cell(1,3);
    CA_Field_2 = cell(1,3);

    CA_Field_1{1} = CA_Grad_cell{NC_idx1, 1};
    CA_Field_1{2} = CA_Grad_cell{NC_idx1, 2};
    CA_Field_1{3} = CA_Field_cell{NC_idx1};

    CA_Field_2{1} = CA_Grad_cell{NC_idx2, 1};
    CA_Field_2{2} = CA_Grad_cell{NC_idx2, 2};
    CA_Field_2{3} = CA_Field_cell{NC_idx2};

    CA_Field_Tag = {'HA & LG', 'MA & MG', 'LA & HG'};
end

% --- CA field selection ---
if CA_cond == 1          % HA & LG
    CA_field = HALG_CA;
    CA_name = 'HALG';
    CA_title = 'HA & LG';
    CA_titleA = '';
    grad_CA_X = HALG_grad_x;
    grad_CA_Y = HALG_grad_y;
elseif CA_cond == 2      % MA & MG
    CA_field = MAMG_CA;
    CA_name = 'MAMG';
    CA_title = 'MA & MG';
    CA_titleA = '';
    grad_CA_X = MAMG_grad_x;
    grad_CA_Y = MAMG_grad_y;
elseif CA_cond == 3      % LA & HG
    CA_field = LAHG_CA;
    CA_name = 'LAHG';
    CA_title = 'LA & HG';
    CA_titleA = '';
    grad_CA_X = LAHG_grad_x;
    grad_CA_Y = LAHG_grad_y;
elseif CA_cond == 4      % Vessel degradation
    CA_field = HALG_CA;
    CA_name = 'AA';
    CA_titleA = ' (Vessel degradation)';
    grad_CA_X = HALG_grad_x;
    grad_CA_Y = HALG_grad_y;
elseif CA_cond == 5      % Tumor angiogenesis
    CA_field = LAHG_CA;
    CA_name = 'TA';
    CA_titleA = ' (Tumor angiogenesis)';
    grad_CA_X = LAHG_grad_x;
    grad_CA_Y = LAHG_grad_y;
end

%% ECM Setup
La2 = La * 1.5;
Lb2 = Lb * 1.5;

Y_boundary2 = Lb2 * sqrt(1 - (x_grid.^2) / (La2^2));
InnerBin2 = (abs(y_grid) < Y_boundary2);

E_density = zeros(size(P_density));        % Exhausted leaders

StiffDen = ones(size(d_data));
StiffDen(InnerBin2) = 0;
StiffDen = imgaussfilt(StiffDen, 1/sc);
StiffDen(InnerBin) = 0;
StiffDen = imgaussfilt(StiffDen, 0.2/sc);

StiffDen = StiffDen * ECM_den;
ECM_den_str = sprintf('%.1f', ECM_den);

ECM_name = ['ECM', ECM_den_str];
ECM_title = ['ECM = ', ECM_den_str];

%% Kinetic Parameters
A_density = zeros(size(d_data));   % Active leader density
D_density = zeros(size(d_data));   % Dead tumor cell density

% -- Chemotaxis coefficients --
kvP = 0.05;          % Chemotaxis of plain followers
kvC = 0.2;           % Chemotaxis of leader-competent followers
kvA = 2;             % Chemotaxis of active leaders

% -- ECM repulsion coefficients --
kSP = 1;             % ECM repulsion to plain followers
kSC = 1;             % ECM repulsion to leader-competent followers
kSA = 1;             % ECM repulsion to active leaders
kSE = 1;             % ECM repulsion to exhausted leaders

% -- Recovery and diffusion --
E_recov = 0.5;       % Exhausted leaders to become prospective leaders

Df_P = 0.25;         % Diffusion of plain followers
Df_C = 0.20;         % Diffusion of leader-competent followers
Df_A = 0.05;         % Diffusion of active leaders
Df_E = 0.01;         % Diffusion of exhausted leaders

% -- Death and decay rates --
P_death_rate = 0.03;  % Death rate of plain follower
C_death_rate = 0.03;  % Death rate of leader-competent follower
E_death_rate = 0.5;   % Death rate of exhausted leader
D_lysis_rate = 0.02;  % Lysis rate of dead tumor cell

% -- Growth rate --
P_gr = 0.05;          % Basal plain follower growth rate

%% Leader Suppression Filter
r1 = 3 * RepR / sc;
r2 = 10 * RepR / sc;  % Inner and outer radii for suppression of new active leader emergence

A_sp = 2;              % Active leader's suppression to prevent emergence of new active leaders
decayA = 100;         % Basal decay rate of active leaders

DG_P = 200;           % Proteolysis material's (MMP) ECM degradation rate
DG_L = 0.05;          % Lactate's ECM degradation rate

H = fspecial('disk', SenseR/sc);   % Disk filter (average) for stiffness

filter_size = 20 / sc;
[x, y] = meshgrid(linspace(-filter_size/2, filter_size/2, filter_size));

G_filtered0 = ((sqrt(x.^2 + y.^2) > r1) & (sqrt(x.^2 + y.^2) < r2));
G_filtered0 = double(G_filtered0);
G_filtered1 = imfilter(G_filtered0, H);

[kn_Y, kn_X] = size(G_filtered1);

pad_Y = n_Y + kn_Y - 1;
pad_X = n_X + kn_X - 1;

fft_kernel = fft2(G_filtered1, pad_Y, pad_X);
fft_kernel = gpuArray(fft_kernel);   % Kernel for leader suppression cue

%% Leader Transition Thresholds
TransStiffMin = 0.1;
TransStiffMax = 0.3;

%% Chemotaxis Condition (CTX)
switch CTX_cond
    % 1: CA-pos
    % 2: Neutral
    % 3: CA-neg
    case 1
        C = 40;
        [CA_W_x, CA_W_y] = BoundaryDen(CA_Weight(CA_field, C), BD_X, BD_Y);

        QC_x = CA_W_x .* grad_CA_X;
        QC_y = CA_W_y .* grad_CA_Y;

        CTX_name = 'CA-PosCTX';
        CTX_title = 'CA-positive chemotaxis';
    case 2
        QC_x = grad_CA_X;
        QC_y = grad_CA_Y;

        CTX_name = 'NeuCTX';
        CTX_title = 'Neutral chemotaxis';
    case 3
        C = 20;
        [CA_W_x, CA_W_y] = BoundaryDen(2.5 - CA_Weight(CA_field, C), BD_X, BD_Y);

        QC_x = CA_W_x .* grad_CA_X;
        QC_y = CA_W_y .* grad_CA_Y;

        CTX_name = 'CA-NegCTX';
        CTX_title = 'CA-negative chemotaxis';
end

%% GPU Transfer
FigCount = 0;

CA_field_G = gpuArray(CA_field);
A_density_G = gpuArray(A_density);
P_density_G = gpuArray(P_density);
C_density_G = gpuArray(C_density);
StiffDen_G = gpuArray(StiffDen);
E_density_G = gpuArray(E_density);
D_density_G = gpuArray(D_density);

x_grid = gpuArray(x_grid);
y_grid = gpuArray(y_grid);

%% Duration Recording (Initial)
if DurRec == 1
    dur_ct = 1;
    Rec_dT = 0.5;

    DurRecMat = zeros(6, round(T/Rec_dT) + 1);
    % Row 1: time, Row 2: total mass, Row 3: live mass,
    % Row 4: invasion index, Row 5: area, Row 6: starness

    StarShapeAnalyze

    DurRecMat(2,dur_ct) = TumMass_tot;
    DurRecMat(3,dur_ct) = TumMass_liv;
    DurRecMat(4,dur_ct) = Inv_score;
    DurRecMat(5,dur_ct) = TumArea;
    DurRecMat(6,dur_ct) = PeakScore;
end

%% Figure Initialization
if VisMod == 1
    fig = figure('Name','Progress','Visible','on', 'position', [400 200 1600 1000], 'color', 'white');

    if CA_Change == 1
        CA_title = [CA_Field_Tag{NC_idx1}, ' to ', CA_Field_Tag{NC_idx2}, CA_titleA];
    end

    Fig_title = [CA_title, ', ', CTX_title, ', ', ECM_title];
    f_name = [ParamTag, ' ', CA_name, ' ', CTX_name, ' ', ECM_name, ' (M21)'];

    minUran = 0.00001;
    minUranD = 0.01;

    tiledlayout(2,4,'Padding','compact','TileSpacing','compact');

    % -- Tile 1: Plain followers --
    h1 = nexttile;
    im1 = imagesc(zeros(size(A_density))); axis image off;
    title('Plain followers'); colormap(im1.Parent,'pink');
    cb1 = colorbar('southoutside');

    % -- Tile 2: Leader-competent followers --
    h2 = nexttile;
    im2 = imagesc(zeros(size(A_density))); axis image off;
    title('Leader-competent followers'); colormap(im2.Parent,'pink');
    cb2 = colorbar('southoutside');

    % -- Tile 3: ECM --
    h3 = nexttile;
    im3 = imagesc(zeros(size(A_density))); axis image off;
    title('ECM'); colormap(im3.Parent,'pink');
    cb3 = colorbar('southoutside'); box on
    hold on
    xLimits = xlim; yLimits = ylim;
    rectangle('Position', [xLimits(1), yLimits(1), diff(xLimits), diff(yLimits)], 'EdgeColor', 'k', 'LineWidth', 0.5);
    hold off

    % -- Tile 4: Chemoattractant --
    hN = nexttile;
    imN = imagesc(CA_field); axis image off;
    title('Chemoattractant'); colormap(imN.Parent,'pink');
    cbN = colorbar('southoutside');
    hold on
    xLimits = xlim; yLimits = ylim;
    rectangle('Position', [xLimits(1), yLimits(1), diff(xLimits), diff(yLimits)], 'EdgeColor', 'k', 'LineWidth', 0.5);
    hold off

    % -- Tile 5: Active leaders (fixed scale) --
    h4 = nexttile;
    im4 = imagesc(zeros(size(A_density))); axis image off;
    title('Active leaders (fixed scale)'); colormap(im4.Parent,'pink');
    cb4 = colorbar('southoutside');

    % -- Tile 6: Active leaders (liberal scale) --
    h5 = nexttile;
    im5 = imagesc(zeros(size(A_density))); axis image off;
    title('Active leaders (liberal scale)'); colormap(im5.Parent,'pink');
    cb5 = colorbar('southoutside');

    % -- Tile 7: Leader transition suppression --
    h6 = nexttile;
    im6 = imagesc(zeros(size(A_density))); axis image off;
    title('Leader transition suppression'); colormap(im6.Parent,'pink');
    cb6 = colorbar('southoutside');

    % -- Tile 8: Dead cells --
    hE = nexttile;
    imE = imagesc(zeros(size(A_density))); axis image off;
    title('Dead cells'); colormap(imE.Parent,'pink');
    cbE = colorbar('southoutside');

    % -- Color axis limits --
    hN.CLim = [0 70];         % CA
    h3.CLim = [0 3];          % ECM density
    h4.CLim = [0 0.00001];    % Active leader (fixed)
    h5.CLim = [0 minUran];
    h6.CLim = [0 minUran];
    hE.CLim = [0 minUranD];

    % -- Colorbar exponent formatting --
    cb3.Ruler.Exponent = 0;
    cb4.Ruler.Exponent = 0;
    cb5.Ruler.Exponent = 0;
    cb6.Ruler.Exponent = 0;
    cbN.Ruler.Exponent = 0;
    cbE.Ruler.Exponent = 0;

    % -- Video output --
    if FileSav
        Vdata = VideoWriter(f_name,'MPEG-4');
        Vdata.FrameRate = 2;
        Vdata.Quality = 100;
        open(Vdata)
    end
end

%% Running Numerical Computation
for t = 0:dt:T

    % --- CA field transition ---
    if (CA_Change == 1) && (t >= T_nc1) && (t <= T_nc2)
        [QC_x, QC_y, CA_field_G] = CA_Trans(CA_Field_1, CA_Field_2, t, T_nc1, T_nc2, CTX_cond, BD_X, BD_Y);
    end

    % --- Derived fields ---
    CA_Rep = CA_field_G - 20;          % CA for the cell replication
    CA_Rep = max(CA_Rep, 0);

    CA_field_net = min(CA_field_G, 40);  % CA associated with active leader decay

    A_supp_cue = radial_filter_fft(A_density_G, fft_kernel, kn_Y, kn_X);
    A_supp_cue(A_supp_cue < 1.0e-10) = 0;

    Proteo = imgaussfilt(A_density_G, 0.7/sc, 'FilterSize', [71 71]);             % Proteolysis material (MMP) density
    Lactate = imgaussfilt((P_density_G + C_density_G), 0.4/sc, 'FilterSize', [51 51]);  % Lactate density (weakly degrading ECM)

    StiffSense = imfilter(StiffDen_G, H);

    TotalDen = A_density_G + C_density_G + P_density_G + E_density_G + D_density_G;

    % --- Spatial gradients ---
    [M_x, M_y] = BoundaryGrad(StiffSense, sc, BD_X, BD_Y);

    [A_x, A_y] = BoundaryGrad(A_density_G, sc, BD_X, BD_Y);
    [P_x, P_y] = BoundaryGrad(P_density_G, sc, BD_X, BD_Y);
    [C_x, C_y] = BoundaryGrad(C_density_G, sc, BD_X, BD_Y);
    [E_x, E_y] = BoundaryGrad(E_density_G, sc, BD_X, BD_Y);

    % --- Boundary densities ---
    [A_BD_den_X, A_BD_den_Y] = BoundaryDen(A_density_G, BD_X, BD_Y);
    [P_BD_den_X, P_BD_den_Y] = BoundaryDen(P_density_G, BD_X, BD_Y);
    [C_BD_den_X, C_BD_den_Y] = BoundaryDen(C_density_G, BD_X, BD_Y);
    [E_BD_den_X, E_BD_den_Y] = BoundaryDen(E_density_G, BD_X, BD_Y);
    [D_BD_den_X, D_BD_den_Y] = BoundaryDen(D_density_G, BD_X, BD_Y);

    T_BD_den_X = A_BD_den_X + P_BD_den_X + C_BD_den_X + E_BD_den_X + D_BD_den_X;
    T_BD_den_Y = A_BD_den_Y + P_BD_den_Y + C_BD_den_Y + E_BD_den_Y + D_BD_den_Y;

    [M_BD_den_X, M_BD_den_Y] = BoundaryDen(StiffDen_G, BD_X, BD_Y);

    % --- Movement impedance ---
    MoveImpIdx_BD_X = (0.5 * (T_BD_den_X) + 5 * M_BD_den_X);
    MoveImpIdx_BD_Y = (0.5 * (T_BD_den_Y) + 5 * M_BD_den_Y);

    ImpMovA_X = 2 ./ (1 + exp(1.5 * MoveImpIdx_BD_X));  % Impeded movement of active leaders
    ImpMovA_Y = 2 ./ (1 + exp(1.5 * MoveImpIdx_BD_Y));

    ImpMovP_X = 2 ./ (1 + exp(3 * MoveImpIdx_BD_X));
    ImpMovP_Y = 2 ./ (1 + exp(3 * MoveImpIdx_BD_Y));

    ImpMovE_X = ImpMovP_X;    ImpMovE_Y = ImpMovP_Y;
    ImpMovC_X = ImpMovA_X;    ImpMovC_Y = ImpMovA_Y;

    % --- Flux computation ---
    J_Ax = -ImpMovA_X .* (A_BD_den_X .* (-kvA * QC_x + kSA * M_x) + Df_A * A_x);  % X-flux of active leaders
    J_Ay = -ImpMovA_Y .* (A_BD_den_Y .* (-kvA * QC_y + kSA * M_y) + Df_A * A_y);  % Y-flux of active leaders

    J_Px = -ImpMovP_X .* (P_BD_den_X .* (-kvP * QC_x + kSP * M_x) + Df_P * P_x);
    J_Py = -ImpMovP_Y .* (P_BD_den_Y .* (-kvP * QC_y + kSP * M_y) + Df_P * P_y);

    J_Cx = -ImpMovC_X .* (C_BD_den_X .* (-kvC * QC_x + kSC * M_x) + Df_C * C_x);
    J_Cy = -ImpMovC_Y .* (C_BD_den_Y .* (-kvC * QC_y + kSC * M_y) + Df_C * C_y);

    J_Ex = -ImpMovE_X .* (E_BD_den_X .* (kSE * M_x) + Df_E * E_x);
    J_Ey = -ImpMovE_Y .* (E_BD_den_Y .* (kSE * M_y) + Df_E * E_y);

    CTX_size = sqrt(Collapse_X(QC_x).^2 + Collapse_Y(QC_y).^2);

    P_mov_cost = 0.1 * kvP * CTX_size;
    C_mov_cost = 0.05 * kvC * CTX_size;

    thr1 = 1e-10;

    J_Ax_tot = J_Ax;    J_Ay_tot = J_Ay;
    J_Px_tot = J_Px;    J_Py_tot = J_Py;
    J_Cx_tot = J_Cx;    J_Cy_tot = J_Cy;
    J_Ex_tot = J_Ex;    J_Ey_tot = J_Ey;

    A_div = Flux_Div(J_Ax_tot, J_Ay_tot, sc);
    P_div = Flux_Div(J_Px_tot, J_Py_tot, sc);
    C_div = Flux_Div(J_Cx_tot, J_Cy_tot, sc);
    E_div = Flux_Div(J_Ex_tot, J_Ey_tot, sc);

    A_density_Ga = A_density_G - dt * A_div;
    P_density_Ga = P_density_G - dt * P_div;
    C_density_Ga = C_density_G - dt * C_div;
    E_density_Ga = E_density_G - dt * E_div;

    % --- Reaction step ---
    DenEffect = (1 - TotalDen / 1.5);  % Effect of the maximal confluence

    P_density_G = P_density_Ga + dt * ((P_gr * CA_Rep .* DenEffect - P_death_rate - P_mov_cost) .* P_density_Ga) ...
        + C_density_Ga * Tr_CP * dt - P_density_Ga * Tr_PC * dt;

    StiffDenTrans = (StiffDen_G - TransStiffMin);
    StiffDenTrans = max(StiffDenTrans, 0);
    StiffDenTrans(StiffDenTrans > (TransStiffMax - TransStiffMin)) = (TransStiffMax - TransStiffMin);
    StiffDenTrans = StiffDenTrans * (1 / (TransStiffMax - TransStiffMin));

    A_sus = decayA * (StiffSense + 0.1) ./ (CA_field_net + 10);  % Lifetime of the active leaders

    A_density_G = A_density_Ga + dt * (max(C_density_Ga .* StiffDenTrans - A_sp * A_supp_cue, 0) - A_sus .* A_density_Ga);
    A_density_G = max(A_density_G, 0);

    C_density_G = C_density_Ga + dt * (-max(C_density_Ga .* StiffDenTrans - A_sp * A_supp_cue, 0)) + ...
        + P_density_Ga * Tr_PC * dt - C_density_Ga .* (Tr_CP + C_death_rate + C_mov_cost) * dt + (E_recov * CA_field_net .* E_density_Ga) * dt;

    E_density_G = E_density_Ga + dt * (A_sus .* A_density_Ga - (E_recov * CA_field_net + E_death_rate) .* E_density_Ga);

    D_density_G = D_density_G + (P_death_rate + P_mov_cost) .* P_density_Ga * dt + ...
        (C_death_rate + C_mov_cost) .* C_density_Ga * dt + ...
        E_death_rate .* E_density_Ga * dt - D_lysis_rate * D_density_G * dt;

    StiffDen_G = StiffDen_G + dt * (-DG_P * Proteo - DG_L * Lactate);
    StiffDen_G = max(StiffDen_G, 0);

    %% Duration Recording (In-loop)
    if DurRec == 1 && rem(t, Rec_dT) == 0

        dur_ct = dur_ct + 1;
        % Row 1: time, Row 2: total mass, Row 3: live mass,
        % Row 4: invasion index, Row 5: area, Row 6: starness

        StarShapeAnalyze

        DurRecMat(1,dur_ct) = t;
        DurRecMat(2,dur_ct) = TumMass_tot;
        DurRecMat(3,dur_ct) = TumMass_liv;
        DurRecMat(4,dur_ct) = Inv_score;
        DurRecMat(5,dur_ct) = TumArea;
        DurRecMat(6,dur_ct) = PeakScore;

        if PeakScore > 100
            pause
        end

        if TimeLapse == 1 && ismember(t, LapsePoints)
            TL_rec{1, TL_idx} = t;
            TL_rec{2, TL_idx} = densityMat_tot_C0;
            TL_rec{3, TL_idx} = densityMat_Liv_C0;
            TL_rec{4, TL_idx} = Gyrus_data;
            TL_rec{5, TL_idx} = Sulcus_data;

            TL_idx = TL_idx + 1;
        end
    end

    %% Visualization
    if VisMod == 1

        if rem(FigCount, SHWINT) == 0

            T = gather(cat(3, P_density_G, C_density_G, A_density_G, E_density_G, StiffDen_G, A_supp_cue, D_density_G));

            set(im1, 'CData', T(:,:,1));
            set(im2, 'CData', T(:,:,2));
            set(im3, 'CData', T(:,:,5));
            set(im4, 'CData', T(:,:,3));
            set(im5, 'CData', T(:,:,3));
            set(im6, 'CData', T(:,:,6));
            set(imE, 'CData', T(:,:,7));

            sgtitle(fig, sprintf('%s,\\it t\\rm = %.2f', Fig_title, t));

            if max(A_density_G, [], 'all') < minUran
                h5.CLim = [0 minUran];
            else
                h5.CLim = [0 inf];
            end

            if max(A_supp_cue, [], 'all') < minUran
                h6.CLim = [0 minUran];
            else
                h6.CLim = [0 inf];
            end

            if max(D_density_G, [], 'all') < minUranD
                hE.CLim = [0 minUranD];
            else
                hE.CLim = [0 inf];
            end

            if CA_Change == 1
                Tn = gather(CA_field_G);
                set(imN, 'CData', Tn);
            end

            drawnow
            pause(0.01)

            if FileSav
                frameV = getframe(fig);
                writeVideo(Vdata,frameV)
            end
        end

        FigCount = FigCount + 1;
    end
end

%% Cleanup
if FileSav
    close(Vdata)
end
