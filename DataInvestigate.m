% DataInvestigate.m
%
% Data file: RawData.mat  (variable: DataCell, 282 x 4 cell)
%
% Each row of DataCell is one simulated condition:
%   {i,1}  ParamTag     - 4-character condition code, 'ABCD':
%                           A    = chemoattractant (CA) configuration
%                           B    = chemotaxis (CTX) strategy
%                           C,D  = ECM density index (ECM density = 0.1 * index)
%   {i,2}  DurRecMat    - 6 x 101 time record at t = 0, 0.5, ..., 50; rows:
%                           1 = time, 2 = total tumor mass, 3 = live tumor mass,
%                           4 = invasion score, 5 = tumor area, 6 = peak score
%   {i,3}  Gyrus_data   - 2 x nGyri, rows: 1 = angle theta, 2 = radius
%   {i,4}  Sulcus_data  - 2 x nSulci, rows: 1 = angle theta, 2 = radius
%
%
% CA configuration (ParamTag character 1)
%   1 = high average, mild gradient (HA & LG)
%   2 = medium average, medium gradient (MA & MG)
%   3 = low average, steep gradient (LA & HG)
%   4 = HA & LG -> LA & HG over time (anti-angiogenesis)
%   5 = LA & HG -> HA & LG over time (rapid angiogenesis)
%
% CTX strategy (ParamTag character 2)
%   1 = CA-positive (higher CA enhances chemotaxis)
%   2 = normal (gradient-based)
%   3 = CA-negative (higher CA reduces chemotaxis)
%
% Marker convention used throughout:
%   hue         = CTX strategy   (blue / green / red)
%   shape       = CA config      (square / circle / triangle)
%   edge shade  = ECM density    (darker = denser)

%% ------------------------------------------------------------------ Settings

DataFile = 'RawData.mat';

ViewMode = 3;
% 1: Result vs. ECM density   (fixed CA and CTX, final time)
% 2: Heatmap over ECM and time (fixed CA and CTX)
% 3: Result vs. time          (fixed CA, CTX, and ECM)
% 4: Result vs. CA config     (fixed ECM)

CmpRes = 2;
% Which quantity to plot:
%   1 = total tumor mass, 2 = live tumor mass, 3 = invasion score,
%   4 = tumor area,       5 = peak score

unit_conv = 10^(-5);   % invasion score: ug*mm^2 -> g*cm^2

% Row of DurRecMat holding the selected quantity (row 1 is time).
RecRow = CmpRes + 1;

% ViewMode 1: averaging window, in simulation time units. Editable; the default
% (50 to 50) averages over the final time point alone.
t_init = 50;
t_end  = 50;

% ViewMode 4: ECM density index to hold fixed. Editable; 1..30 = density
% 0.1..3.0. The nine static CA (1..3) x CTX (1..3) combinations are plotted here.
tg_ecm = 10;

Data_labels = {'Total tumor mass', 'Live tumor mass', 'Invasion score', ...
               'Tumor area', 'Peak score'};

%% ---------------------------------------------------------------- Load data

if ~exist('DataCell', 'var')
    if ~exist(DataFile, 'file')
        error('Cannot find %s. Run this script from the folder containing it.', DataFile);
    end
    load(DataFile, 'DataCell');
end

[nr, ~] = size(DataCell);

%% ------------------------------------------- Condition / result index tables

% All nine CA x CTX combinations, in (CA, CTX) order.
tg_cds_T = zeros(9, 2);
cnt0 = 1;
for gi = 1:3
    for gj = 1:3
        tg_cds_T(cnt0, 1) = gi;
        tg_cds_T(cnt0, 2) = gj;
        cnt0 = cnt0 + 1;
    end
end

% tg_cds selects WHICH conditions are drawn, and is meant to be edited: the
% values below reproduce the published panels, but any valid combination can be
% substituted to inspect a different slice of the 282 deposited conditions.
% (ViewMode 4 ignores tg_cds; it is driven by tg_ecm in the Settings block.)
%
% Each ROW of tg_cds is one condition group, and its COLUMNS are:
%   ViewMode 1: [CA, CTX]       - ECM density is swept automatically; one curve
%                                 per row.
%   ViewMode 2: [CA, CTX]       - exactly one row; ECM and time span the mesh.
%   ViewMode 3: [CA, CTX, ECM]  - one curve per row.
%
% Valid values (see the header for what each code means):
%   CA  = 1..5, but 4 and 5 (the time-varying angiogenesis cases) were only run
%         at ECM index 10 and 20. Use CA = 1..3 for ECM sweeps (ViewMode 1, 2).
%   CTX = 1..3.
%   ECM = 1..30, i.e. ECM density 0.1..3.0 in steps of 0.1.
%
% tg_cds_T lists all nine CA x CTX combinations in (CA, CTX) order, so rows 4:6
% are CA = 2 paired with each CTX in turn. Passing the whole of tg_cds_T to
% ViewMode 1 plots all nine combinations at once.

switch ViewMode
    case 1
        tg_cds = tg_cds_T(4:6, :);          % CA = 2, all three CTX strategies
    case 2
        tg_cds = [2, 2];                    % CA = 2, CTX = 2
    case 3
        tg_cds = [3, 1, 20; ...             % columns: CA, CTX, ECM index
                  3, 2, 20; ...
                  3, 3, 20];
end

%% ------------------------------------------------- Summary table of all runs

% WholeResults columns:
%   1 = CA config, 2 = CTX strategy, 3 = ECM index
%   4 = total tumor mass, 5 = live tumor mass, 6 = invasion score,
%   7 = tumor area, 8 = peak score  (all at the final recorded time)
%   9 = number of peaks (gyri)
%  10 = mean sulcus (valley) radius, 11 = mean gyrus (peak) radius

WholeResults = zeros(nr, 11);

for dct = 1:nr

    Condition_str = DataCell{dct, 1};

    WholeResults(dct, 1) = str2double(Condition_str(1));    % CA config
    WholeResults(dct, 2) = str2double(Condition_str(2));    % CTX strategy
    WholeResults(dct, 3) = str2double(Condition_str(3:4));  % ECM index

    Record_MAT = DataCell{dct, 2};
    Final_rec  = Record_MAT(:, end);

    Gyrus_data  = DataCell{dct, 3};
    Sulcus_data = DataCell{dct, 4};

    Gyrus_r  = Gyrus_data(2, :);
    Sulcus_r = Sulcus_data(2, :);

    WholeResults(dct, 4:8) = Final_rec(2:6)';

    WholeResults(dct, 9)  = numel(Gyrus_r);
    WholeResults(dct, 10) = mean(Sulcus_r);
    WholeResults(dct, 11) = mean(Gyrus_r);

end

Cds = WholeResults(:, 1:3);   % conditions only

%% ------------------------------------------------------- Plot style per run

CLR        = zeros(size(Cds, 1), 3);   % marker face colour
MarkerList = cell(size(Cds, 1), 1);    % marker shape
MarkerSize = zeros(size(Cds, 1), 1);   % marker size
BD_List    = cell(size(Cds, 1), 1);    % marker edge colour

bd_str = 0.4;   % line / edge width

sc_ftr = 0.5;
SqMkSize = 6.5 * sc_ftr;
CiMkSize = 5   * sc_ftr;
TrMkSize = 5   * sc_ftr;

dk_colors = [  0   0 255;    % blue  - CA-positive
               0 250   0;    % green - gradient-based
             200   0  20] / 255;   % red - CA-negative

for i = 1:size(Cds, 1)

    hue_idx   = Cds(i, 2);   % CTX strategy -> hue
    shape_idx = Cds(i, 1);   % CA config    -> shape

    CLR(i, :) = dk_colors(hue_idx, :);

    switch shape_idx
        case 1                                  % HA & LG
            MarkerList{i} = "square";  MarkerSize(i) = SqMkSize;
        case 2                                  % MA & MG
            MarkerList{i} = "o";       MarkerSize(i) = CiMkSize;
        case 3                                  % LA & HG
            MarkerList{i} = "^";       MarkerSize(i) = TrMkSize;
        case 4                                  % anti-angiogenesis (ends LA & HG)
            MarkerList{i} = "^";       MarkerSize(i) = TrMkSize;
        case 5                                  % rapid angiogenesis (ends HA & LG)
            MarkerList{i} = "square";  MarkerSize(i) = SqMkSize;
    end

    % ECM density -> edge shade (denser ECM = darker edge)
    BD_List{i} = 0.8 * [1 1 1] * (1 - Cds(i, 3) * 0.1 / 3);

end

%% ------------------------------------------------------------------- Figures

switch ViewMode

    case 1   % ---------------------------------- Result vs. ECM density

        figure('Units', 'pixels', 'Position', [400 600 330 200]);

        for cnum = 1:size(tg_cds, 1)

            tg_CA  = tg_cds(cnum, 1);
            tg_CTX = tg_cds(cnum, 2);

            tg_pos = find(Cds(:, 1) == tg_CA & Cds(:, 2) == tg_CTX);

            ECM_range = 0.1 * WholeResults(tg_pos, 3);

            Avr_res = zeros(numel(tg_pos), 1);

            for j0 = 1:numel(tg_pos)

                REC_DATA = DataCell{tg_pos(j0), 2};
                t_vals   = REC_DATA(1, :);

                t_trg = find(t_vals >= t_init & t_vals <= t_end);

                Avr_res(j0) = sum(REC_DATA(RecRow, t_trg)) / numel(t_trg);

            end

            if CmpRes == 3
                Avr_res = Avr_res * unit_conv;
            end

            plot(ECM_range, Avr_res, 'Color', CLR(tg_pos(1), :), 'LineWidth', bd_str)
            hold on

            % Mark only round ECM values, to keep the line readable.
            for w0 = 1:numel(Avr_res)
                idx_val = tg_pos(w0);
                if ismember(ECM_range(w0), 0.5:0.5:3)
                    plot(ECM_range(w0), Avr_res(w0), 'Marker', MarkerList{idx_val}, ...
                        'MarkerFaceColor', CLR(idx_val, :), 'MarkerSize', MarkerSize(idx_val), ...
                        'MarkerEdgeColor', BD_List{idx_val}, 'LineWidth', bd_str)
                end
            end

        end

        xlabel('ECM density')
        ylabel(Data_labels{CmpRes})
        box off

    case 2   % ---------------------------------- Heatmap over ECM and time

        figure('Units', 'pixels', 'Position', [400 600 330 250]);

        tg_CA  = tg_cds(1);
        tg_CTX = tg_cds(2);

        tg_pos = find(Cds(:, 1) == tg_CA & Cds(:, 2) == tg_CTX);

        ECM_range = 0.1 * WholeResults(tg_pos, 3);

        HM_mat = [];

        for j0 = 1:numel(tg_pos)
            REC_DATA = DataCell{tg_pos(j0), 2};
            t_vals   = REC_DATA(1, :);
            HM_mat   = [HM_mat; REC_DATA(RecRow, :)];
        end

        if CmpRes == 3
            HM_mat = HM_mat * unit_conv;
        end

        [Xm, Ym] = meshgrid(t_vals, ECM_range);
        mesh(Xm, Ym, HM_mat)

        xlabel('Time')
        ylabel('ECM density')
        zlabel(Data_labels{CmpRes})

        if CmpRes == 3
            set(gca, 'ZScale', 'log')
        end

    case 3   % ---------------------------------- Result vs. time

        figure('Units', 'pixels', 'Position', [750 700 260 170]);

        for cnum = 1:size(tg_cds, 1)

            tg_CA  = tg_cds(cnum, 1);
            tg_CTX = tg_cds(cnum, 2);
            tg_ECM = tg_cds(cnum, 3);

            tg_pos = find(Cds(:, 1) == tg_CA & Cds(:, 2) == tg_CTX & Cds(:, 3) == tg_ECM);

            for j0 = 1:numel(tg_pos)

                idx_val = tg_pos(j0);

                REC_DATA  = DataCell{idx_val, 2};
                t_vals    = REC_DATA(1, :);
                data_vals = REC_DATA(RecRow, :);

                if CmpRes == 3
                    data_vals = data_vals * unit_conv;
                end

                plot(t_vals, data_vals, 'LineWidth', bd_str, 'Color', CLR(idx_val, :))
                hold on

                if ismember(tg_CA, [1 2 3])

                    % Static CA config: one marker shape throughout.
                    sel_t_idx = find(ismember(t_vals, 5:5:50));

                    plot(t_vals(sel_t_idx), data_vals(sel_t_idx), 'Marker', MarkerList{idx_val}, ...
                        'MarkerFaceColor', CLR(idx_val, :), 'MarkerSize', MarkerSize(idx_val), ...
                        'MarkerEdgeColor', BD_List{idx_val}, 'LineWidth', bd_str, 'Color', 'none')

                else

                    % Time-varying CA config: marker switches shape at the
                    % crossover, so early and late regimes are distinguishable.
                    if tg_CA == 4        % HA & LG -> LA & HG
                        MK_A = 'square'; MS_A = SqMkSize;
                        MK_B = '^';      MS_B = TrMkSize;
                    else                 % tg_CA == 5: LA & HG -> HA & LG
                        MK_A = '^';      MS_A = TrMkSize;
                        MK_B = 'square'; MS_B = SqMkSize;
                    end

                    sel_t_idx_A = find(ismember(t_vals,  5:5:25));
                    sel_t_idx_B = find(ismember(t_vals, 30:5:50));

                    plot(t_vals(sel_t_idx_A), data_vals(sel_t_idx_A), 'Marker', MK_A, ...
                        'MarkerFaceColor', CLR(idx_val, :), 'MarkerSize', MS_A, ...
                        'MarkerEdgeColor', BD_List{idx_val}, 'LineWidth', bd_str, 'Color', 'none')

                    plot(t_vals(sel_t_idx_B), data_vals(sel_t_idx_B), 'Marker', MK_B, ...
                        'MarkerFaceColor', CLR(idx_val, :), 'MarkerSize', MS_B, ...
                        'MarkerEdgeColor', BD_List{idx_val}, 'LineWidth', bd_str, 'Color', 'none')

                end

            end

        end

        xlabel('Time')
        ylabel(Data_labels{CmpRes})

        ax = gca;
        ax.TickLength = [0.02 0.04];
        xticks(0:10:100)
        box off

        if CmpRes == 3
            set(ax, 'YScale', 'log')
            ylim([5*10^(-3) 5])
            yticks([0.01, 0.1, 1])
        end

    case 4   % ---------------------------------- Result vs. CA config

        figure('Units', 'pixels', 'Position', [750 700 150 230]);

        tg_pos = find(Cds(:, 3) == tg_ecm);

        ResVals     = zeros(3, 3);   % (CTX, CA) -> value at final time
        ResVals_idx = zeros(3, 3);   % (CTX, CA) -> row of DataCell

        for k0 = 1:9

            pos_val = tg_pos(k0);

            tg_CA  = Cds(pos_val, 1);
            tg_CTX = Cds(pos_val, 2);

            tg_Cell = DataCell{pos_val, 2};

            ResVals(tg_CTX, tg_CA)     = tg_Cell(RecRow, end);
            ResVals_idx(tg_CTX, tg_CA) = pos_val;

        end

        if CmpRes == 3
            ResVals = ResVals * unit_conv;
        end

        % One line per CTX strategy, nudged sideways so they do not overlap.
        for k0 = 1:3
            idx_val = ResVals_idx(k0, 1);
            devX    = 0.08 * (2 - Cds(idx_val, 2));
            plot((1:3) + devX, ResVals(k0, :), 'Color', CLR(idx_val, :), 'LineWidth', bd_str)
            hold on
        end

        for k0 = 1:9

            idx_val = tg_pos(k0);

            tg_CA  = Cds(idx_val, 1);
            tg_CTX = Cds(idx_val, 2);
            devX   = 0.08 * (2 - tg_CTX);

            plot(tg_CA + devX, ResVals(tg_CTX, tg_CA), 'Marker', MarkerList{idx_val}, ...
                'MarkerFaceColor', CLR(idx_val, :), 'MarkerSize', MarkerSize(idx_val), ...
                'MarkerEdgeColor', BD_List{idx_val}, 'LineWidth', bd_str)

        end

        xlim([0.6, 3.4])
        xticks(1:3)
        xticklabels({'HA & LG', 'MA & MG', 'LA & HG'})
        ylabel(Data_labels{CmpRes})
        box off

        if CmpRes == 3
            set(gca, 'YScale', 'log')
            ylim([5*10^(-3) 5])
        end

        if CmpRes == 2
            ylim([0 700])
            yticks(0:200:700)
        end

        ax = gca;
        ax.TickLength = [0.02 0.05];

end
