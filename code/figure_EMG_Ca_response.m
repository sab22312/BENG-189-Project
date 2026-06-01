% Plot the EMG-to-calcium activation response used by the PAD/Ca extension.
%
% Run from the code folder:
%   cd code
%   figure_EMG_Ca_response

clear variables
clear; close all;

outdir = 'figure_EMG_Ca_response';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

data_resting = readtable('../raw_data/Initial_state.xlsx', 'Sheet', 'Summary');
table_emg = readtable('../raw_data/Emg_for_fitting_DPF.xlsx');
x_emg = table_emg{:,1};
y_emg = table_emg{:,2};
[hill_a, hill_b, hill_c] = classic_hill(x_emg, y_emg);

param_table = readtable('params/params.xlsx');
params = param_table.estimate;

data_Pcr = readtable('../raw_data/Pcr_for_fitting_DPF.xlsx');
cycle_index_exp = data_Pcr{:,1};
cycles = 1:1:max(cycle_index_exp);
iemg_profile = ((hill_a*(cycles.^hill_b))./(cycles + hill_c))/100;

pad = defaultCalciumConfig();
response = simulateCalciumResponse(params, cycles, iemg_profile, data_resting, pad);

fig = figure(1);
clf;
set(fig, 'Color', 'w', 'InvertHardcopy', 'off', 'Units', 'Inches', 'Position', [1 1 8.25 3.75]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
yyaxis left;
plot(x_emg, y_emg, '.', 'Color', [0.55 0.55 0.55], 'MarkerSize', 10); hold on;
plot(cycles, iemg_profile*100, 'k', 'LineWidth', 2);
ylabel('Normalized EMG (%)');
ylim([0 max(iemg_profile*100)*1.15]);

yyaxis right;
plot(cycles, response.Activation, 'Color', [0.10 0.35 0.75], 'LineWidth', 2);
ylabel('Ca activation (0-1)');
ylim([0 1]);

xlabel('Cycle Index');
title('EMG input drives Ca activation');
xlim([0 max(cycles)]);
styleCurrentAxes();

nexttile;
scatter(iemg_profile*100, response.Activation, 22, cycles, 'filled'); hold on;
plot(iemg_profile*100, response.Activation, 'Color', [0.10 0.35 0.75], 'LineWidth', 1.5);
xlabel('Normalized EMG (%)');
ylabel('Ca activation (0-1)');
title('EMG vs Ca activation response');
ylim([0 1]);
cb = colorbar;
cb.Label.String = 'Cycle Index';
styleCurrentAxes();

exportgraphics(fig, fullfile(outdir, 'EMG_Ca_activation_response.png'), 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, fullfile(outdir, 'EMG_Ca_activation_response.pdf'), 'BackgroundColor', 'w', 'ContentType', 'vector');

summary = table(cycles(:), iemg_profile(:)*100, response.Ca(:), response.Activation(:), ...
    'VariableNames', {'CycleIndex', 'NormalizedEMG_percent', 'NormalizedCa', 'CaActivation'});
writetable(summary, fullfile(outdir, 'EMG_Ca_activation_response.xlsx'));

disp(['Saved EMG-Ca response figure to ', fullfile(pwd, outdir)]);

function pad = defaultCalciumConfig()
pad = struct();
pad.name = 'Healthy perfusion';
pad.perfusion = 1.0;
pad.ca_rest = 0.05;
pad.ca_max = 1.0;
pad.k_release = 12.0;
pad.k_uptake = 6.0;
pad.k_ca = 0.35;
pad.n_ca = 3.0;
pad.k_pi_ca = 20.0;
pad.k_h_ca = 0.00020;
pad.atp_ref = 8.2;
pad.uptake_floor = 0.30;
pad.pi_accel = 0.50;
pad.h_accel = 0.90;
pad.pcr_loss = 0.20;
pad.atp_loss = 0.03;
pad.recovery_floor = 0.15;
pad.glycolysis_gain = 0.25;
end

function response = simulateCalciumResponse(params, cycles, iemg_profile, data_resting, pad)
TmpC = 37;
MgATP = 8.2;
MgADP = data_resting{1,2}*10^-3;
Pi = data_resting{4,2};
Pcr = data_resting{2,2};
SL0 = 3.23;
pH = data_resting{3,2};
H = 1e3*10^-pH;
N0 = 1;
pad.pi_ref = Pi;
init = [zeros(1,9), N0, SL0, Pi, MgADP, Pcr, H, MgATP, pad.ca_rest];

cycle_time = 10/6.33;
tspan = 0:0.1:cycle_time;
n = length(tspan);
m = length(cycles);

response.Ca = zeros(m,1);
response.Activation = zeros(m,1);

for i = 1:m
    SL_set = 3.23;
    dSL_set = -0.68;
    iemg = iemg_profile(i);
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 5e-3);
    [T, Y] = ode15s(@Model_XB_human_QC_Ca_PAD, tspan, init, options, TmpC, SL_set, params, iemg, dSL_set, Pcr, H, pad);

    init(10) = Y(n,10);
    init(12) = max(Y(n,12), 0);
    init(13) = max(Y(n,13), 0);
    init(14) = max(Y(n,14), 0);
    init(15) = max(Y(n,15), eps);
    init(16) = max(Y(n,16), eps);
    init(17) = min(max(Y(n,17), pad.ca_rest), pad.ca_max);

    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Ca, activation] = ...
        Model_XB_human_QC_Ca_PAD(T(n), Y(n,:), TmpC, SL_set, params, iemg, dSL_set, Pcr, H, pad);
    response.Ca(i) = Ca;
    response.Activation(i) = activation;
end
end

function styleCurrentAxes()
ax = gca;
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.75 0.75 0.75], ...
    'MinorGridColor', [0.85 0.85 0.85], 'Box', 'off');
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
set(get(ax, 'Title'), 'Color', 'k');
set(get(ax, 'XLabel'), 'Color', 'k');
set(get(ax, 'YLabel'), 'Color', 'k');
end
