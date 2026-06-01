% Simulate the original EMG-driven cross-bridge model with an added calcium
% activation state and PAD-like ischemic metabolite effects.
%
% Run from the code folder:
%   cd code
%   figure_PAD_calcium_extension

clear variables
clear; close all;

outdir = 'figure_PAD_calcium_extension';
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

healthy = defaultPadConfig();
healthy.name = 'Healthy perfusion';
healthy.perfusion = 1.00;

moderate = defaultPadConfig();
moderate.name = 'Moderate PAD';
moderate.perfusion = 0.65;

severe = defaultPadConfig();
severe.name = 'Severe PAD';
severe.perfusion = 0.45;

scenarios = [healthy moderate severe];
colors = [0.05 0.05 0.05; 0.10 0.35 0.75; 0.75 0.15 0.10];

results = cell(numel(scenarios), 1);
for s = 1:numel(scenarios)
    results{s} = simulatePadScenario(params, cycles, iemg_profile, data_resting, scenarios(s));
end

fig = figure(1);
clf;
set(fig, 'Color', 'w', 'Units', 'Inches', 'Position', [1 1 9.5 6.5]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(x_emg, y_emg, '.', 'Color', [0.55 0.55 0.55], 'MarkerSize', 9); hold on;
plot(cycles, iemg_profile*100, 'k', 'LineWidth', 2);
xlabel('Cycle Index');
ylabel('Normalized EMG (%)');
title('EMG input');
xlim([0 max(cycles)]);
styleCurrentAxes();

nexttile;
for s = 1:numel(scenarios)
    plot(cycles, results{s}.Ca, 'LineWidth', 2, 'Color', colors(s,:)); hold on;
end
xlabel('Cycle Index');
ylabel('Normalized Ca');
title('Calcium activation state');
xlim([0 max(cycles)]);
styleCurrentAxes();

nexttile;
for s = 1:numel(scenarios)
    plot(cycles, results{s}.Pi, 'LineWidth', 2, 'Color', colors(s,:)); hold on;
end
xlabel('Cycle Index');
ylabel('Pi (mM)');
title('Inorganic phosphate');
xlim([0 max(cycles)]);
styleCurrentAxes();

nexttile;
for s = 1:numel(scenarios)
    plot(cycles, results{s}.pH, 'LineWidth', 2, 'Color', colors(s,:)); hold on;
end
xlabel('Cycle Index');
ylabel('pH');
title('Acidosis');
xlim([0 max(cycles)]);
styleCurrentAxes();

nexttile;
for s = 1:numel(scenarios)
    plot(cycles, results{s}.ATP, 'LineWidth', 2, 'Color', colors(s,:)); hold on;
end
xlabel('Cycle Index');
ylabel('MgATP (mM)');
title('ATP availability');
xlim([0 max(cycles)]);
styleCurrentAxes();

nexttile;
for s = 1:numel(scenarios)
    plot(cycles, results{s}.Force, 'LineWidth', 2, 'Color', colors(s,:)); hold on;
end
xlabel('Cycle Index');
ylabel('Force (N)');
title('Force generation');
xlim([0 max(cycles)]);
lgd = legend({scenarios.name}, 'Location', 'best', 'Box', 'off');
set(lgd, 'Color', 'w', 'TextColor', 'k');
styleCurrentAxes();

set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
set(findall(fig, 'Type', 'axes'), 'Color', 'w');
set(findall(fig, '-property', 'XColor'), 'XColor', 'k');
set(findall(fig, '-property', 'YColor'), 'YColor', 'k');
set(findall(fig, '-property', 'TextColor'), 'TextColor', 'k');
exportgraphics(fig, fullfile(outdir, 'PAD_calcium_extension.png'), 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, fullfile(outdir, 'PAD_calcium_extension.pdf'), 'BackgroundColor', 'w', 'ContentType', 'vector');

summary = table();
for s = 1:numel(scenarios)
    r = results{s};
    scenarioSummary = table( ...
        string(scenarios(s).name), ...
        scenarios(s).perfusion, ...
        r.Force(1), ...
        r.Force(end), ...
        100*(r.Force(end)/r.Force(1)), ...
        r.Pi(end), ...
        r.pH(end), ...
        r.ATP(end), ...
        r.Ca(end), ...
        r.Activation(end), ...
        'VariableNames', {'Scenario', 'Perfusion', 'InitialForce_N', 'FinalForce_N', 'FinalForce_percent_initial', 'FinalPi_mM', 'FinalpH', 'FinalMgATP_mM', 'FinalNormalizedCa', 'FinalCaActivation'});
    summary = [summary; scenarioSummary]; %#ok<AGROW>
end
writetable(summary, fullfile(outdir, 'PAD_calcium_extension_summary.xlsx'));
disp(summary);
disp(['Saved results to ', fullfile(pwd, outdir)]);

function pad = defaultPadConfig()
pad = struct();
pad.name = 'Scenario';
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

function result = simulatePadScenario(params, cycles, iemg_profile, data_resting, pad)
TmpC = 37;
MgATP = 8.2;
MgADP = data_resting{1,2}*10^-3;
Pi = data_resting{4,2};
Pcr = data_resting{2,2};
SL0 = 3.23;
pH = data_resting{3,2};
H = 1e3*10^-pH;
N0 = 1;
Ca0 = pad.ca_rest;
pad.pi_ref = Pi;
init = [zeros(1,9), N0, SL0, Pi, MgADP, Pcr, H, MgATP, Ca0];

cycle_time = 10/6.33;
tspan = 0:0.1:cycle_time;
n = length(tspan);
m = length(cycles);

result.Pi = zeros(m,1);
result.ADP = zeros(m,1);
result.PCr = zeros(m,1);
result.H = zeros(m,1);
result.pH = zeros(m,1);
result.ATP = zeros(m,1);
result.Ca = zeros(m,1);
result.Activation = zeros(m,1);
result.Force = zeros(m,1);

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

    result.Pi(i) = init(12);
    result.ADP(i) = init(13);
    result.PCr(i) = init(14);
    result.H(i) = init(15);
    result.pH(i) = -log10(init(15)*10^-3);
    result.ATP(i) = init(16);
    result.Ca(i) = init(17);

    [~, Ftotal, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, activation] = ...
        Model_XB_human_QC_Ca_PAD(T(n), Y(n,:), TmpC, SL_set, params, iemg, dSL_set, Pcr, H, pad);
    result.Force(i) = Ftotal;
    result.Activation(i) = activation;
end
end
