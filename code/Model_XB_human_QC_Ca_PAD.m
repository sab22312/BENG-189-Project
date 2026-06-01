function [dYdT, Ftotal, F_active, F_passive, F_active_s, sov_thick, B_process, C_process, dCK, dK3, dCK_r, dPi_cons, dH_cons, dGly, dAdk, Ca, activation, perfusion] = Model_XB_human_QC_Ca_PAD(t, y, TmpC, SLset, par, iemg, dSL_set, Pcr0, H0, pad)
% Model_XB_human_QC_Ca_PAD extends the original cross-bridge model with a
% calcium activation state and a perfusion-dependent ischemia layer.
%
% State 1:16 are identical to Model_XB_human_QC.
% State 17 is normalized cytosolic calcium.
%
% The original model treated EMG as direct thin-filament activation. Here,
% EMG drives calcium release, calcium determines activation, and activation
% is passed into the original cross-bridge model. Existing metabolite states
% feed back into calcium handling: elevated Pi and H+ inhibit release, while
% reduced MgATP limits uptake.

if nargin < 10 || isempty(pad)
    pad = struct();
end

perfusion = getFieldWithDefault(pad, 'perfusion', 1.0);
perfusion = min(max(perfusion, 0), 1);
ischemia = 1 - perfusion;

ca_rest = getFieldWithDefault(pad, 'ca_rest', 0.05);
ca_max = getFieldWithDefault(pad, 'ca_max', 1.0);
k_release = getFieldWithDefault(pad, 'k_release', 12.0);
k_uptake = getFieldWithDefault(pad, 'k_uptake', 6.0);
k_ca = getFieldWithDefault(pad, 'k_ca', 0.35);
n_ca = getFieldWithDefault(pad, 'n_ca', 3.0);
k_pi_ca = getFieldWithDefault(pad, 'k_pi_ca', 20.0);
k_h_ca = getFieldWithDefault(pad, 'k_h_ca', 0.00020);
atp_ref = getFieldWithDefault(pad, 'atp_ref', 8.2);
uptake_floor = getFieldWithDefault(pad, 'uptake_floor', 0.30);

pi_accel = getFieldWithDefault(pad, 'pi_accel', 0.50);
h_accel = getFieldWithDefault(pad, 'h_accel', 0.90);
pcr_loss = getFieldWithDefault(pad, 'pcr_loss', 0.20);
atp_loss = getFieldWithDefault(pad, 'atp_loss', 0.03);
recovery_floor = getFieldWithDefault(pad, 'recovery_floor', 0.15);
glycolysis_gain = getFieldWithDefault(pad, 'glycolysis_gain', 0.25);

Ca = y(17);
Pi = max(y(12), 0);
H = max(y(15), 0);
MgATP = max(y(16), 0);

pi_ca_factor = 1/(1 + max(Pi - getFieldWithDefault(pad, 'pi_ref', Pi), 0)/k_pi_ca);
h_ca_factor = 1/(1 + max(H - H0, 0)/k_h_ca);
atp_ca_factor = min(max(MgATP/atp_ref, uptake_floor), 1);

effective_release = k_release*pi_ca_factor*h_ca_factor;
effective_uptake = k_uptake*atp_ca_factor;
dCa = effective_release*max(iemg, 0)*(ca_max - Ca) - effective_uptake*(Ca - ca_rest);
activation = Ca^n_ca/(Ca^n_ca + k_ca^n_ca);
activation = min(max(activation, 0), 1);

par_pad = par;
recovery_scale = recovery_floor + (1 - recovery_floor)*perfusion;
par_pad(15) = par_pad(15)*recovery_scale;              % Pi consumption/clearance
par_pad(17) = par_pad(17)*recovery_scale;              % PCr recovery through CK reverse
par_pad(18) = par_pad(18)*(1 + glycolysis_gain*ischemia); % compensatory glycolytic drive

[dBase, Ftotal, F_active, F_passive, F_active_s, sov_thick, B_process, C_process, dCK, dK3, dCK_r, dPi_cons, dH_cons, dGly, dAdk] = ...
    Model_XB_human_QC(t, y(1:16), TmpC, SLset, par_pad, activation, dSL_set, Pcr0, H0);

% PAD is represented as impaired metabolite recovery plus extra accumulation
% proportional to ATPase/cross-bridge cycling demand.
xb_demand = max(dK3, 0);
dBase(12) = dBase(12) + ischemia*pi_accel*xb_demand; % Pi accumulation
dBase(15) = dBase(15) + ischemia*h_accel*xb_demand;  % H+ accumulation
dBase(14) = dBase(14) - ischemia*pcr_loss*max(dCK, 0);
dBase(16) = dBase(16) - ischemia*atp_loss*xb_demand;

dYdT = [dBase; dCa];
end

function value = getFieldWithDefault(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
