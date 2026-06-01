# How to Run the PAD Calcium Extension

This extension keeps the original repository files intact and adds two new MATLAB files:

- `code/Model_XB_human_QC_Ca_PAD.m`
- `code/figure_PAD_calcium_extension.m`

## Run in MATLAB

1. Open MATLAB.
2. Set the current folder to this repository.
3. Run:

```matlab
cd code
figure_PAD_calcium_extension
```

The script reads the repo's EMG and metabolite data, fits the EMG input the same way as the original figure scripts, runs healthy/moderate PAD/severe PAD simulations, and saves outputs to:

```text
code/figure_PAD_calcium_extension/
```

Expected files:

- `PAD_calcium_extension.png`
- `PAD_calcium_extension.pdf`
- `PAD_calcium_extension_summary.xlsx`

## What the Figure Shows

The figure is organized around the proposed mechanism:

```text
EMG input -> calcium activation -> metabolite accumulation -> force decline
```

The PAD scenarios use lower perfusion values to reduce metabolite recovery and increase accumulation during repeated contractions.

## Where to Tune PAD Severity

Open `code/figure_PAD_calcium_extension.m` and edit:

```matlab
moderate.perfusion = 0.65;
severe.perfusion = 0.45;
```

Lower perfusion means stronger ischemic impairment.

## Where to Tune Calcium

Open `code/figure_PAD_calcium_extension.m` and edit the defaults in:

```matlab
function pad = defaultPadConfig()
```

The most important calcium parameters are:

- `k_release`: how strongly EMG releases calcium
- `k_uptake`: calcium removal rate
- `k_ca`: calcium sensitivity of activation
- `n_ca`: Hill coefficient for calcium activation
- `k_pi_ca`: Pi increase required to substantially inhibit calcium release
- `k_h_ca`: H+ increase required to substantially inhibit calcium release
- `atp_ref`: resting ATP reference for calcium uptake
- `uptake_floor`: minimum calcium uptake retained when ATP falls

## Where PAD Enters the Model

PAD effects are implemented in `code/Model_XB_human_QC_Ca_PAD.m`.

The model:

- converts EMG into calcium-mediated activation
- passes that activation into the original `Model_XB_human_QC`
- reduces metabolite recovery under impaired perfusion
- adds extra Pi and H+ accumulation proportional to cross-bridge demand
- feeds elevated Pi and H+ back into calcium release
- feeds reduced ATP back into calcium uptake

The resulting feedback loop is:

```text
EMG -> calcium activation -> cross-bridge cycling -> metabolite dynamics
          ^                                          |
          |------------- Pi, H+, and ATP ------------|
```
