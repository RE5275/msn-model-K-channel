import os
import numpy as np
import matplotlib.pyplot as plt
from neuron import h  # Import h from neuron
from scipy import stats  # For statistical tests
import csv  # For writing to CSV files
import pandas as pd

# =========================
# NEURON setup
# =========================
# Load standard run system (needed for h.tstop, h.run, etc.)
h.load_file("stdrun.hoc")

# --- DLL Load Handling ---
dll_path = os.path.abspath("nrnmech.dll")
try:
    h.nrn_load_dll(dll_path)
    print("Mechanisms loaded from:", dll_path)
except RuntimeError:
    print("Mechanisms already loaded or failed to load. Ensure mechanisms were compiled.")

# --- Simulation Parameters ---
h.dt = 0.1      # ms
h.tstop = 700   # total simulation time in ms
h.v_init = -70  # initial membrane potential

print("stdrun.hoc loaded; using DLL at:", dll_path)
print("dt=0.1 ms, tstop=700 ms, v_init=-70 mV")

# =========================
# Cell & stimulus
# =========================
soma = h.Section(name='soma')
soma.L = soma.diam = 20  # um
soma.Ra = 100            # Ohm*cm
soma.cm = 1              # uF/cm2

# Passive props
soma.insert('pas')
soma.g_pas = 0.0001
soma.e_pas = -65

# Custom K+ mechanism
soma.insert('k_mod_pd')
soma.ek = -90  # mV

# Current clamp
stim = h.IClamp(soma(0.5))
stim.delay = 100  # ms
stim.dur = 500    # ms
# stim.amp will be set per trial for variability

# Mechanism instance on the segment
chan = soma(0.5).k_mod_pd

# =========================
# Time windows
# =========================
# Your original windows for local stats/plots (based on stim)
STIM_START = stim.delay
STIM_END   = stim.delay + stim.dur

# Your steady-state window used in your earlier script
STEADY_STATE_START = 500
STEADY_STATE_END   = 650

# Exact windows required for the stats CSV
VM_MEAN_WIN = (100.0, 600.0)   # Vm mean window (ms)
SS_MEAN_WIN = (400.0, 600.0)   # g_mod, cAMP, Ik mean window (ms)

# =========================
# CSV outputs
# =========================
CSV_FILE_NAME   = "simulation_results.csv"      # your original CSV
CSV_STATS_FILE  = "sim_metrics_trials.csv"      # stats-friendly CSV

csv_headers = [
    'Scenario', 'Trial_Number',
    'Vm_Avg_mV', 'IK_Avg_mA_cm2', 'Gmod_Avg_S_cm2', 'cAMP_Avg', 'Stim_Amp_nA'
]

csv_stats_headers = [
    "mechanism","condition","trial",
    "Vm_mean_100_600_mV",
    "cAMP_mean_400_600_uM",
    "gmod_mean_400_600_S_per_cm2",
    "Ik_mean_400_600_mA_per_cm2"
]

def init_csvs():
    sim_out   = open(CSV_FILE_NAME, 'w', newline='')
    stats_out = open(CSV_STATS_FILE, 'w', newline='')
    sim_w     = csv.writer(sim_out)
    stats_w   = csv.writer(stats_out)
    sim_w.writerow(csv_headers)
    stats_w.writerow(csv_stats_headers)
    return sim_out, stats_out, sim_w, stats_w

def close_csvs(sim_out, stats_out):
    sim_out.close()
    stats_out.close()

def write_stats_trial(stats_writer, condition_label, trial_idx_1based,
                      t_arr, v_arr, ik_arr, gmod_arr, camp_arr):
    """Compute windowed means from single-trial vectors and append one row to sim_metrics_trials.csv."""
    vm_mask = (t_arr >= VM_MEAN_WIN[0]) & (t_arr <= VM_MEAN_WIN[1])
    ss_mask = (t_arr >= SS_MEAN_WIN[0]) & (t_arr <= SS_MEAN_WIN[1])

    def safe_mean(arr, mask):
        return float(np.mean(arr[mask])) if np.any(mask) else float("nan")

    Vm_mean_100_600_mV      = safe_mean(v_arr,    vm_mask)
    Ik_mean_400_600_mA_cm2  = safe_mean(ik_arr,   ss_mask)
    gmod_mean_400_600_S     = safe_mean(gmod_arr, ss_mask)
    cAMP_mean_400_600_uM    = safe_mean(camp_arr, ss_mask)

    stats_writer.writerow([
        "k_mod_pd",                       # mechanism
        condition_label,                  # condition (label)
        int(trial_idx_1based),            # trial
        f"{Vm_mean_100_600_mV:.6f}",
        f"{cAMP_mean_400_600_uM:.6f}",
        f"{gmod_mean_400_600_S:.6f}",
        f"{Ik_mean_400_600_mA_cm2:.6f}"
    ])

# =========================
# Simulation helpers
# =========================
def run_condition_and_collect_metrics(label, A2A_on_val, DA_level_val,
                                      num_trials, csv_writer, csv_stats_writer):
    """
    Runs a scenario multiple times, records single-trial vectors using your method,
    writes a row to both CSVs per trial, and returns arrays for local stats.
    """
    print(f"\n--- Running Scenario: {label} ({num_trials} trials) ---")

    vm_avg_trials   = []
    ik_avg_trials   = []
    gmod_avg_trials = []
    camp_avg_trials = []

    for i in range(num_trials):
        # Set per-trial parameters
        h.finitialize(h.v_init)
        h.cvode_active(0)  # fixed time step

        # Slight perturbation for variability (centered at 0.2 nA)
        current_stim_amp = 0.2 + np.random.normal(0, 0.005)
        stim.amp = current_stim_amp

        # Set RANGE params via channel instance
        chan.A2A_on = A2A_on_val
        chan.DA_level = DA_level_val
        chan.gkbar = 0.05
        chan.pd_cAMP_gain = 5.0
        chan.base_cAMP = 1

        # Set GLOBAL params via h
        h.nq_k_mod_pd = 4
        h.prod_rate_k_mod_pd = 1.0
        h.deg_rate_k_mod_pd = 0.1
        h.basal_cAMP_prod_k_mod_pd = 0.1

        # Record *fresh* vectors for this trial to avoid contamination
        t_vec   = h.Vector(); t_vec.record(h._ref_t)
        v_vec   = h.Vector(); v_vec.record(soma(0.5)._ref_v)
        ik_vec  = h.Vector(); ik_vec.record(soma(0.5)._ref_ik_k_mod_pd)  # mechanism-specific IK
        g_vec   = h.Vector(); g_vec.record(chan._ref_g_mod)
        c_vec   = h.Vector(); c_vec.record(chan._ref_cAMP)

        # Run
        h.run()

        # Convert to numpy
        t_arr   = np.array(t_vec)
        v_arr   = np.array(v_vec)
        ik_arr  = np.array(ik_vec)
        g_arr   = np.array(g_vec)
        c_arr   = np.array(c_vec)

        # Local windows (your earlier script's windows)
        stim_mask = (t_arr >= STIM_START) & (t_arr <= STIM_END)
        ss_mask   = (t_arr >= STEADY_STATE_START) & (t_arr <= STEADY_STATE_END)

        def safe_mean(arr, mask):
            return float(np.mean(arr[mask])) if np.any(mask) else float("nan")

        current_vm_avg   = safe_mean(v_arr,  stim_mask)
        current_ik_avg   = safe_mean(ik_arr, stim_mask)
        current_gmod_avg = safe_mean(g_arr,  ss_mask)
        current_camp_avg = safe_mean(c_arr,  ss_mask)

        vm_avg_trials.append(current_vm_avg)
        ik_avg_trials.append(current_ik_avg)
        gmod_avg_trials.append(current_gmod_avg)
        camp_avg_trials.append(current_camp_avg)

        # Write original CSV row
        csv_writer.writerow([
            label, i + 1,
            f"{current_vm_avg:.3f}",
            f"{current_ik_avg:.3f}",
            f"{current_gmod_avg:.5f}",
            f"{current_camp_avg:.3f}",
            f"{current_stim_amp:.4f}"
        ])

        # Write stats-friendly CSV row (exact required columns/windows)
        write_stats_trial(csv_stats_writer, label, i + 1, t_arr, v_arr, ik_arr, g_arr, c_arr)

    return {
        'label': label,
        'vm_avg':   np.array(vm_avg_trials),
        'ik_avg':   np.array(ik_avg_trials),
        'gmod_avg': np.array(gmod_avg_trials),
        'camp_avg': np.array(camp_avg_trials),
    }

def run_single_scenario_and_plot_traces(label, A2A_on_val, DA_level_val, axs):
    """Run one scenario and overlay 4 traces (Vm, IK, g_mod, cAMP) on given axes."""
    print(f"--- Plotting Single Scenario: {label} ---")
    h.finitialize(h.v_init)
    h.cvode_active(0)

    # Set params
    chan.A2A_on = A2A_on_val
    chan.DA_level = DA_level_val
    chan.gkbar = 0.05
    chan.pd_cAMP_gain = 5.0
    chan.base_cAMP = 1

    h.nq_k_mod_pd = 4
    h.prod_rate_k_mod_pd = 1.0
    h.deg_rate_k_mod_pd = 0.1
    h.basal_cAMP_prod_k_mod_pd = 0.1

    stim.amp = 0.2  # base amplitude for clean single-trace view

    # Fresh vectors
    t_vec   = h.Vector(); t_vec.record(h._ref_t)
    v_vec   = h.Vector(); v_vec.record(soma(0.5)._ref_v)
    ik_vec  = h.Vector(); ik_vec.record(soma(0.5)._ref_ik_k_mod_pd)
    g_vec   = h.Vector(); g_vec.record(chan._ref_g_mod)
    c_vec   = h.Vector(); c_vec.record(chan._ref_cAMP)

    h.run()

    axs[0].plot(t_vec, v_vec,   label=label)
    axs[1].plot(t_vec, ik_vec,  label=label)
    axs[2].plot(t_vec, g_vec,   label=label)
    axs[3].plot(t_vec, c_vec,   label=label)

# =========================
# Main execution
# =========================
NUM_TRIALS = 20

# Init CSV files
sim_out, stats_out, csv_writer, csv_stats_writer = init_csvs()

print("\n--- Running Scenarios for Statistical Analysis and CSV Logging ---")
all_results = [
    run_condition_and_collect_metrics("PD + A2A ON",      1, 0, NUM_TRIALS, csv_writer, csv_stats_writer),
    run_condition_and_collect_metrics("PD + A2A OFF",     0, 0, NUM_TRIALS, csv_writer, csv_stats_writer),
    run_condition_and_collect_metrics("Control + A2A ON", 1, 1, NUM_TRIALS, csv_writer, csv_stats_writer),
    run_condition_and_collect_metrics("Control Baseline", 0, 1, NUM_TRIALS, csv_writer, csv_stats_writer),
]

# Close CSVs
close_csvs(sim_out, stats_out)
print(f"\nAll simulation data has been saved to '{CSV_FILE_NAME}'.")
print(f"Stats-ready metrics saved to '{CSV_STATS_FILE}'.")

# =========================
# Quick integrity printout
# =========================
try:
    _df = pd.read_csv(CSV_STATS_FILE)
    print("\nCounts per metric × condition (from stats CSV):")
    print(_df.groupby("condition")[[
        "Vm_mean_100_600_mV",
        "cAMP_mean_400_600_uM",
        "gmod_mean_400_600_S_per_cm2",
        "Ik_mean_400_600_mA_per_cm2"
    ]].apply(lambda x: x.notna().sum()))
except Exception as e:
    print("Could not read stats CSV for integrity printout:", e)

# =========================
# Local statistical analysis + simple bar plots
# =========================
print("\n--- Performing Local Statistical Analysis and Plotting Summary Bar Charts ---")
labels = [res['label'] for res in all_results]
vm_data = [res['vm_avg'] for res in all_results]
ik_data = [res['ik_avg'] for res in all_results]
gmod_data = [res['gmod_avg'] for res in all_results]
camp_data = [res['camp_avg'] for res in all_results]

metrics = {
    "Membrane Potential (Vm_avg)": {"data": vm_data, "ylabel": "Avg Vm (mV)"},
    "Potassium Current (ik_avg)": {"data": ik_data, "ylabel": "Avg IK (mA/cm²)"},
    "Modulated K+ Conductance (g_mod_avg)": {"data": gmod_data, "ylabel": "Avg g_mod (S/cm²)"},
    "cAMP Dynamics (cAMP_avg)": {"data": camp_data, "ylabel": "Avg cAMP"}
}

fig_stats, axs_stats = plt.subplots(len(metrics), 1, figsize=(10, 15), sharex=False)
axs_stats = np.atleast_1d(axs_stats)

for i, (metric_name, metric_info) in enumerate(metrics.items()):
    ax = axs_stats[i]
    data_list = metric_info['data']
    ylabel = metric_info['ylabel']

    means = [np.nanmean(d) for d in data_list]
    sems  = [np.nanstd(d) / np.sqrt(len(d)) if len(d) > 0 else 0 for d in data_list]

    x_pos = np.arange(len(labels))
    ax.bar(x_pos, means, yerr=sems, capsize=5)
    ax.set_ylabel(ylabel)
    ax.set_title(metric_name)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # ANOVA for each metric (ignore NaNs)
    clean_data_list = [d[~np.isnan(d)] for d in data_list]
    if sum(len(d) > 1 for d in clean_data_list) >= 2:
        f_stat, p_value = stats.f_oneway(*clean_data_list)
        ax.text(0.02, 0.95, f'ANOVA p={p_value:.3f}', transform=ax.transAxes, fontsize=9,
                va='top', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5))
        print(f"\n--- ANOVA: {metric_name} ---")
        print(f"F-statistic: {f_stat:.3f}, p-value: {p_value:.3f}")

plt.tight_layout()
plt.show()

print("\n--- Running Selected Scenarios for Individual Trace Plots ---")
fig_traces, axs_traces = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

run_single_scenario_and_plot_traces("PD + A2A ON",      1, 0, axs_traces)
run_single_scenario_and_plot_traces("PD + A2A OFF",     0, 0, axs_traces)
run_single_scenario_and_plot_traces("Control + A2A ON", 1, 1, axs_traces)
run_single_scenario_and_plot_traces("Control Baseline", 0, 1, axs_traces)

axs_traces[0].set_ylabel("Vm (mV)")
axs_traces[0].set_title("Membrane Potential (All Conditions)")
axs_traces[0].legend(loc='upper right', fontsize='small')
axs_traces[0].grid(True)

axs_traces[1].set_ylabel("K+ Current (ik)")
axs_traces[1].set_title("K+ Current")
axs_traces[1].grid(True)

axs_traces[2].set_ylabel("g_mod (S/cm²)")
axs_traces[2].set_title("Modulated K+ Conductance")
axs_traces[2].grid(True)

axs_traces[3].set_ylabel("cAMP")
axs_traces[3].set_title("cAMP Dynamics")
axs_traces[3].set_xlabel("Time (ms)")
axs_traces[3].grid(True)

plt.tight_layout()
plt.show()

print("\nAll simulations and plots complete.")

