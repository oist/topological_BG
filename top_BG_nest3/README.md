# Topological Basal Ganglia Model — NEST 3

This directory contains the translation of the topological BG model for [NEST 3](https://www.nest-simulator.org/). The biological parameters, simulation modes, and overall structure are identical to the root (NEST 2.20) version; only the PyNEST API calls have been updated for NEST 3 compatibility.

See the [main README](../README.md) for full model description and citation.

---

## Repository Structure

| File | Description |
|---|---|
| `simParams.py` | Global simulation settings: time step, duration, number of channels, which simulation mode is active (`resting_state`, `action_selection`, or `plasticity`), number of CPUs, and scale factor. **Edit this file to switch between simulation modes.** |
| `bgParams.py` | Biological parameters for all BG nuclei and connections: neuron counts, in-degrees, synaptic weights, delays, asymmetry factors (structural and strength), D1/D2 pathway overlap, and STDP parameters. `IeMSN` is set to 26.0 pA in this version. |
| `fetch_params.py` | Utility to load `simParams.py` and `bgParams.py` at runtime and return them as dictionaries. |
| `nest_routine.py` | Core PyNEST functions (NEST 3 API): kernel initialization, layer creation with 2D/3D topology, focused/diffuse connectivity, Poisson input generators, spike recorders, and firing rate utilities. |
| `ini_all.py` | High-level BG instantiation: creates all nucleus layers, wires them together, and sets up dopamine-modulated STDP synapses for the plasticity mode. |
| `stim_all_model.py` | Main simulation script. Reads parameters, instantiates the network, and runs the selected simulation mode. Saves firing rates and spike data to `./log/`. |
| `go.slurm` | SLURM batch script for running the simulation on an HPC cluster (10 CPUs, 2-day wall time). |

---

## Requirements

- [NEST simulator](https://www.nest-simulator.org/) 3.x with PyNEST
- Python 3.8+
- NumPy

A `log/` directory must exist in the working directory before running:

```bash
mkdir -p log
```

---

## Running the Simulations

All three modes are run with the same command from inside `top_BG_nest3/`:

```bash
python3 stim_all_model.py
```

The active mode is selected by setting `"on": True` for exactly one mode in `simParams.py`. The other two must have `"on": False`.

---

### Resting State

Simulates baseline BG dynamics with no external stimulation. Records mean firing rates for all nuclei.

In `simParams.py`:
```python
"resting_state":    {"on": True,  ...},
"action_selection": {"on": False, ...},
"plasticity":       {"on": False, ...},
```

Output files written to `./log/`: `mean_fr.json`, `At.json`, `performance.txt`, nucleus position files (`{nucleus}.txt`), and one NEST spike file per nucleus.

---

### Action Selection

Drives two cortical channels (CSN/PTN) with a grid of Poisson input rates and records BG output (GPi, MSN, etc.) to evaluate channel competition.

In `simParams.py`:
```python
"resting_state":    {"on": False, ...},
"action_selection": {"on": True,  ...},
"plasticity":       {"on": False, ...},
```

Output files written to `./log/`: `GPi_fr.txt`, `MSN_fr.txt`, `MSN_d1_fr.txt`, `MSN_d2_fr.txt`, `PTN_fr.txt`, `CSN_fr.txt`, `performance.txt`, nucleus position files, and spike files.

---

### Plasticity (Learning)

Runs a conditioning and discrimination protocol with dopamine-modulated STDP on corticostriatal synapses. The simulation proceeds through three phases:

| Phase | Episodes | Description |
|---|---|---|
| CS+ conditioning | 0 – 100 | Channel 1 strong, DA bursts reinforce D1/D2 synapses |
| Generalization probing | 101 – 130 | No stimulus (silence); probe episodes at +5, +10, +15, +20 present CS+ (T1) or CS− (T2) patterns |
| CS− discrimination | 131 – 229 | Channel 2 strong, DA dips; stimulus pattern controlled by `./log/do_tests.txt` (values: `F`, `T1`, `T2`) |

In `simParams.py`:
```python
"resting_state":    {"on": False, ...},
"action_selection": {"on": False, ...},
"plasticity":       {"on": True,  ...},
```

Output files written to `./log/`: `GPi_fr.txt`, `MSN_fr.txt`, `MSN_d1_fr.txt`, `MSN_d2_fr.txt`, `FSI_fr.txt`, `GPe_fr.txt`, `STN_fr.txt`, `CS_plus.txt`, `CS_minus.txt`, `current_task.txt`, `weights_d1.txt`, `weights_d2.txt`, `performance.txt`, nucleus position files, and spike files.

---

### Running on a Cluster (SLURM)

```bash
sbatch go.slurm
```

The job requests 10 CPUs and up to 2 days of wall time. Edit `go.slurm` to adjust resources or paths for your cluster environment.
