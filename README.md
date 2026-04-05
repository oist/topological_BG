# Topological Basal Ganglia Model

A biologically detailed spiking neural network model of the basal ganglia (BG) implemented in [PyNEST](https://www.nest-simulator.org/). The model reproduces topological connectivity across BG nuclei and supports three simulation modes: resting state, action selection, and dopamine-dependent learning (plasticity).

---

## Repository Structure

| File | Description |
|---|---|
| `simParams.py` | Global simulation settings: time step, duration, number of channels, which simulation mode is active (`resting_state`, `action_selection`, or `plasticity`), number of CPUs, and scale factor. **Edit this file to switch between simulation modes.** |
| `bgParams.py` | Biological parameters for all BG nuclei and connections: neuron counts, in-degrees, synaptic weights, delays, asymmetry factors (structural and strength), D1/D2 pathway overlap, and STDP parameters. |
| `fetch_params.py` | Utility to load `simParams.py` and `bgParams.py` at runtime and return them as dictionaries. |
| `nest_routine.py` | Core PyNEST functions: kernel initialization, layer creation with 2D/3D topology, focused/diffuse connectivity, Poisson input generators, spike detectors, and firing rate utilities. |
| `ini_all.py` | High-level BG instantiation: creates all nucleus layers, wires them together, and sets up dopamine-modulated STDP synapses for the plasticity mode. |
| `stim_all_model.py` | Main simulation script. Reads parameters, instantiates the network, and runs the selected simulation mode. Saves firing rates and spike data to `./log/`. |
| `task.py` | Standalone module implementing a 2D point-mass motor task and a Q-learning agent (conceptually related to BG reinforcement learning; not used in the main NEST pipeline). |
| `go.slurm` | SLURM batch script for running the simulation on an HPC cluster (10 CPUs, 2-day wall time). |
| `figures/plots.py` | Standalone plotting script that reproduces all paper figures (Figures 1, 3–8) from pre-saved simulation data. Run from inside the `figures/` directory. |
| `figures/run4_*/` | Pre-saved simulation output folders used by `plots.py`. Each folder name encodes the parameter set: κ (kappa), λ (lambda), η (mod), and condition (`_p_no_conv` = Parkinson, `_s` = Schizophrenia). |

---

## Requirements

- [NEST simulator](https://www.nest-simulator.org/) 2.20 with PyNEST
- Python 3.7+
- NumPy

A `log/` directory must exist in the working directory before running:

```bash
mkdir -p log
```

---

## Running the Simulations

All three modes are run with the same command:

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

Output files written to `./log/`: `mean_fr.json`, `At.json`, `performance.txt`, nucleus position files (`{nucleus}.txt`), and one NEST spike file per nucleus (`{nucleus}-0.gdf`, two columns: neuron GID and spike time).

---

### Action Selection

Drives two cortical channels (CSN/PTN) with a grid of Poisson input rates and records BG output (GPi, MSN, etc.) to evaluate channel competition.

In `simParams.py`:
```python
"resting_state":    {"on": False, ...},
"action_selection": {"on": True,  ...},
"plasticity":       {"on": False, ...},
```

Output files written to `./log/`: `GPi_fr.txt`, `MSN_fr.txt`, `MSN_d1_fr.txt`, `MSN_d2_fr.txt`, `PTN_fr.txt`, `CSN_fr.txt`, `performance.txt`, nucleus position files, and spike files (`{nucleus}-0.gdf`).

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

Output files written to `./log/`: `GPi_fr.txt`, `MSN_fr.txt`, `MSN_d1_fr.txt`, `MSN_d2_fr.txt`, `FSI_fr.txt`, `GPe_fr.txt`, `STN_fr.txt`, `CS_plus.txt`, `CS_minus.txt`, `current_task.txt`, `weights_d1.txt`, `weights_d2.txt`, `performance.txt`, nucleus position files, and spike files (`{nucleus}-0.gdf`).

---

### Running on a Cluster (SLURM)

```bash
sbatch go.slurm
```

The job requests 10 CPUs and up to 2 days of wall time. Edit `go.slurm` to adjust resources or paths for your cluster environment.

---

## Citation

If you use this code or data in your work, please cite:

> Gutierrez CE, Lienard J, Girard B, Urakubo H, Doya K (2025). *Topological basal ganglia model with dopamine-modulated spike-timing-dependent plasticity reproduces reinforcement learning, discriminatory learning, and neuropsychiatric disorders.* bioRxiv. https://doi.org/10.1101/2025.11.10.687760

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
