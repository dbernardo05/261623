# MacGregor neuron network simulations

The subfolder names `uniform`, `smallworld`, and `scalefree` refer to the network connectivity for the simulations contained in those folders. The code in the `uniform` folder has additional comments on the code, when compared to the other two. Each folder has their own READMEs as well.

The name of the executable is `macgregor`. In each of the folders, you will also find the executables, scripts, and outputs for independent simulations with different values for the parameter of interest (neuron neighborhood for uniform, percentage of rewired connections in small world, and slope of line from scale free definitions for the scale free configuration).

Each of the folders contains:

1. `Makefile`
2. `parameters.txt` (used by the program to read in the various parameters)
3. `main.cpp`
4. `neuron.cpp` and `neuron.h`, containing all of the code related to individual neurons
5. `layer.cpp` and `layer.h`, containing all of the code related to the network layer

OpenMP is used to parallelize the for loops.

## Building and running

Each topology has two build targets. `make` builds the corrected `macgregor`. `make reference` builds `macgregor-reference` from the locally stored ModelDB commit `8d5c835`, with the parameter-pointer startup fix and the recording backport. The reference target uses the original neuron update, random draws and wiring algorithm. Use GCC with OpenMP. On this Mac, build with `COMPILER=g++-16` so the reference uses GNU's C++ standard library.

From the repository root, this example builds both small-world executables and runs the reference in its own output directory. Substitute `uniform` or `scalefree` for the other topologies. Each run writes into its current directory.

```sh
make -C smallworld COMPILER=g++-16
make -C smallworld reference COMPILER=g++-16
mkdir -p runs/reference-01
(cd runs/reference-01 && ../../smallworld/macgregor-reference < ../../smallworld/parameters.txt > electrode.txt)
```

### Default recordings

Both builds use a 100 × 100 grid and run 60,000 steps by default. `--steps N` sets the exclusive end step. Steps and neuron IDs start at zero. The neuron at `(row, column)` has ID `row * columns + column`. A cell is inhibitory when both coordinates are multiples of five. The TSV files have column headers.

| File | Contents and timing |
|---|---|
| `field_potentials.tsv` | One row per step with `step`, `legacy` and `fixed`. Both electrodes sample the same membrane voltages **after** the update. Use `legacy` for scoring on the original electrode. |
| `spikes.txt` | The existing headerless `step excitatory_count inhibitory_count` rows, sampled **before** the update. |
| `spike_events.tsv` | One `step, neuron` row for each neuron with `S == 1` at the same sample as `spikes.txt`. A cell that stays active for consecutive steps appears in each bin. Grouping events by step and cell type reproduces the population counts. |
| `edges.tsv` | One row per final directed synapse, after wiring and rewiring, with `pre_neuron`, `post_neuron`, `weight` and `excitatory` (1 or 0 for the presynaptic cell). Parallel synapses have separate rows. |

For squared distance `r²` from the central grid cell and `σ = 2`, the legacy electrode uses `exp(−r² σ² / 2)` and the fixed electrode uses `exp(−r² / (2σ²))`. Both retain the original factor `1 / sqrt(2πσ²)` and sum in grid order. Standard output keeps its two columns, potential and step, with the fixed electrode for `macgregor` and the legacy electrode for `macgregor-reference`. Fresh scale-free runs also print the total connection count at startup.

### Optional output

- `--per-neuron-output` enables the three dense startup maps (`glutamatemap.txt`, `initialglutamatemap.txt`, `weightmap.txt`) and seven traces (`potential.txt`, `releasedglutamate.txt`, `releasableglutamate.txt`, `availableglutamate.txt`, `gaba.txt`, `spikepatterns.txt`, `kconductance.txt`). The traces sample before the update and store one value per neuron in grid order. Small-world records every step, uniform records steps 5,001–64,999 and scale-free records steps 5,001–9,999, subject to the run end. Scale-free also prints its incoming counts. A full small-world run with this flag can write roughly 35 GB. Both builds leave dense files untouched by default.
- `--epoch-steps N` writes `epochs.tsv` with excitatory and inhibitory spike totals and the mean of each electrode. `neuron_epochs.tsv` gives each cell's spike count and mean post-update `glutamate`, the model's count of received vesicles in its PSP history. Epochs use `[first_step, end_step)` and include the final partial epoch. Divide a cell's count by the epoch length for its mean activity per step.
- `--snapshot-steps N` writes `snapshot-N.bin`, `snapshot-2N.bin` and subsequent multiples after those numbers of completed steps. Each file contains a full restart state. Snapshot size depends on edge count and occupied synaptic history slots. Both epochs and snapshots are off until a cadence is supplied.

### Restarting a snapshot

From a run directory, the following example checkpoints the corrected small-world model and continues the same trajectory. `--steps` is the absolute end step on both invocations. The snapshot supplies the parameters on restart.

```sh
OMP_NUM_THREADS=1 ../../smallworld/macgregor --steps 1200 --epoch-steps 500 --snapshot-steps 1200 < ../../smallworld/parameters.txt
OMP_NUM_THREADS=1 ../../smallworld/macgregor --resume snapshot-1200.bin --steps 2000 --epoch-steps 500
```

Use the matching topology, model build and standard library on a POSIX system. The resumed files have a `from-1200-` prefix. The first spike sample is step 1200. Startup maps describe the initial construction and are written only on a fresh run with dense output enabled. Keeping the epoch cadence carries a partial epoch through the restart. When combining segments, use the continued epoch row in place of a partial row with the same start. Changing the cadence begins a new epoch at the restart step.

Snapshots store the input parameters, neuron voltages and conductances, firing and refractory state, ordered synapses, vesicle pools, delayed PSP histories, the random state driving the model and any partial epoch. Restarts restore the saved OpenMP team size. The corrected build also saves the normal generator and its cached variate for each OpenMP worker. Use one OpenMP thread for bit-for-bit replay. The existing shared `drand48()` calls make parallel trajectories scheduling-dependent.

For reading voltage images directly, `MGSNAP01` starts with an 8-byte magic string and six native-endian 32-bit integers (byte-order marker `0x01020304`, topology 1/2/3 for uniform/small-world/scale-free, reference flag, rows, columns, next step). The 18 parameter doubles and 32-bit `NRmax` follow in input order. The next `rows * columns` records are in grid order, each with seven binary64 values (`E`, `Th`, `Gk`, `Ge`, `Gi`, `glutamate`, `gaba`) and two 32-bit integers (`S`, `P`). Fields are packed without padding.

The older standalone sources `scalefree/mainRN.cpp` and `scalefree/main_scalefree.cpp` also produce the compact spike, edge and paired field recordings when rebuilt. They keep their original electrode widths and run lengths. `--per-neuron-output` enables the potential trace in `mainRN.cpp` and the potential, glutamate and spike-pattern traces in `main_scalefree.cpp`. The bundled experiment scripts invoke archived binaries that predate these changes.

Run `CXX=g++-16 python3 tests/test_recording.py` from the repository root to check both builds against a bounded reference simulation and verify checkpoint replay.

## Bug fixes (November 2026)

At commit `81265ef`, fixes 1–5 apply to all three variants. Fixes 6–8 apply to `scalefree/layer.cpp` (`addSynapse`).

| # | Bug | Fix |
|---|-----|-----|
| 1 | **Crash on startup.** `Param* param` was never initialized, so `parse()` wrote through a wild pointer (`main.cpp`). | Declare `Param param` and pass `&param`. |
| 2 | **Membrane noise was a constant.** A new, identically seeded random generator was created on every `computePotential` call, so every neuron got the same value at every step (`neuron.cpp`). | One `mt19937` generator per OpenMP thread (`static thread_local`), seeded once. |
| 3 | **Firing was a hard threshold.** Integer division (`i/100`) made 99 of the 100 `firing_Prob` entries zero (`neuron.cpp`). | Use `i/100.0`, giving the intended soft threshold. |
| 4 | **Recorded potential was dominated by one neuron.** Operator precedence in `/(double)2*sigma*sigma` gave the centre neuron 62% of the electrode weight instead of 4% (`layer.cpp`, `recordPotential`). | Divide by `(2*sigma*sigma)`. |
| 5 | **Vesicle release corrupted the PSP buffer.** Releasing *k* vesicles spread them over *k* past time steps and shifted the synapse's whole history *k* steps back (`neuron.cpp`, `glutamaterelease`). | Add all vesicles to the current step: `psps[0] += Glu_Rel*weight`. The vesicle counter now counts `psps[j]/weight`, so `releasedglutamate.txt` still reports vesicles. |
| 6 | **Scale-free inputs came only from the rows above.** Candidates were accepted in row-by-row order until the degree cap was reached, so the network was wired top to bottom. | Collect all candidates that pass the distance test, then keep a random `d` of them. |
| 7 | **Wrong degree distribution.** `d = (k+1)*minconn` made every degree a multiple of `minconn` (up to 777, above `maxconn`), and the sampling followed exponent λ+1 instead of λ. | Inverse-CDF sampling of P(d) ∝ d^−λ for `minconn` ≤ d ≤ `maxconn`. |
| 8 | **Out-of-bounds read** one element past the end of `scalefreeprob`. | Removed by the new sampler. |

**Effect on results.** The executables and outputs in the experiment subfolders were produced before these fixes. Use `macgregor-reference` for comparisons to the deposited model and `macgregor` for the corrected dynamics. In particular, the mean scale-free in-degree drops from about 36 to about 15.

**Known remaining issues** (not yet fixed):

- `drand48()` is called inside OpenMP loops, which is a data race, so runs aren't reproducible.
- The conductance time constants (`Te1`–`Ti2`) in `parameters.txt` have no effect.
- Synapses onto inhibitory neurons ignore the release parameters.
- The connection kernels in `uniform` and `smallworld` use integer division.
- The small-world rewiring loop skips edges.
