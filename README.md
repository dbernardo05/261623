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

```sh
make
./macgregor < parameters.txt
```

The program writes out various parameters and outputs as can be seen in the code. The code for visualizing and analyzing the program output can be found in the Matlab and Python folders of the root directory.

## Bug fixes (November 2026)

Fixes 1–5 apply to all three variants. Fixes 6–8 apply to `scalefree/layer.cpp` (`addSynapse`).

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

**Effect on results.** The executables and outputs in the experiment subfolders were produced before these fixes. Rebuild with `make` before running new simulations, and don't compare new results directly with old ones. In particular, the mean scale-free in-degree drops from about 36 to about 15.

**Known remaining issues** (not yet fixed):

- `drand48()` is called inside OpenMP loops, which is a data race, so runs aren't reproducible.
- The conductance time constants (`Te1`–`Ti2`) in `parameters.txt` have no effect.
- Synapses onto inhibitory neurons ignore the release parameters.
- The connection kernels in `uniform` and `smallworld` use integer division.
- The small-world rewiring loop skips edges.
