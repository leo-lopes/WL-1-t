# 2D Ising Model — Wang-Landau Monte Carlo

A Fortran implementation of the **Wang-Landau algorithm** for the 2D Ising model on a
square lattice with periodic boundary conditions. Instead of running separate Monte
Carlo simulations at many temperatures, Wang-Landau estimates the **density of states**
`g(E)` directly from a single random walk in energy space — the full thermodynamics
(energy, specific heat, entropy) at *any* temperature then follows by reweighting, in
one post-processing step.

## Why Wang-Landau

Standard Metropolis Monte Carlo needs an independent run per temperature and suffers
from critical slowing down near the phase transition. Wang-Landau sidesteps both
problems: the random walk is driven to visit every energy level with equal probability,
which self-corrects the sampling exactly where a fixed-temperature simulation would
get stuck (near `T_c`).

## Algorithm

1. Start with a flat guess for the entropy, `S(E) = 0` for every accessible energy
   level, and a modification factor `f = f0 = 1`.
2. Propose a single-spin flip; accept it with probability
   `min(1, exp(S(E_old) - S(E_new)))` — i.e. the walk is biased *away* from
   energies whose density of states is already overestimated.
3. Whichever state is landed on (old or new), update `S(E) += f` and the visit
   histogram `H(E) += 1`.
4. Periodically check the histogram for **flatness** (every bin within 80% of the
   mean). Once flat, reset `H` and refine the estimate: `f -> f/2`.
5. Repeat until `f` drops below a stopping threshold (`1e-7`) — at that point `S(E)`
   has converged to the true (relative) density of states.

The code also contains a `1/t` variant of the algorithm (Belardinelli & Pereyra,
2007), which replaces the geometric `f -> f/2` schedule with `f = 1/t` once the first
flatness criterion is met, for better asymptotic convergence; toggle it via the
`method` variable in the main program (`"wl"` or `"it"`).

## Requirements

- A Fortran compiler (`gfortran` recommended, any reasonably modern version)

## Build

```bash
gfortran -O2 -o ising_wl ising_wl.f90
```

(adjust the source filename to whatever you save the module + program as)

## Run

The program writes its output to `./data/L<L>/`, and it does **not** create that
directory itself — make it first:

```bash
mkdir -p data/L10
./ising_wl
```

## Configuration

All parameters are set at the top of the `ising_wl` program (no command-line
arguments):

| Parameter | Meaning | Default |
|---|---|---|
| `L` | Linear lattice size (`L x L` spins) | `10` |
| `mc_steps` | Monte Carlo sweeps per inner loop iteration | `1000` |
| `ic` | Initial condition: `"rand"`, `"allp"` (all +1), `"allm"` (all -1) | `"rand"` |
| `method` | `"wl"` (standard Wang-Landau) or `"it"` (1/t variant) | `"wl"` |
| `f0` / `fstop` | Initial / final modification factor | `1.0` / `1e-7` |
| `hist_pr` | Flatness threshold (fraction of the mean histogram count) | `0.8` |
| `dT`, `Tf` | Temperature step and upper bound for the post-run sweep (starts at `T = 0.5`) | `1e-4`, `5.0` |

## Output

| File | Contents |
|---|---|
| `data/L<L>/wl_err_vs_t.dat` | `time`, `err`, `f` — a convergence diagnostic (spread of `S(E)` across energy bins) logged periodically during the run |
| `data/L<L>/wl_S.dat` | `S(E)`, `g(E) = exp(S(E))`, `E` — the converged (normalized) entropy and density of states, one row per accessible energy level |
| `data/L<L>/wl_avgE_Cv.dat` | `<E>`, `Cv`, `T` — average energy and specific heat as a function of temperature, obtained by reweighting `g(E)` over the configured temperature range |

`wl_avgE_Cv.dat` is the payoff: plotting `Cv` against `T` should show a peak growing
and sharpening around the exact 2D Ising critical temperature,
`T_c = 2 / ln(1 + sqrt(2)) ≈ 2.269` (Onsager, 1944) — a good sanity check for the
simulation, and a natural place to extend into a finite-size scaling study by
re-running at several `L`.

## References

- Wang, F. & Landau, D.P. (2001). *Efficient, multiple-range random walk algorithm to
  calculate the density of states.* Physical Review Letters, 86(10), 2050.
- Belardinelli, R.E. & Pereyra, V.D. (2007). *Fast algorithm to calculate density of
  states.* Physical Review E, 75(4), 046701.
- Onsager, L. (1944). *Crystal statistics. I. A two-dimensional model with an
  order-disorder transition.* Physical Review, 65(3-4), 117.
