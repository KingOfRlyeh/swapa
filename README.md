# B₃ bond-edge lower bound

Two programs for the self-assembling-DNA project:

- **`swapcolor.py`** — finds swappability and the minimum *bond-edge coloring* β(G) = χ(H) of a graph (the 2-switch obstruction, with brute-force tie-breaking).
- **`b3bound.py`** — takes that coloring and **tightens it into a lower bound on B₃(G)** using the construction matrix. This is the one with the CLI below. It imports `swapcolor.py`, so keep both files in the same folder.

`B₃(G)` = the minimum number of bond-edge types needed to build `G` as a Scenario-3 self-assembling complex (no smaller complex and no same-size non-isomorphic complex allowed).

---

## What the tool computes

For a graph `G` it prints two numbers:

| column | meaning |
|--------|---------|
| **chi(H)** | Lower bound from the **swap obstruction alone**: two non-adjacent edges must get different bond-edge types when *every* 2-switch on them changes the graph. This is `swapcolor`'s β. |
| **B_3 >=** | The **tightened** lower bound. Starting from chi(H), it also demands the coloring be *realizable*: some orientation of the colored edges must give a pot whose construction matrix realizes `G` **with no smaller complex** (mₚ = n). The smallest number of types for which that is possible. |

Both are **lower bounds**: `chi(H) ≤ B_3 >= ≤ B₃(G)`.

### Reading the `exact` column
- **`yes`** — at every number of types below the reported value, the search was exhaustive and found no realizable coloring, so the bound is **certified**.
- **`no`** — a search cap cut the enumeration short. The printed number is a **floor that may rise**; re-run with a larger `--color-cap` or `--k-max`.

### Strength verdict (printed per graph)
- **TIGHTENED** — the construction matrix pushed the bound above chi(H). The realizability obstruction did real work.
- **BASELINE** — a valid pot already existed at chi(H) types, so no improvement over the pure swap bound.
- **PARTIAL** — the search was truncated (`exact = no`); treat the value as a floor.

### Scope / honesty note
This bound captures the **swap** obstruction and the **no-smaller-complex** obstruction. It does **not** capture the Scenario-3-only "no same-size non-isomorphic complex" obstruction (the loop/multi-edge arguments, Lemmas 2–3 of Ellis-Monaghan et al.). Consequences:
- **Tight** for cycles `Cₙ` and complete bipartite `K_{m,n}`.
- **Loose** for graphs whose B₃ is driven by that missing obstruction — e.g. `Kₙ`, where the true `B₃ = n−1` but this tool returns roughly `B₂ ≈ 1–2`.

---

## Usage

Three ways to run it.

### 1. One-liners with flags
```
python3 b3bound.py --family cycle --range 3 8
python3 b3bound.py --family sun   --range 3 6 --verbose
python3 b3bound.py --edges "0-1 1-2 2-0 2-3"
```

### 2. Interactive (no flags)
```
python3 b3bound.py
```
Gives a numbered menu of families, then prompts for the `n` range and a trace toggle. Pick `custom` to paste an edge list.

### 3. Help
```
python3 b3bound.py --help
```

---

## Flags

| flag | argument | meaning | default |
|------|----------|---------|---------|
| `--family` | one of the names below | which graph family to test | — |
| `--range` | `start` **or** `start end` **or** `start end step` | which values of `n` to run. `--range 3 8` → n = 3,4,5,6,7,8. `--range 3 9 2` → 3,5,7,9. `--range 5` → just n = 5. | family minimum |
| `--edges` | quoted edge list, e.g. `"0-1 1-2 2-0"` | test one custom graph instead of a family. Vertices are integers, edges are `u-v`, space- or comma-separated. | — |
| `--k-max` | integer | highest number of bond-edge types to search up to | `|E(G)|` |
| `--color-cap` | integer | max proper colorings enumerated per type-count. Raise it if you see `exact = no`. | `500000` |
| `--orient-cap` | integer | if `|E| ≤` this, try **all** 2^{\|E\|} edge orientations (exact); above it, sample orientations randomly (may under-count, giving a safe floor). | `14` |
| `--verbose` | (none) | print the per-`k` search trace: for each candidate type-count, whether a realizable pot was found or ruled out. | off |

### The two caps, in plain terms
The tool searches over (colorings × orientations). Both caps bound that search:
- **`--color-cap`** limits how many edge-colorings are tried per type-count. Hitting it flips `exact` to `no`.
- **`--orient-cap`** decides exact-vs-sampled orientations. Graphs with ≤ 14 edges are handled exactly; larger graphs fall back to sampling, so their bound is trustworthy as a floor but a "no realizable pot" conclusion is not certified.

Both only ever make the reported bound **more conservative** (smaller or floored) — they never inflate it.

---

## Graph families (`--family`)

`n` is the family's symmetry index (so `mobius 3` → M₆, `sun 3` → 6 vertices).

| name | graph | vertices | min n |
|------|-------|----------|-------|
| `cycle` | Cₙ | n | 3 |
| `complete` | Kₙ | n | 1 |
| `sunlet` | Cₙ with a pendant per vertex | 2n | 3 |
| `sun` | complete n-sun (clique + outer vertices) | 2n | 3 |
| `mobius` | Möbius ladder M₂ₙ | 2n | 2 |
| `ladder` | ladder Lₙ | 2n | 2 |
| `circular_ladder` | prism CLₙ | 2n | 3 |
| `star` | star Sₙ | n+1 | 1 |
| `web` | CLₙ with a pendant per outer vertex | 3n | 3 |
| `custom` | your own edge list (via `--edges` or the interactive menu) | — | — |

---

## Example output

```
$ python3 b3bound.py --family cycle --range 4 6

graph          verts edges  chi(H)  B_3 >=  exact
-------------------------------------------------
cycle(4)           4     4       1       2    yes
cycle(5)           5     5       1       3    yes
cycle(6)           6     6       1       3    yes
```

Here chi(H) = 1 (cycle swaps are mostly harmless on their own) but the
construction matrix raises each bound — matching the known `B₃(Cₙ) = ⌈n/2⌉`
(2, 3, 3 for n = 4, 5, 6). Every row is `TIGHTENED` and `exact = yes`.

---

## Performance

The search is exponential in `|E(G)|` (colorings × orientations × exact rational
construction matrix). It is comfortable up to roughly **6–7 edges**; larger
graphs get slow or fall back to sampled orientations. For bigger inputs, lower
`--k-max` (if you only need to confirm a bound up to some value) or accept the
sampled-orientation floor.

## Requirements
`python3`, `networkx`, `numpy`, `sympy`, and `swapcolor.py` in the same directory.