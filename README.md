# b3bound — optimal-pot explorer for B₃ bond-edge types

Menu-driven. Pick a graph, then explore it: colorings, orientations, the lower
bound, failing construction matrices, and (opt-in) isomorphism validation and
recoloring. Every expensive isomorphism step asks first.

Keep `swapcolor.py` (the triangle machinery + graph families) next to
`b3bound.py`.

`B₃(G)` = fewest bond-edge types to self-assemble `G` uniquely (Scenario 3).

## Run

```
python3 b3bound.py                              # interactive menu
python3 b3bound.py --family cycle --range 4 6
python3 b3bound.py --family circular_ladder --range 4
python3 b3bound.py --edges "0-1 1-2 2-0 2-3"
```

`--range` is `start`, `start end`, or `start end step`. `--edges` is `u-v`
tokens. With no flags you get a family menu.

## The menu (per graph)

```
lower bound:  B_3 >= k bond-edge types
  1. show lower bound + summary
  2. visualize colorings                      -> coloring_*.png
  3. visualize valid orientations             -> orient_*.png
  4. dump failing construction matrices       -> failing_construction_matrices.txt
  5. list candidate pots
  6. validate candidates  (isomorphism -- slow, asks first)
  7. try next color up    (exhausts the recolor pile -- slow, asks first)
  0. done
```

All output for a graph goes to a folder `b3_<label>/`.

- **colorings** — every minimum coloring the triangle method admits (there can
  be several; each is an image with edges colored by bond-edge type).
- **orientations** — for each coloring, every loop-rule-valid orientation, drawn
  as a directed graph (arrow = unhatted→hatted).
- **failing construction matrices** — orientations the construction matrix
  rejects (they realize a *smaller* complex, m_P < n). The file lists the
  orientation, the tile types, the Z matrix, and m_P. Useful raw data; there can
  be many, hence a file.
- **candidates** — orientations that pass loop + multiedge + m_P = n. Not yet
  proven optimal.
- **validate** — for each candidate, every allowed same-type arrow swap must
  keep the graph isomorphic to G. Survivors are **optimal tilings**, written as
  `optimal_tiling_*.txt` (pot set form + oriented edges) and `.png`. This is the
  only isomorphism-heavy step; it asks before running.
- **next color up** — if there are no candidates (or none validate), this digs
  through the whole set-aside pile, recolors each failing edge to a fresh color
  (k+1), and re-runs. Also asks first.

## Reading the bound

The number is a **lower bound** (looseness is fine). No candidates at k proves
B₃ ≥ k+1. A failed recolor round doesn't prove more than that — it only tries
one-edge recolorings; the systematic k-loop is future work.

## Requirements

`python3`, `networkx`, `numpy`, `matplotlib`; `swapcolor.py` alongside.