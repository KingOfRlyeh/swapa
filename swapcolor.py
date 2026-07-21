import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
from networkx.algorithms.isomorphism import GraphMatcher

# ============================================================
#  Families (n = symmetry index)
# ============================================================

def sunlet_graph(n):
    G = nx.cycle_graph(n)
    for i in range(n): G.add_edge(i, n + i)
    return G

def sun_graph(n):
    G = nx.complete_graph(n)
    for i in range(n): G.add_edge(n+i, i); G.add_edge(n+i, (i+1) % n)
    return G

def mobius_ladder_graph(n):
    G = nx.cycle_graph(2*n)
    for i in range(n): G.add_edge(i, i+n)
    return G

def web_graph(n):
    G = nx.Graph()
    for i in range(n):
        G.add_edge(i, (i+1) % n); G.add_edge(n+i, n+(i+1) % n)
        G.add_edge(i, n+i);       G.add_edge(n+i, 2*n+i)
    return G

def kayak_paddle_graph(c1, c2=None, l=1):
    """Kayak paddle, tuples read as cycle-cycle-path: a C_c1 blade, a C_c2
    blade, joined by a shaft path of length l. So kayak(3,4,2) is C_3 and C_4
    joined by a 2-edge path. Nodes: blade A = 0..c1-1, blade B = c1..c1+c2-1,
    and l-1 internal shaft vertices c1+c2..c1+c2+l-2 (c1+c2+l-1 verts,
    c1+c2+l edges). Defaults c2=c1, l=1 (single bridge edge)."""
    if c2 is None:
        c2 = c1
    G = nx.cycle_graph(c1)                                       # blade A: 0..c1-1
    G.add_edges_from((c1+i, c1+(i+1) % c2) for i in range(c2))   # blade B: c1..c1+c2-1
    prev = 0                                                     # shaft: attach at A-vertex 0
    for j in range(l - 1):                                       # ... through internal verts
        v = c1 + c2 + j
        G.add_edge(prev, v); prev = v
    G.add_edge(prev, c1)                                         # ... to B-vertex c1
    return G

def kayak_params(t):
    """Normalize a kayak tuple to (c1, c2, l) in cycle-cycle-path order.
    (n,) -> (n,n,1) [bridge]; (a,b) -> (a,b,1) [two cycles, bridge]."""
    t = tuple(t)
    if len(t) == 1: return (t[0], t[0], 1)
    if len(t) == 2: return (t[0], t[1], 1)
    return (t[0], t[1], t[2])

# Families whose members are picked by tuple parameters, not a single index n.
TUPLE_FAMILIES = {"kayak"}

FAMILY_GENS = {
    "cycle":           (nx.cycle_graph,           "C_n  (n verts)",      3),
    "complete":        (nx.complete_graph,        "K_n  (n verts)",      1),
    "sunlet":          (sunlet_graph,             "C_n + pendants (2n)", 3),
    "sun":             (sun_graph,                "complete n-sun (2n)", 3),
    "mobius":          (mobius_ladder_graph,      "M_2n (2n verts)",     2),
    "ladder":          (nx.ladder_graph,          "L_n  (2n verts)",     2),
    "circular_ladder": (nx.circular_ladder_graph, "CL_n (2n verts)",     3),
    "star":            (nx.star_graph,            "S_n  (n+1 verts)",    1),
    "web":             (web_graph,                "CL_n + pendants (3n)",3),
    "kayak":           (kayak_paddle_graph,       "cycle-cycle-path, tuples (c1,c2,l)", 3),
}

# ============================================================
#  Delta mechanism (triangle-count; no isomorphism tests here)
# ============================================================

def ekey(u, v):
    return (u, v) if u <= v else (v, u)

def alternating_4cycles(G):
    """Yield (v,w,s,t, delta). d_ij=(A^2)[i,j];
    delta = (d_vt+d_sw)-(d_vw+d_st)-2(A[v,s]+A[w,t])."""
    nodes = list(G.nodes()); idx = {v: k for k, v in enumerate(nodes)}
    A = nx.to_numpy_array(G, nodelist=nodes, dtype=np.int64)
    D = A @ A
    edges = [(idx[u], idx[v]) for u, v in G.edges()]
    for (v, w), (a, b) in combinations(edges, 2):
        for s, t in ((a, b), (b, a)):
            if len({v, w, s, t}) < 4: continue
            if A[v, t] or A[s, w]:    continue
            delta = int((D[v, t] + D[s, w]) - (D[v, w] + D[s, t])
                        - 2 * (A[v, s] + A[w, t]))
            yield nodes[v], nodes[w], nodes[s], nodes[t], delta

def switch(G, v, w, s, t):
    """2-switch: drop vw,st ; add vt,sw. Used ONLY to resolve Delta==0 ties."""
    H = G.copy()
    H.remove_edges_from([(v, w), (s, t)]); H.add_edges_from([(v, t), (s, w)])
    return H

def analyze(G):
    """Cheap pass. Returns (targets, H, ambiguous).
    H: conflict graph on E(G); edge iff EVERY switch of that pair is certified
       (or later confirmed) non-isomorphic -- a pair has up to two distinct
       2-switches, and the edges are forced apart only if BOTH are non-iso.
       One non-iso witness alone is NOT sufficient: if the other switch is
       isomorphic, that pair may still share a bond type via that reconnection.
    ambiguous: {frozenset({e1,e2}): [ (v,w,s,t), ... ]} listing EVERY Delta==0
       reconnection of a pair that has at least one Delta==0 switch (a pair with
       all switches Delta!=0 is decided immediately; nothing left to check)."""
    H = nx.Graph(); H.add_nodes_from(ekey(*e) for e in G.edges())
    pair_switches = {}
    targets = 0
    for v, w, s, t, delta in alternating_4cycles(G):
        targets += 1
        e1, e2 = ekey(v, w), ekey(s, t)
        pair_switches.setdefault(frozenset((e1, e2)), []).append((v, w, s, t, delta))
    ambiguous = {}
    for pair, switches in pair_switches.items():
        zero = [(v, w, s, t) for (v, w, s, t, delta) in switches if delta == 0]
        if zero:
            ambiguous[pair] = zero        # need iso check before deciding this pair
        else:
            H.add_edge(*tuple(pair))      # every switch already Delta != 0: forced
    return targets, H, ambiguous

# ============================================================
#  Coloring of the conflict graph (bond-edge assignment)
# ============================================================

def enumerate_colorings(H, k, cap=2000):
    """Partitions of V(H) into <=k independent sets, symmetry-broken over
    color permutations (each partition once)."""
    nodes = sorted(H.nodes(), key=lambda x: (-H.degree(x), x))
    adj = {u: set(H[u]) for u in nodes}
    out, col = [], {}
    def bt(i, used):
        if len(out) >= cap: return
        if i == len(nodes):
            classes = {}
            for e, c in col.items(): classes.setdefault(c, set()).add(e)
            out.append(frozenset(frozenset(s) for s in classes.values()))
            return
        u = nodes[i]; forb = {col[x] for x in adj[u] if x in col}
        for c in range(min(used + 1, k)):
            if c in forb: continue
            col[u] = c; bt(i + 1, max(used, c + 1)); del col[u]
    bt(0, 0)
    return out

def greedy_chi(H):
    """Exact chi(H) via existence backtracking (cap=1), no enumeration."""
    if H.number_of_nodes() == 0: return 1
    for k in range(1, H.number_of_nodes() + 1):
        if enumerate_colorings(H, k, cap=1): return k
    return H.number_of_nodes()

def min_coloring(H):
    """(chi, one minimum partition as list of frozensets)."""
    if H.number_of_nodes() == 0: return 1, []
    k = greedy_chi(H)
    cs = enumerate_colorings(H, k, cap=1)
    part = list(cs[0]) if cs else [frozenset([v]) for v in H.nodes()]
    return k, part

def chromatic_and_colorings(H, cap=2000):
    if H.number_of_nodes() == 0: return 1, []
    k = greedy_chi(H)
    return k, [p for p in enumerate_colorings(H, k, cap) if len(p) == k] or \
              enumerate_colorings(H, k, cap)

def _color_map(part):
    return {e: i for i, cl in enumerate(part) for e in cl}

# ============================================================
#  Lazy refinement: brute-check only same-colored ambiguous pairs
# ============================================================

def lazy_refine(G, H, ambiguous):
    """Color H; for each Delta==0 pair placed in the SAME color class, run the
    brute-force switch test on EVERY reconnection of that pair (not just the
    first) -- forced conflict requires ALL of a pair's switches to be non-iso;
    if even one is isomorphic, the pair is a genuine free swap via that
    reconnection. Repeat to fixpoint. At fixpoint chi(H) == chi(H*) (true
    bond-edge number). Returns stats."""
    resolved = set(); checks = 0; added = 0
    while True:
        k, part = min_coloring(H); color = _color_map(part)
        changed = False
        for pair, switches in ambiguous.items():
            if pair in resolved: continue
            e1, e2 = tuple(pair)
            if color[e1] == color[e2]:
                iso_flags = []
                for (v, w, s, t) in switches:          # test every reconnection
                    checks += 1
                    iso_flags.append(nx.is_isomorphic(G, switch(G, v, w, s, t)))
                if not any(iso_flags):                 # all non-iso => must differ
                    H.add_edge(e1, e2); added += 1
                    changed = True
                else:
                    resolved.add(pair)                 # at least one iso: free
        if not changed:
            return dict(chi=k, coloring=part, checks=checks,
                        free_swaps=len(resolved), new_conflicts=added)

def complete_conflicts(G, H, ambiguous):
    """Resolve EVERY ambiguous pair by brute force -> true conflict graph H*.
    A pair is forced apart only if ALL of its switches are non-isomorphic;
    if any switch is isomorphic, the pair is free via that specific
    reconnection. Needed to enumerate ALL minimum colorings correctly.
    Returns (H*, stats, safe_switches): safe_switches maps each FREE pair to
    the list of its isomorphic (v,w,s,t) reconnections -- the only ones an
    orientation may realize if the pair ends up same-colored."""
    H2 = H.copy(); free = added = checks = 0
    safe_switches = {}
    for pair, switches in ambiguous.items():
        results = []
        for (v, w, s, t) in switches:
            checks += 1
            results.append(((v, w, s, t), nx.is_isomorphic(G, switch(G, v, w, s, t))))
        if not any(iso for _, iso in results):
            H2.add_edge(*tuple(pair)); added += 1
        else:
            free += 1
            safe_switches[pair] = [sw for sw, iso in results if iso]
    return H2, dict(checks=checks, free_swaps=free, new_conflicts=added), safe_switches

# ============================================================
#  Automorphism dedupe
# ============================================================

def automorphisms(G, cap=20000):
    auts = []
    for m in GraphMatcher(G, G).isomorphisms_iter():
        auts.append(m)
        if len(auts) >= cap: break
    return auts

def dedupe_under_aut(colorings, auts):
    def canon(p): return tuple(sorted(tuple(sorted(cl)) for cl in p))
    seen, reps = set(), []
    for p in colorings:
        orbit = min(canon([[ekey(a[u], a[v]) for (u, v) in cl] for cl in p])
                    for a in auts)
        if orbit not in seen:
            seen.add(orbit); reps.append(sorted([sorted(cl) for cl in p]))
    return reps

# ============================================================
#  Layouts + plotting
# ============================================================

def _ring(nodes, r, offset=0.0):
    m = len(nodes)
    return {v: (r*np.cos(2*np.pi*i/m + offset), r*np.sin(2*np.pi*i/m + offset))
            for i, v in enumerate(nodes)}

def canonical_layout(key, G, n):
    if key in ("cycle", "complete"):
        return _ring(list(range(n)), 1.0)
    if key == "star":
        pos = {0: (0.0, 0.0)}; pos.update(_ring(list(range(1, n+1)), 1.0)); return pos
    if key == "mobius":
        return _ring(list(range(2*n)), 1.0)
    if key == "sunlet":
        pos = _ring(list(range(n)), 1.0); pos.update(_ring(list(range(n, 2*n)), 1.7)); return pos
    if key == "sun":
        pos = _ring(list(range(n)), 1.0)
        pos.update(_ring(list(range(n, 2*n)), 1.8, offset=np.pi/n)); return pos
    if key == "circular_ladder":
        pos = _ring(list(range(n)), 1.0); pos.update(_ring(list(range(n, 2*n)), 1.8)); return pos
    if key == "web":
        pos = _ring(list(range(n)), 1.0)
        pos.update(_ring(list(range(n, 2*n)), 1.8))
        pos.update(_ring(list(range(2*n, 3*n)), 2.6)); return pos
    if key == "ladder":
        return {**{i: (i, 0.0) for i in range(n)}, **{n+i: (i, 1.0) for i in range(n)}}
    if key == "kayak":
        if not isinstance(n, (tuple, list)):
            return nx.kamada_kawai_layout(G)          # scalar caller: no params to place
        c1, c2, l = kayak_params(n)
        pos = {v: (x - 2.4, y) for v, (x, y) in _ring(list(range(c1)), 1.0).items()}
        pos.update({v: (x + 2.4, y) for v, (x, y) in _ring(list(range(c1, c1+c2)), 1.0).items()})
        for j in range(l - 1):                        # shaft verts strung between blades
            pos[c1 + c2 + j] = (-1.4 + (j + 1) * 2.8 / l, 0.0)
        return pos
    return nx.kamada_kawai_layout(G)

PALETTE = plt.cm.tab10.colors

def plot_colorings(key, G, n, colorings, chi, save=None, max_plots=12):
    shown = colorings[:max_plots]
    cols = math.ceil(math.sqrt(len(shown))); rows = math.ceil(len(shown)/cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 4*rows), squeeze=False)
    for ax in axes.flat: ax.axis("off")
    pos = canonical_layout(key, G, n)
    for ax, part in zip(axes.flat, shown):
        for ci, cl in enumerate(part):
            nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[tuple(e) for e in cl],
                                   edge_color=[PALETTE[ci % 10]], width=2.4)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=80, node_color="#333")
    fig.suptitle(f"{key} n={n}: bond-edge colorings, beta={chi}  "
                 f"(showing {len(shown)}/{len(colorings)})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if save: fig.savefig(save, dpi=130); print(f"saved {save}")
    else: plt.show()

# ============================================================
#  Reporting
# ============================================================

def report(key, G, n, quiet=False):
    targets, H, ambiguous = analyze(G)
    if targets == 0:
        return dict(targets=0, verdict="no targets", beta=1, n_zero=0,
                    checks=0, free=0, reps=[], conflicts=0)
    if quiet:
        st = lazy_refine(G, H.copy(), dict(ambiguous))     # true beta, minimal checks
        verdict = ("unswappable" if not ambiguous else
                   ("∃ swap" if st["free_swaps"] else "unswappable*"))
        return dict(targets=targets, verdict=verdict, beta=st["chi"],
                    n_zero=len(ambiguous), checks=st["checks"],
                    free=st["free_swaps"], reps=None, conflicts=H.number_of_edges())
    Hstar, st, _ = complete_conflicts(G, H, ambiguous)      # full resolution
    beta, colorings = chromatic_and_colorings(Hstar)
    reps = dedupe_under_aut(colorings, automorphisms(G)) if colorings else []
    if not ambiguous:
        verdict = "unswappable"
    elif st["free_swaps"]:
        verdict = "∃ swap"
    else:
        verdict = "unswappable"        # all Delta==0 switches proved non-iso too
    return dict(targets=targets, verdict=verdict, beta=max(beta, 1),
                n_zero=len(ambiguous), checks=st["checks"], free=st["free_swaps"],
                reps=reps, conflicts=Hstar.number_of_edges())

# ============================================================
#  CLI
# ============================================================

def ask(p):
    try: return input(p).strip()
    except EOFError: return ""

def parse_range(s):
    p = s.replace(",", " ").split()
    if len(p) == 1: return range(int(p[0]), int(p[0]) + 1)
    if len(p) == 2: return range(int(p[0]), int(p[1]) + 1)
    return range(int(p[0]), int(p[1]) + 1, int(p[2]))

def parse_tuples(s):
    """Parse whitespace-separated tuples like '(3,3,1) (4,5,2)'. A bare
    '3,3,1' with no parentheses is accepted as one tuple. Returns int tuples."""
    import re
    groups = re.findall(r"\(([^)]*)\)", s)
    if not groups and s.strip():
        groups = [s]
    out = []
    for g in groups:
        vals = [int(x) for x in g.replace(",", " ").split()]
        if vals:
            out.append(tuple(vals))
    return out

def kayak_valid_params(raw_tuples, minn):
    """Normalize raw kayak tuples to (c1, l, c2), warning on and dropping any
    out of range (a cycle below minn, or a path length below 1)."""
    out = []
    for t in raw_tuples:
        p = kayak_params(t)                       # (c1, c2, l)
        if p[0] >= minn and p[1] >= minn and p[2] >= 1:
            out.append(p)
        else:
            print(f"  skipping {p}: need cycles c1,c2 >= {minn} and path l >= 1")
    return out

def parse_custom_graph(s):
    G = nx.Graph()
    for tok in s.replace(",", " ").split():
        u, v = tok.split("-"); G.add_edge(int(u), int(v))
    return G

def main():
    names = list(FAMILY_GENS)
    print("families (n = symmetry index):")
    for i, nm in enumerate(names, 1):
        print(f"  {i:>2}. {nm:<16} {FAMILY_GENS[nm][1]}")
    print(f"  {len(names)+1:>2}. custom           paste an edge list")

    sel = ask("\nselect family [number or name]: ")
    if sel == str(len(names) + 1) or sel == "custom":
        G = parse_custom_graph(ask("edge list (e.g. 0-1 1-2 2-0): "))
        key, items = "custom", [(G.number_of_nodes(), G)]
    else:
        key = names[int(sel) - 1] if sel.isdigit() and 1 <= int(sel) <= len(names) else sel
        if key not in FAMILY_GENS: print("unknown family."); return
        gen, desc, minn = FAMILY_GENS[key]
        if key in TUPLE_FAMILIES:
            raw = parse_tuples(ask(f"parameters for '{key}' as tuples "
                                   f"(cycle,cycle,path), e.g. (3,4,2)  [min cycle {minn}]: "))
            items = [(p, gen(*p)) for p in kayak_valid_params(raw, minn)]
        else:
            rng = parse_range(ask(f"n range for '{key}' (min {minn}) -- 'start end [step]': "))
            items = [(n, gen(n)) for n in rng if n >= minn]
        print(f"\n{key}: {desc}")

    if not items:
        print("no valid graphs to analyze."); return

    quiet = ask("suppress graph output, bond-edge number only? [y/N] ").lower() in ("y", "yes")

    if quiet:
        print(f"\n{'n':>9}{'verdict':>14}{'beta':>6}{'Δ=0':>6}{'checks':>8}")
        print("-" * 43)
        for n, G in items:
            r = report(key, G, n, quiet=True)
            print(f"{str(n):>9}{r['verdict']:>14}{r['beta']:>6}{r['n_zero']:>6}{r['checks']:>8}")
        print("\nbeta = bond-edge number chi(H*), refined by brute-forcing only the"
              "\nDelta==0 pairs that a minimum coloring places in one class."
              "\nunswappable* = beta certified but no free swap surfaced (lazy).")
        return

    print(f"\n{'n':>9}{'verts':>7}{'edges':>7}{'targets':>9}{'Δ=0':>6}"
          f"{'checks':>8}{'conflicts':>11}{'verdict':>13}{'beta':>6}{'#col':>7}")
    print("-" * 85)
    results = {}
    for n, G in items:
        r = report(key, G, n); results[n] = (G, r)
        print(f"{str(n):>9}{G.number_of_nodes():>7}{G.number_of_edges():>7}{r['targets']:>9}"
              f"{r['n_zero']:>6}{r['checks']:>8}{r['conflicts']:>11}{r['verdict']:>13}"
              f"{r['beta']:>6}{len(r['reps']):>7}")
    print("\nΔ=0: targets the triangle test could not decide."
          "\nchecks: brute-force switch tests actually run (only on Δ=0 pairs)."
          "\nconflicts: edges of the true conflict graph H* (Δ!=0 plus Δ=0-non-iso)."
          "\nbeta = chi(H*); #col = minimum colorings up to relabeling and Aut(G).")

    prompt = ("\nplot colorings for (k,m,l) = ? (blank to quit): " if key in TUPLE_FAMILIES
              else "\nplot colorings for n = ? (blank to quit): ")
    while True:
        s = ask(prompt)
        if not s: break
        if key in TUPLE_FAMILIES:
            tp = parse_tuples(s)
            n = kayak_params(tp[0]) if tp else None
        else:
            n = int(s)
        if n not in results: print("not computed."); continue
        G, r = results[n]
        if not r["reps"]: print("no colorings to plot."); continue
        plot_colorings(key, G, n, r["reps"], r["beta"])

if __name__ == "__main__":
    main()