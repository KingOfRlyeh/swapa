"""
b3bound.py -- menu-driven optimal-pot explorer for B_3 bond-edge types.

Run it, pick a graph, then use the menu to see colorings, orientations, the
lower bound, failing construction matrices, and (opt-in) isomorphism validation
and recoloring. Heavy isomorphism work is always behind a prompt.

    python3 b3bound.py                       interactive
    python3 b3bound.py --family cycle --range 4 6
    python3 b3bound.py --edges "0-1 1-2 2-0 2-3"
"""
import itertools, os
from fractions import Fraction
from math import gcd

import networkx as nx
import swapcolor as sc

# ============================================================
#  Stage 1: triangle method (colorings)
# ============================================================

def triangle_colorings(G, aut_cap=20000):
    """Full conflict graph and its minimum colorings, matching swapcolor.

    The triangle test (Delta != 0) certifies distinct bond-edge types, but is
    SILENT when Delta == 0 -- a swap can preserve triangle count yet still change
    the graph (e.g. alter 4-cycle count or disconnect it). Those Delta == 0 pairs
    are resolved by an isomorphism check (sc.complete_conflicts), exactly as
    swapcolor does, so the coloring here has no residual swap failures.

    Also returns safe_switches: for every pair the analysis found FREE (allowed
    to share a bond type) despite having at least one non-isomorphic switch, the
    list of that pair's isomorphic (v,w,s,t) reconnections -- i.e. the only
    orientations that pair may realize. Graph-level, independent of which
    coloring/orientation is later chosen; lets stage 3 reject a doomed
    orientation with a cheap lookup instead of a fresh isomorphism check."""
    _, H, ambiguous = sc.analyze(G)
    Hstar, _, safe_switches = sc.complete_conflicts(G, H, ambiguous)   # resolve Delta==0 pairs
    k, parts = sc.chromatic_and_colorings(Hstar)
    auts = sc.automorphisms(G, cap=aut_cap)
    reps = sc.dedupe_under_aut(parts, auts) if parts else []
    colorings = []
    for part in reps:
        cd = {}
        for i, cl in enumerate(part):
            for e in cl:
                cd[sc.ekey(*e)] = i
        colorings.append(cd)
    return k, colorings, Hstar, safe_switches, auts

# ============================================================
#  Stage 2: orientations (loop rule)
# ============================================================

def _coloring_stabilizer(G, coloring, auts):
    """The subgroup of `auts` (automorphisms of G) that maps this coloring's
    color-class PARTITION onto itself -- as a set of edge-sets, ignoring which
    numeric color id lands where, since bond-type labels are arbitrary. Always
    includes the identity."""
    class_sets = [frozenset(cl) for cl in color_classes(G, coloring).values()]
    stab = []
    for a in auts:
        ok = True
        for cl in class_sets:
            img = frozenset(sc.ekey(a[u], a[v]) for (u, v) in cl)
            if img not in class_sets:
                ok = False; break
        if ok:
            stab.append(a)
    return stab


def _component_action(a, comps):
    """How stabilizer automorphism `a` acts on the free-component flip-vector:
    component i's chosen side maps onto component perm[i]'s side, XORed with
    flipbits[i] (0 = preserves which side is the source, 1 = swaps it). Well
    defined because a connected bipartite graph's 2-coloring is unique up to a
    single global swap, so one sample vertex per component fixes the bit for
    the whole component."""
    edge_sets = [frozenset(ce) for ce, _ in comps]
    perm = [None] * len(comps)
    flipbits = [0] * len(comps)
    for i, (cedges, side) in enumerate(comps):
        image = frozenset(sc.ekey(a[u], a[v]) for (u, v) in cedges)
        j = edge_sets.index(image)
        perm[i] = j
        v = cedges[0][0]
        flipbits[i] = 0 if side[v] == comps[j][1][a[v]] else 1
    return perm, flipbits


def _flip_orbits(m, actions):
    """Orbit representatives of {0,1}^m (as bit-tuples) under the group action
    given by `actions` (a list of (perm, flipbits) pairs, one per generator)."""
    total = 1 << m
    visited = [False] * total
    reps = []
    for start in range(total):
        if visited[start]:
            continue
        reps.append(start)
        visited[start] = True
        frontier = [start]
        while frontier:
            cur = frontier.pop()
            bits = [(cur >> k) & 1 for k in range(m)]
            for perm, flipbits in actions:
                new_bits = [0] * m
                for i in range(m):
                    new_bits[perm[i]] = bits[i] ^ flipbits[i]
                val = 0
                for k, b in enumerate(new_bits):
                    val |= b << k
                if not visited[val]:
                    visited[val] = True
                    frontier.append(val)
    return [tuple((r >> k) & 1 for k in range(m)) for r in reps]


def orientations(G, coloring, auts=None):
    """Yield loop-valid orientations, one representative per pot-equivalence class.

    Loop rule: each color class must be bipartite; each connected component has
    two directions. But a SATURATED class -- all its edges share one vertex (a
    star) -- has both directions equivalent: flipping the star is just the
    a<->a' relabel of that bond type, a non-distinction. So a saturated class is
    pinned (arrows point OUT of the center) and contributes x1, not x2. This
    removes a factor of 2 per saturated class from the orientation count.

    If `auts` (automorphisms of G, e.g. from swapcolor.automorphisms) is given,
    a further reduction applies: this coloring's STABILIZER within `auts` (the
    automorphisms mapping its color-class partition onto itself) acts on the
    free-component flip-vector space, and flip-vectors in the same orbit yield
    orientations related by an actual automorphism of G -- so they succeed or
    fail every later check identically and produce isomorphic pots. Only one
    representative per orbit is yielded. Safe to omit (auts=None): falls back
    to enumerating every flip combination, as before."""
    by_color = {}
    for e in G.edges():
        k = sc.ekey(*e); by_color.setdefault(coloring[k], []).append(k)
    fixed = {}          # pinned edges (saturated classes)
    comps = []          # (component_edges, side) for classes that keep their freedom
    for c, edges in by_color.items():
        Gc = nx.Graph(); Gc.add_edges_from([tuple(e) for e in edges])
        try:
            side = nx.bipartite.color(Gc)
        except nx.NetworkXError:
            return                                    # odd cycle: no loop-free orientation
        center = saturated_vertex(edges)              # common to ALL edges of the class?
        if center is None and len(edges) == 1:
            center = edges[0][0]                       # single edge: either endpoint sources it
        if center is not None:
            for (u, v) in edges:                      # pin: center is the tail (arrow out)
                fixed[sc.ekey(u, v)] = (u, v) if u == center else (v, u)
            continue
        for comp in nx.connected_components(Gc):
            sub = Gc.subgraph(comp)
            comps.append(([sc.ekey(u, v) for u, v in sub.edges()], side))

    m = len(comps)
    flip_tuples = None
    if auts and m > 0:
        stab = _coloring_stabilizer(G, coloring, auts)
        if len(stab) > 1:                              # nontrivial: worth the BFS
            actions = [_component_action(a, comps) for a in stab]
            flip_tuples = _flip_orbits(m, actions)
    if flip_tuples is None:
        flip_tuples = list(itertools.product((0, 1), repeat=m))

    for flips in flip_tuples:
        orient = dict(fixed)
        for (cedges, side), flip in zip(comps, flips):
            for (u, v) in cedges:
                orient[sc.ekey(u, v)] = (u, v) if side[u] == flip else (v, u)
        yield orient

# ============================================================
#  Stage 3: multiedge rule
# ============================================================

def multiedge_offenders(G, coloring, orient):
    """ALL nonadjacent same-type arrow pairs whose head-swap would form a
    multiedge (completely enumerated -- not stopped at the first). Returns a list
    of offending edge pairs (e1, e2)."""
    Eset = {frozenset(e) for e in G.edges()}
    by_color = {}
    for e in G.edges():
        k = sc.ekey(*e); by_color.setdefault(coloring[k], []).append(orient[k])
    out = []
    for c, arcs in by_color.items():
        for (t1, h1), (t2, h2) in itertools.combinations(arcs, 2):
            if len({t1, h1, t2, h2}) < 4:
                continue
            if frozenset((t1, h2)) in Eset or frozenset((t2, h1)) in Eset:
                out.append((sc.ekey(t1, h1), sc.ekey(t2, h2)))
    return out

# ============================================================
#  Stage 3b: precomputed non-isomorphism rule (cheap, no fresh iso check)
# ============================================================

def iso_offenders(G, coloring, orient, safe_switches):
    """ALL nonadjacent same-type arrow pairs whose realized head-swap is a
    CERTIFIED non-isomorphism, per the graph-level Delta/isomorphism analysis
    (sc.complete_conflicts). A pair only ends up in safe_switches when it has
    at least one non-isomorphic switch alongside its isomorphic one(s); an
    orientation is doomed exactly when it realizes one of the unsafe switches
    instead. This is a lookup against work already done once per graph -- no
    fresh nx.is_isomorphic call -- so it lets stage 3 reject orientations that
    would otherwise only be caught later (and much more expensively) by
    validate_candidate. Pairs absent from safe_switches need no check here:
    either every switch was Delta != 0 (forced apart, never same-colored) or
    both possible reconnections multiedge -- already caught above."""
    by_color = {}
    for e in G.edges():
        k = sc.ekey(*e); by_color.setdefault(coloring[k], []).append(orient[k])
    out = []
    for c, arcs in by_color.items():
        for (t1, h1), (t2, h2) in itertools.combinations(arcs, 2):
            if len({t1, h1, t2, h2}) < 4:
                continue
            pair = frozenset((sc.ekey(t1, h1), sc.ekey(t2, h2)))
            safe = safe_switches.get(pair)
            if safe is None:
                continue
            safe_edges = {frozenset((sc.ekey(v, t), sc.ekey(s, w))) for (v, w, s, t) in safe}
            realized = frozenset((sc.ekey(t1, h2), sc.ekey(t2, h1)))
            if realized not in safe_edges:
                out.append((sc.ekey(t1, h1), sc.ekey(t2, h2)))
    return out


def color_classes(G, coloring):
    by_color = {}
    for e in G.edges():
        k = sc.ekey(*e); by_color.setdefault(coloring[k], []).append(k)
    return by_color


def saturated_vertex(edges):
    """The unique vertex incident to EVERY edge of the class (the star center),
    or None if the class is not a single star."""
    common = set(edges[0])
    for e in edges[1:]:
        common &= set(e)
    return next(iter(common)) if len(common) == 1 else None


def outlier_edges(G, coloring):
    """For each color class that is 'a star plus one stray edge', the stray edge
    -- the one whose removal leaves the rest sharing a common vertex (saturated).
    Returns list of (color, edge). These are the cheapest promotions to try."""
    out = []
    for c, edges in color_classes(G, coloring).items():
        if len(edges) < 2:
            continue
        if saturated_vertex(edges) is not None:
            continue                                     # already saturated: no outlier
        for i, cand in enumerate(edges):
            rest = edges[:i] + edges[i + 1:]
            if rest and saturated_vertex(rest) is not None:
                out.append((c, cand)); break
    return out

# ============================================================
#  Stage 4: construction matrix (m_P)
# ============================================================

def _lcm(a, b): return a * b // gcd(a, b)

def _solve_exact(A, b, nvars):
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    rows_ = len(M); r = 0; pivots = []
    for c in range(nvars):
        piv = next((i for i in range(r, rows_) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for i in range(rows_):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        pivots.append(c); r += 1
        if r == rows_: break
    for i in range(rows_):
        if all(M[i][j] == 0 for j in range(nvars)) and M[i][nvars] != 0:
            return None
    if len(pivots) < nvars: return None
    x = [Fraction(0)] * nvars
    for i, c in enumerate(pivots): x[c] = M[i][nvars]
    return x

def build_pot_Z(G, coloring, orient):
    colors = sorted(set(coloring.values()))
    vnet = {v: {c: 0 for c in colors} for v in G.nodes()}
    vmult = {v: [] for v in G.nodes()}
    for e in G.edges():
        k = sc.ekey(*e); c = coloring[k]; tail, head = orient[k]
        vnet[tail][c] += 1; vmult[tail].append((c, +1))
        vnet[head][c] -= 1; vmult[head].append((c, -1))
    sig_net = {}
    for v in G.nodes():
        sig = tuple(sorted(vmult[v]))
        sig_net[sig] = tuple(vnet[v][c] for c in colors)
    tiles = list(sig_net.keys())
    Z = [[sig_net[t][i] for t in tiles] for i in range(len(colors))]
    return Z, tiles, colors

def min_realizable_size(Z, p, n_cap):
    """m_P: the size of the smallest complex the pot can build = the smallest
    nonnegative INTEGER tile multiset n != 0 with Z n = 0 (every cohesive end
    matched). We search m = 1, 2, ... for a nonneg integer null vector summing
    to m; the first hit is m_P. (All-ones is always a null vector -- it is G
    itself -- so m_P <= n_cap = |V(G)| always.)

    This replaces an earlier vertex-only search that considered only supports
    with a UNIQUE solution; that skipped complexes on underdetermined supports --
    e.g. G when its tiles are all distinct -- and could report m_P > n."""
    rows = Z
    nrows = len(rows)

    def _search(m):
        # nonneg integer vectors of length p summing to m with Z n = 0, with
        # partial-sum pruning on each bond-type row.
        result = []
        counts = [0] * p

        def rec(j, remaining, partial):
            if j == p:
                if remaining == 0 and all(partial[i] == 0 for i in range(nrows)):
                    result.append(tuple(counts))
                    return True
                return False
            for v in range(remaining + 1):
                counts[j] = v
                np_ = [partial[i] + rows[i][j] * v for i in range(nrows)]
                if rec(j + 1, remaining - v, np_):
                    return True
            counts[j] = 0
            return False

        return result if rec(0, m, [0] * nrows) else None

    for m in range(1, n_cap + 1):
        hit = _search(m)
        if hit:
            n = hit[0]
            supp = tuple(j for j in range(p) if n[j] > 0)
            sol = [Fraction(n[j], m) for j in supp]
            return m, sol, supp
    return n_cap, None, None

# ============================================================
#  Stage 5: isomorphism validation (opt-in)
# ============================================================

def _wl(H):
    """Weisfeiler-Lehman hash, warning silenced (we only ever compare hashes
    computed within a single run/version, so the v3.5 reproducibility note is
    irrelevant here)."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return nx.weisfeiler_lehman_graph_hash(H)

def _all_null_vectors(Z, n):
    """Every nonnegative integer tile-count vector R (len = #tiles) with
    sum R = n and Z R = 0 -- i.e. every way the pot assembles a complete complex
    of order n. (min_realizable_size finds only the first; a valid Scenario-3
    pot must have this be unique up to the pairings below giving only G.)"""
    p = len(Z[0]); nrows = len(Z); sols = []; counts = [0] * p

    def rec(j, remaining, partial):
        if j == p:
            if remaining == 0 and all(x == 0 for x in partial):
                sols.append(tuple(counts))
            return
        for v in range(remaining + 1):
            counts[j] = v
            rec(j + 1, remaining - v, [partial[i] + Z[i][j] * v for i in range(nrows)])
        counts[j] = 0

    rec(0, n, [0] * nrows)
    return sols

def _distinct_matchings(tails, heads):
    """Distinct ways to bond one bond type's unhatted ends (`tails`, tile-instance
    ids with multiplicity) to its hatted ends (`heads`). Two matchings that only
    permute ends on the SAME tile give the same complex, so identical tails are
    forced to take non-decreasing heads -- that dedup keeps stars from blowing up
    to k! copies of one graph."""
    tails = sorted(tails)

    def rec(i, remaining, acc):
        if i == len(tails):
            yield tuple(acc); return
        t = tails[i]
        lo = acc[-1][1] if (i > 0 and tails[i - 1] == t) else None
        for h in sorted(set(remaining)):
            if lo is not None and h < lo:
                continue
            rem = list(remaining); rem.remove(h)
            yield from rec(i + 1, rem, acc + [(t, h)])

    yield from rec(0, list(heads), [])

def _pot_complexes(tiles, colors, R):
    """Every complete complex the pot assembles from the tile multiset R: the
    product over bond types of that type's distinct end-pairings. Yields
    nx.MultiGraph on the tile instances (loops/parallel edges may appear)."""
    inst = []
    for sig, cnt in zip(tiles, R):
        inst.extend([sig] * cnt)
    tails = {c: [] for c in colors}; heads = {c: [] for c in colors}
    for idx, sig in enumerate(inst):
        for (c, s) in sig:
            (tails if s > 0 else heads)[c].append(idx)
    per_color = [list(_distinct_matchings(tails[c], heads[c])) for c in colors]
    for combo in itertools.product(*per_color):
        M = nx.MultiGraph(); M.add_nodes_from(range(len(inst)))
        for cm in combo:
            for (t, h) in cm:
                M.add_edge(t, h)
        yield M

def validate_candidate(G, coloring, orient, cap=500000):
    """Complete Scenario-3 certificate for one candidate pot: EVERY complex it can
    assemble at order n = |V(G)| is isomorphic to G, so {G} = C_min(P) (given the
    caller already has m_P = n). This is stronger than the old swap-only test: a
    pot can have m_P = n and pass every single 2-switch yet still build a
    non-isomorphic order-n complex from a DIFFERENT tile-count mix (e.g. a second
    integer null vector of the construction matrix). We therefore enumerate all
    tile-count vectors and all end-pairings.

    Returns True (certified), False (a non-isomorphic / looped / multi-edged
    complex exists), or None (enumeration passed `cap` complexes undecided --
    treat as unverified, not as pass)."""
    # Fast necessary pre-filter: G's own realization must be swap-stable (Lemma 4).
    # A failing single 2-switch is already a non-isomorphic order-n complex.
    by_color = {}
    for e in G.edges():
        k = sc.ekey(*e); by_color.setdefault(coloring[k], []).append(orient[k])
    for c, arcs in by_color.items():
        for (t1, h1), (t2, h2) in itertools.combinations(arcs, 2):
            if len({t1, h1, t2, h2}) < 4: continue
            H = G.copy()
            H.remove_edge(t1, h1); H.remove_edge(t2, h2)
            H.add_edge(t1, h2); H.add_edge(t2, h1)
            if not nx.is_isomorphic(G, H):
                return False

    # Complete check: every order-n complex the pot can build must be iso to G.
    n = G.number_of_nodes()
    Z, tiles, colors = build_pot_Z(G, coloring, orient)
    ghash = _wl(G)
    seen = 0
    for R in _all_null_vectors(Z, n):
        for M in _pot_complexes(tiles, colors, R):
            seen += 1
            if seen > cap:
                return None
            if any(u == v for u, v in M.edges()):
                return False                                   # loop -> not iso to loopless G
            S = nx.Graph(M)
            if S.number_of_edges() != M.number_of_edges():
                return False                                   # multiedge
            if _wl(S) != ghash:
                return False                                   # cheap certificate mismatch
            if not nx.is_isomorphic(G, S):
                return False                                   # WL-collision guard
    return True

# ============================================================
#  Full analysis of one graph (stages 1-4, cached)
# ============================================================

class _Progress:
    """Minimal CLI progress bar (no dependencies). Writes to stderr with \\r."""
    def __init__(self, total, label="working", width=30):
        import time
        self.total = max(total, 1); self.label = label; self.width = width
        self.i = 0; self.start = time.time(); self.active = total > 1

    def tick(self, step=1):
        self.i += step
        if not self.active:
            return
        import sys, time
        frac = min(self.i / self.total, 1.0)
        fill = int(self.width * frac)
        bar = "#" * fill + "-" * (self.width - fill)
        sys.stderr.write(f"\r  [{bar}] {self.i}/{self.total} {self.label} "
                         f"({time.time() - self.start:.1f}s)")
        sys.stderr.flush()

    def close(self):
        if self.active:
            import sys
            sys.stderr.write("\n"); sys.stderr.flush()


def _analyze_colorings(G, colorings, safe_switches, auts=None, show=False):
    n = G.number_of_nodes()
    jobs = [(coloring, list(orientations(G, coloring, auts))) for coloring in colorings]
    prog = _Progress(sum(len(o) for _, o in jobs), "orientations") if show else None
    per = []
    for coloring, orients in jobs:
        cands, mat_fails, me_fails, iso_fails = [], [], [], []
        for orient in orients:
            if prog:
                prog.tick()
            offenders = multiedge_offenders(G, coloring, orient)      # ALL of them
            if offenders:
                me_fails.append((orient, offenders)); continue
            offenders = iso_offenders(G, coloring, orient, safe_switches)  # cheap, precomputed
            if offenders:
                iso_fails.append((orient, offenders)); continue
            Z, tiles, colors = build_pot_Z(G, coloring, orient)
            mp, sol, supp = min_realizable_size(Z, len(tiles), n)
            if mp == n:
                cands.append((orient, Z, tiles, colors))
            else:
                mat_fails.append((orient, Z, tiles, colors, mp, sol, supp))
        per.append(dict(coloring=coloring, orients=orients, candidates=cands,
                        matrix_fails=mat_fails, multiedge_fails=me_fails, iso_fails=iso_fails))
    if prog:
        prog.close()
    return per


def analyze_graph(G):
    n = G.number_of_nodes()
    k, colorings, Hstar, safe_switches, auts = triangle_colorings(G)
    print("analyzing orientations...")
    per = _analyze_colorings(G, colorings, safe_switches, auts, show=True)
    ncand = sum(len(p["candidates"]) for p in per)
    return dict(n=n, k=k, colorings=colorings, per=per, Hstar=Hstar,
                safe_switches=safe_switches, auts=auts,
                proven=(k if ncand > 0 else k + 1), partial=False)

# ============================================================
#  Rendering + text
# ============================================================

PALETTE = None
def _pal():
    global PALETTE
    if PALETTE is None:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        PALETTE = plt.cm.tab10.colors
    return PALETTE

def _layout(G, family_key):
    n = G.number_of_nodes()
    if family_key in sc.FAMILY_GENS:
        # kayak is tuple-parameterized, not derivable from vertex count -> kamada
        fam_n = {"cycle": n, "complete": n, "star": n - 1, "mobius": n // 2,
                 "sunlet": n // 2, "sun": n // 2, "circular_ladder": n // 2,
                 "ladder": n // 2, "web": n // 3}.get(family_key, n)
        try:
            return sc.canonical_layout(family_key, G, fam_n)
        except Exception:
            pass
    return nx.kamada_kawai_layout(G)

def _letter(c): return chr(ord('a') + c)

def render_coloring(G, coloring, path, family_key=None, title=""):
    import matplotlib.pyplot as plt
    pal = _pal(); pos = _layout(G, family_key)
    fig, ax = plt.subplots(figsize=(4.5, 4.5)); ax.axis("off")
    for e in G.edges():
        k = sc.ekey(*e)
        nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[e],
                               edge_color=[pal[coloring[k] % 10]], width=2.4)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=90, node_color="#333")
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_color="white")
    ax.set_title(title); fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)

def render_orientation(G, coloring, orient, path, family_key=None, title=""):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    pal = _pal(); pos = _layout(G, family_key)
    fig, ax = plt.subplots(figsize=(4.5, 4.5)); ax.axis("off")
    D = nx.DiGraph(); D.add_nodes_from(G.nodes())
    for e in G.edges():
        k = sc.ekey(*e); t, h = orient[k]; D.add_edge(t, h, c=coloring[k])
    nx.draw_networkx_nodes(D, pos, ax=ax, node_size=90, node_color="#333")
    nx.draw_networkx_labels(D, pos, ax=ax, font_size=8, font_color="white")
    for t, h, d in D.edges(data=True):
        nx.draw_networkx_edges(D, pos, ax=ax, edgelist=[(t, h)],
                               edge_color=[pal[d["c"] % 10]], width=2.2,
                               arrows=True, arrowstyle="-|>", arrowsize=15, node_size=90)
    letters = sorted({coloring[sc.ekey(*e)] for e in G.edges()})
    ax.legend(handles=[Line2D([], [], color=pal[c % 10], lw=2.2, label=_letter(c))
                       for c in letters], loc="lower center", ncol=len(letters),
              frameon=False, fontsize=8)
    ax.set_title(title); fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)

def render_before_after(G, coloring, orient, offender, path, family_key=None, title=""):
    """Two side-by-side panels for one failing swap.
      LEFT  (before): the orientation; the two same-type edges are bold black.
      RIGHT (after):  those two edges removed and their heads re-paired; the new
                      edges are red, and the one landing on an existing edge is
                      the multiedge that rejects the orientation."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    pal = _pal(); pos = _layout(G, family_key)
    e1, e2 = offender
    t1, h1 = orient[e1]; t2, h2 = orient[e2]
    L = _letter(coloring[e1])
    Eset = {frozenset(e) for e in G.edges()}

    fig, (axb, axa) = plt.subplots(1, 2, figsize=(10, 5.2))
    for ax in (axb, axa):
        ax.axis("off")

    # ---- LEFT: before ----
    Db = nx.DiGraph(); Db.add_nodes_from(G.nodes())
    for e in G.edges():
        k = sc.ekey(*e); t, h = orient[k]; Db.add_edge(t, h, c=coloring[k])
    for t, h, d in Db.edges(data=True):
        bold = (t, h) in {(t1, h1), (t2, h2)}
        nx.draw_networkx_edges(Db, pos, ax=axb, edgelist=[(t, h)],
                               edge_color="black" if bold else [pal[d["c"] % 10]],
                               width=3.2 if bold else 1.8, arrows=True, arrowstyle="-|>",
                               arrowsize=16 if bold else 11, node_size=90,
                               alpha=1.0 if bold else 0.35)
    nx.draw_networkx_nodes(Db, pos, ax=axb, node_size=90, node_color="#333")
    nx.draw_networkx_labels(Db, pos, ax=axb, font_size=8, font_color="white")
    axb.set_title(f"before: type {L} edges  {t1}->{h1} , {t2}->{h2}", fontsize=10)

    # ---- RIGHT: after ----
    Da = nx.DiGraph(); Da.add_nodes_from(G.nodes())
    for e in G.edges():
        k = sc.ekey(*e)
        if k in (e1, e2):
            continue
        t, h = orient[k]; Da.add_edge(t, h, c=coloring[k])
    for t, h, d in Da.edges(data=True):
        nx.draw_networkx_edges(Da, pos, ax=axa, edgelist=[(t, h)], edge_color=[pal[d["c"] % 10]],
                               width=1.8, arrows=True, arrowstyle="-|>", arrowsize=11,
                               node_size=90, alpha=0.35)
    multi = None
    for (a, b) in [(t1, h2), (t2, h1)]:
        if frozenset((a, b)) in Eset:
            multi = (a, b)
        axa.annotate("", xy=pos[b], xytext=pos[a],
                     arrowprops=dict(arrowstyle="-|>", color="red", lw=2.8,
                                     connectionstyle="arc3,rad=0.3"))
    nx.draw_networkx_nodes(Da, pos, ax=axa, node_size=90, node_color="#333")
    nx.draw_networkx_labels(Da, pos, ax=axa, font_size=8, font_color="white")
    mtag = f"multiedge on {multi[0]}-{multi[1]}" if multi else "new edges"
    axa.set_title(f"after: add {t1}->{h2} , {t2}->{h1}   ({mtag})", fontsize=10)

    fig.legend(handles=[Line2D([], [], color="black", lw=3.2, label="the two same-type edges"),
                        Line2D([], [], color="red", lw=2.8, label="re-paired (multiedge)")],
               loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    fig.savefig(path, dpi=130); plt.close(fig)


def swap_text(orient, offender, coloring, G):
    e1, e2 = offender
    t1, h1 = orient[e1]; t2, h2 = orient[e2]
    Eset = {frozenset(e) for e in G.edges()}
    which = [f"{a}-{b}" for (a, b) in [(t1, h2), (t2, h1)] if frozenset((a, b)) in Eset]
    return (f"type {_letter(coloring[e1])}: remove {t1}->{h1}, {t2}->{h2}; "
            f"add {t1}->{h2}, {t2}->{h1}  ->  multiedge on {', '.join(which)}")


def pot_text(G, coloring, orient):
    ends = {v: [] for v in G.nodes()}; arcs = []
    for e in G.edges():
        k = sc.ekey(*e); L = _letter(coloring[k]); t, h = orient[k]
        ends[t].append(L); ends[h].append(L + "'"); arcs.append((t, h, L))
    tiles = {}
    for v in G.nodes():
        sig = tuple(sorted(ends[v])); tiles[sig] = tiles.get(sig, 0) + 1
    setf = "P = { " + "  +  ".join(f"{m} x {{{', '.join(s)}}}"
                                   for s, m in sorted(tiles.items())) + " }"
    return setf, arcs

def matrix_text(Z, tiles, colors, mp, n, sol=None, supp=None):
    lines = [f"m_P = {mp}   (need m_P = n = {n};  {'OK' if mp == n else 'FAILS: realizes a smaller complex'})",
             "tile types (columns T0..): " + ", ".join(
                 f"T{j}=" + "{" + ",".join(f"{_letter(c)}{'^' if s>0 else chr(39)}" for c, s in t) + "}"
                 for j, t in enumerate(tiles)),
             "Z (rows = bond types a,b,...; entries = #unhatted - #hatted):"]
    for i, c in enumerate(colors):
        lines.append(f"  {_letter(c)}: [" + " ".join(f"{Z[i][j]:>2}" for j in range(len(tiles))) + " ]")
    if sol is not None and supp is not None:
        terms = ", ".join(f"{s} * T{j}" for j, s in zip(supp, sol))
        lines.append(f"smaller-complex solution (Z r = 0, sum r = 1): r = {terms}")
        lines.append(f"  -> tile counts at size {mp}: " +
                     ", ".join(f"{int(s * mp)} x T{j}" for j, s in zip(supp, sol)))
        lines.append("  (only the support tiles are used; these proportions are the "
                     "extra constraint that pins a solution even when Z r = 0 is "
                     "underdetermined.)")
    return "\n".join(lines)

# ============================================================
#  Menu
# ============================================================

def _confirm(msg):
    try:
        return input(f"{msg} [y/N] ").strip().lower() in ("y", "yes")
    except EOFError:
        return False

def _ask(p):
    try:
        return input(p).strip()
    except EOFError:
        return ""

def explore(label, G, family_key=None):
    print(f"\nAnalyzing {label} ({G.number_of_nodes()} vertices, "
          f"{G.number_of_edges()} edges)...")
    A = analyze_graph(G)
    outdir = f"b3_{label.replace('(', '_').replace(')', '')}"
    os.makedirs(outdir, exist_ok=True)

    while True:
        k = A["k"]; ncand = sum(len(p["candidates"]) for p in A["per"])
        proven = A["proven"]
        print(f"\n===== {label} =====")
        if A.get("partial"):
            print(f"lower bound:  B_3 >= {proven} bond-edge types "
                  f"(recolored to k = {k}; partial coloring set; {ncand} candidate(s))")
        elif ncand > 0:
            print(f"lower bound:  B_3 >= {proven} bond-edge types "
                  f"({ncand} unvalidated candidate(s) at k = {k})")
        else:
            print(f"lower bound:  B_3 >= {proven} bond-edge types "
                  f"(no candidate at k = {k}: every orientation failed)")
        print(f"colorings: {len(A['colorings'])}   candidate pots (m_P=n): {ncand}"
              f"   output dir: {outdir}/")
        print("  1. show lower bound + summary")
        print("  2. visualize colorings")
        print("  3. visualize valid orientations")
        print("  4. dump failing construction matrices to file")
        print("  5. list candidate pots")
        print("  6. validate candidates  (isomorphism checks -- can be slow)")
        print("  7. try next color up     (exhausts the recolor pile -- can be slow)")
        print("  8. visualize edge-swap (multiedge) failures  [loop-valid, fail only the swap]")
        print("  0. done with this graph")
        ch = _ask("select: ")

        if ch == "0" or ch == "":
            return
        elif ch == "1":
            _summary(label, A)
        elif ch == "2":
            _viz_colorings(A, G, family_key, outdir)
        elif ch == "3":
            _viz_orientations(A, G, family_key, outdir)
        elif ch == "4":
            _dump_matrix_fails(A, outdir, label, G)
        elif ch == "5":
            _list_candidates(A, G)
        elif ch == "6":
            _validate(label, A, G, family_key, outdir)
        elif ch == "7":
            newA = _recolor(label, A, G, family_key, outdir)
            if newA is not None:
                A = newA          # advance: bound, k, colorings, viz now reflect k+1
                print(f"state advanced to k = {A['k']}. The menu now reflects this level.")
        elif ch == "8":
            _viz_multiedge_fails(A, G, family_key, outdir)
        else:
            print("?")

def _summary(label, A):
    k = A["k"]
    ncand = sum(len(p["candidates"]) for p in A["per"])
    nmefail = sum(len(p["multiedge_fails"]) for p in A["per"])
    nisofail = sum(len(p.get("iso_fails", [])) for p in A["per"])
    nmatfail = sum(len(p["matrix_fails"]) for p in A["per"])
    print(f"\nB_3({label}) >= {k}.")
    print(f"  {len(A['colorings'])} minimum coloring(s) at k = {k} types.")
    print(f"  candidate pots (passed loop + multiedge + iso + m_P=n): {ncand}")
    print(f"  orientations set aside by the multiedge rule:      {nmefail}")
    print(f"  orientations set aside by the precomputed iso rule: {nisofail}")
    print(f"  orientations rejected by the construction matrix:  {nmatfail}")
    if ncand == 0:
        print(f"  -> no candidate at k = {k}; the coloring must grow, so B_3 >= {k+1}.")
    else:
        print(f"  -> {ncand} candidate(s) exist; validate (menu 6) to confirm B_3 = {k}.")
 
def _viz_multiedge_fails(A, G, family_key, outdir):
    # one before/after per (orientation, offending swap)
    jobs = [(ci, orient, off) for ci, p in enumerate(A["per"], 1)
            for (orient, offs) in p["multiedge_fails"] for off in offs]
    if not jobs:
        print("no orientations fail only the edge-swap (multiedge) rule.")
        return
    if len(jobs) > 40 and not _confirm(f"{len(jobs)} before/after swaps to render; proceed?"):
        return
    prog = _Progress(len(jobs), "rendering")
    for i, (ci, orient, off) in enumerate(jobs, 1):
        prog.tick()
        coloring = A["per"][ci - 1]["coloring"]
        fn = os.path.join(outdir, f"edgeswap_c{ci}_{i}.png")
        render_before_after(G, coloring, orient, off, fn, family_key,
                            f"{outdir[3:]}  edge-swap failure {i}")
    prog.close()
    # also write the swaps as text (readable at a glance)
    txt = os.path.join(outdir, "edge_swaps.txt")
    with open(txt, "w", encoding="utf-8") as f:
        f.write("Loop-valid orientations that fail ONLY the edge-swap rule.\n"
                "Each line: a same-type head-swap that re-pairs into an existing "
                "edge (a multiedge), so the orientation is rejected.\n\n")
        for i, (ci, orient, off) in enumerate(jobs, 1):
            f.write(f"{i:>3}. (coloring {ci})  "
                    + swap_text(orient, off, A["per"][ci - 1]["coloring"], G) + "\n")
    print(f"{len(jobs)} before/after image(s) + {txt} written.")
    print("  left = orientation (swapped edges bold); right = after the swap "
          "(red = the re-paired multiedge).")


def _viz_colorings(A, G, family_key, outdir):
    for ci, coloring in enumerate(A["colorings"], 1):
        p = os.path.join(outdir, f"coloring_{ci}.png")
        render_coloring(G, coloring, p, family_key, f"coloring {ci} (k={A['k']})")
        print(f"  wrote {p}")
    print(f"{len(A['colorings'])} coloring image(s) written.")

def _viz_orientations(A, G, family_key, outdir):
    # a valid orientation passes the loop rule, the multiedge rule, AND the
    # construction matrix (m_P = n) -- i.e. exactly the candidate pots.
    valid = [(ci, orient) for ci, p in enumerate(A["per"], 1)
             for (orient, *_ ) in p["candidates"]]
    if not valid:
        print("no valid orientations at this level "
              "(none passed the multiedge rule and m_P = n).")
        return
    if len(valid) > 40 and not _confirm(f"{len(valid)} valid orientations to render; proceed?"):
        return
    for i, (ci, orient) in enumerate(valid, 1):
        fn = os.path.join(outdir, f"orient_c{ci}_{i}.png")
        render_orientation(G, A["per"][ci - 1]["coloring"], orient, fn, family_key,
                           f"coloring {ci}, valid orientation {i}")
    print(f"{len(valid)} valid orientation image(s) written to {outdir}/.")

def _dump_matrix_fails(A, outdir, label, G):
    path = os.path.join(outdir, "failing_construction_matrices.txt")
    n = A["n"]; total = 0
    with open(path, "w") as f:
        f.write(f"Failing construction matrices for {label} (m_P != n = {n}).\n")
        f.write("These orientations realize a SMALLER complex, so they are rejected.\n\n")
        for ci, p in enumerate(A["per"], 1):
            for oi, (orient, Z, tiles, colors, mp, sol, supp) in enumerate(p["matrix_fails"], 1):
                total += 1
                f.write(f"--- coloring {ci}, matrix-fail {oi} ---\n")
                arcs = "  ".join(f"{t}-{_letter(p['coloring'][sc.ekey(*e)])}->{h}"
                                 for e in G.edges()
                                 for (t, h) in [orient[sc.ekey(*e)]])
                f.write("orientation: " + arcs + "\n")
                f.write(matrix_text(Z, tiles, colors, mp, n, sol, supp) + "\n\n")
    print(f"{total} failing construction matrix/matrices written to {path}.")
# ha
def _list_candidates(A, G):
    tot = 0
    for ci, p in enumerate(A["per"], 1):
        for oi, (orient, Z, tiles, colors) in enumerate(p["candidates"], 1):
            tot += 1
            setf, _ = pot_text(G, p["coloring"], orient)
            print(f"  candidate {tot} (coloring {ci}): {setf}")
    if tot == 0:
        print("  no candidates at this level.")

def _validate(label, A, G, family_key, outdir):
    cands = [(ci, orient) for ci, p in enumerate(A["per"], 1)
             for (orient, *_ ) in p["candidates"]]
    if not cands:
        print("no candidates to validate."); return
    if not _confirm(f"validate {len(cands)} candidate(s) with isomorphism checks? "
                    "this can be expensive"):
        return
    found = 0; rejected = 0; undecided = 0
    prog = _Progress(len(cands), "validating")
    for ci, orient in cands:
        prog.tick()
        coloring = A["per"][ci - 1]["coloring"]
        res = validate_candidate(G, coloring, orient)
        if res is None:
            undecided += 1; continue
        if res is False:
            rejected += 1; continue                # builds a non-iso same-order complex
        found += 1
        setf, arcs = pot_text(G, coloring, orient)
        txt = os.path.join(outdir, f"optimal_tiling_{found}.txt")
        with open(txt, "w") as f:
            f.write(f"OPTIMAL TILING for {label}\n\n{setf}\n\noriented edges:\n" +
                    "\n".join(f"{t} -{L}-> {h}" for t, h, L in sorted(arcs)) + "\n")
        png = os.path.join(outdir, f"optimal_tiling_{found}.png")
        render_orientation(G, coloring, orient, png, family_key,
                           f"optimal tiling {found}: {label}")
        print(f"\n  OPTIMAL TILING {found}:  {setf}")
        print(f"    wrote {txt} , {png}")
    prog.close()
    if rejected:
        print(f"\n{rejected} candidate(s) rejected: they build a non-isomorphic complex "
              "of the same order (m_P = n is necessary but NOT sufficient for Scenario 3).")
    if undecided:
        print(f"{undecided} candidate(s) undecided (enumeration exceeded the cap).")
    if found:
        print(f"\n{found} optimal tiling(s) certified. B_3({label}) = {A['k']} achieved.")
    elif undecided and not rejected:
        print(f"\nno candidate certified within the cap; result inconclusive.")
    else:
        print(f"\nall candidates rejected. B_3({label}) >= {A['k'] + 1}.")

def _recolor(label, A, G, family_key, outdir, color_cap=1e6):
    """Advance one level, k -> k+1, by a COMPLETE enumeration of the proper
    (k+1)-colorings of the conflict graph (deduped under Aut(G)). A valid (k+1)
    coloring need NOT be a single-edge refinement of a k-coloring -- the prism is
    a counterexample -- so nothing less than the full enumeration is sound. Any
    candidate means B_3 could be k+1 (validate to confirm); an exhaustive
    enumeration with no candidate proves B_3 >= k+2. One level per call.
    Returns the advanced analysis dict, or None if declined."""
    k = A["k"]; Hstar = A.get("Hstar"); safe_switches = A.get("safe_switches")
    auts = A.get("auts")
    if Hstar is None or safe_switches is None:
        _, H0, amb = sc.analyze(G); Hstar, _, safe_switches = sc.complete_conflicts(G, H0, amb)
    if auts is None:
        auts = sc.automorphisms(G, cap=20000)
    if not _confirm(f"enumerate ALL proper {k+1}-colorings and analyze them "
                    f"(exhaustive; proves B_3 >= {k+2} if none works)? can be slow"):
        return None
    parts = sc.enumerate_colorings(Hstar, k + 1, cap=color_cap)
    complete = len(parts) < color_cap
    reps = sc.dedupe_under_aut(parts, auts) if parts else []
    colorings = []
    for part in reps:
        cd = {}
        for i, cl in enumerate(part):
            for e in cl:
                cd[sc.ekey(*e)] = i
        if len(set(cd.values())) == k + 1:          # actually uses k+1 types
            colorings.append(cd)
    print(f"analyzing {len(colorings)} coloring(s) at k = {k+1}...")
    per = _analyze_colorings(G, colorings, safe_switches, auts, show=True)
    ncand = sum(len(q["candidates"]) for q in per)
    if ncand > 0:
        proven = k + 1
        print(f"k = {k+1}: {ncand} candidate pot(s). Validate (menu 6) to confirm B_3 = {k+1}.")
    elif complete:
        proven = k + 2
        print(f"k = {k+1}: no candidate in the complete enumeration -> B_3 >= {k+2} PROVEN.")
    else:
        proven = k + 1
        print(f"k = {k+1}: no candidate, but enumeration hit the cap ({color_cap}); "
              "NOT proven -- raise the cap to certify.")
    return dict(n=A["n"], k=k + 1, colorings=colorings, per=per, Hstar=Hstar,
                safe_switches=safe_switches, auts=auts,
                proven=proven, partial=(not complete and ncand == 0))


def matrix_fail_edges(G, coloring, mf):
    """Edges implicated in a construction-matrix failure: those incident to the
    vertices whose tiles appear in the smaller-complex support. Recoloring one of
    these changes a used tile signature and can raise m_P ('increase the
    dimension')."""
    orient, Z, tiles, colors, mp, sol, supp = mf
    if not supp:
        return [sc.ekey(*e) for e in G.edges()]          # unknown: fall back to all
    used_sigs = {tiles[j] for j in supp}
    verts = []
    vmult = {v: [] for v in G.nodes()}
    for e in G.edges():
        k = sc.ekey(*e); c = coloring[k]; t, h = orient[k]
        vmult[t].append((c, +1)); vmult[h].append((c, -1))
    for v in G.nodes():
        if tuple(sorted(vmult[v])) in used_sigs:
            verts.append(v)
    vset = set(verts)
    return [sc.ekey(*e) for e in G.edges() if e[0] in vset or e[1] in vset]

# ============================================================
#  CLI
# ============================================================

def _parse_range(s):
    p = s.replace(",", " ").split()
    if len(p) == 1: return range(int(p[0]), int(p[0]) + 1)
    if len(p) == 2: return range(int(p[0]), int(p[1]) + 1)
    return range(int(p[0]), int(p[1]) + 1, int(p[2]))

def _parse_edges(s):
    G = nx.Graph()
    for tok in s.replace(",", " ").split():
        u, v = tok.split("-"); G.add_edge(int(u), int(v))
    return G

def main():
    import argparse
    fams = list(sc.FAMILY_GENS)
    ap = argparse.ArgumentParser(prog="b3bound",
        description="Menu-driven optimal-pot explorer for B_3 bond-edge types.")
    ap.add_argument("--family", choices=fams)
    ap.add_argument("--range", nargs="+", type=int, metavar="N")
    ap.add_argument("--edges")
    args = ap.parse_args()

    items = []; family_key = None
    if args.edges:
        items = [("custom", _parse_edges(args.edges))]
    elif args.family:
        family_key = args.family
        gen, desc, minn = sc.FAMILY_GENS[args.family]
        if args.family in sc.TUPLE_FAMILIES:
            vals = args.range or [minn, minn, 1]
            raw = [tuple(vals[i:i+3]) for i in range(0, len(vals), 3)]
            items = [(f"{args.family}({p[0]}-{p[1]}-{p[2]})", gen(*p))
                     for p in sc.kayak_valid_params(raw, minn)]
        else:
            rng = _parse_range(" ".join(map(str, args.range))) if args.range \
                else range(minn, minn + 1)
            items = [(f"{args.family}({n})", gen(n)) for n in rng if n >= minn]
    else:
        print("families (n = symmetry index):")
        for i, nm in enumerate(fams, 1):
            print(f"  {i:>2}. {nm:<16} {sc.FAMILY_GENS[nm][1]}")
        print(f"  {len(fams)+1:>2}. custom")
        sel = _ask("\nselect [number/name]: ")
        if sel in (str(len(fams) + 1), "custom"):
            items = [("custom", _parse_edges(_ask("edge list (e.g. 0-1 1-2 2-0): ")))]
        else:
            key = fams[int(sel) - 1] if sel.isdigit() and 1 <= int(sel) <= len(fams) else sel
            if key not in sc.FAMILY_GENS: print("unknown family."); return
            family_key = key
            gen, desc, minn = sc.FAMILY_GENS[key]
            if key in sc.TUPLE_FAMILIES:
                raw = sc.parse_tuples(_ask(f"parameters for '{key}' as tuples "
                                           f"(cycle,cycle,path), e.g. (3,4,2)  [min cycle {minn}]: "))
                items = [(f"{key}({p[0]}-{p[1]}-{p[2]})", gen(*p))
                         for p in sc.kayak_valid_params(raw, minn)]
            else:
                rng = _parse_range(_ask(f"n range for '{key}' (min {minn}): "))
                items = [(f"{key}({n})", gen(n)) for n in rng if n >= minn]

    if not items:
        print("no valid graphs to analyze."); return

    for label, G in items:
        explore(label, G, family_key)

if __name__ == "__main__":
    main()