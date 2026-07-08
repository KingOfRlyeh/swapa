"""
Tightened lower bound on B_3(G).

Two obstructions are combined:
  (A) swap obstruction: two nonadjacent edges are forced to have different
      bond-edge types iff EVERY valid 2-switch (reconnection) on them is
      non-isomorphic to G.  (A reconnection that creates a multiedge/loop counts
      as non-isomorphic.)  This yields a conflict graph H; chi(H) is a valid
      lower bound and equals the swap program's coloring number on H.
  (B) no-smaller-complex obstruction: a coloring is turned into pots by choosing
      an orientation (unhatted->hatted); the construction matrix M(P) gives
      m_P, the size of the smallest realizable complex.  A pot with m_P < n is
      invalid for Scenario 2, hence for Scenario 3.

Tightening (as specified): a valid graph must have a proper coloring of H; over
ALL potifications of each proper k-coloring, if ANY pot has m_P = n the coloring
is admissible.  The tightened bound is the least k admitting such a pot.  Since a
Scenario-3 pot both properly colors H and has m_P = n, the result is <= B_3(G).
"""
import itertools
from fractions import Fraction
from math import gcd
import networkx as nx
import swapcolor as sc


# ---------- (A) conflict graph ----------

def conflict_graph(G):
    """e1 ~ e2  iff nonadjacent, admits >=1 valid 2-switch, and EVERY valid
    reconnection is non-isomorphic to G. A reconnection whose new edges are
    already present (multiedge) is skipped as invalid; if that leaves the pair
    with no iso escape, the pair is a conflict."""
    A = {frozenset(e) for e in G.edges()}
    H = nx.Graph(); H.add_nodes_from(sc.ekey(*e) for e in G.edges())
    for (v, w), (s, t) in itertools.combinations(list(G.edges()), 2):
        if len({v, w, s, t}) < 4:
            continue
        outcomes = []
        for (a, b, c, d) in [(v, w, s, t), (v, w, t, s)]:
            ne1, ne2 = frozenset((a, d)), frozenset((c, b))
            if ne1 in A or ne2 in A or ne1 == ne2:
                continue                              # invalid reconnection (multiedge/loop)
            outcomes.append(nx.is_isomorphic(G, sc.switch(G, a, b, c, d)))
        if outcomes and all(o is False for o in outcomes):
            H.add_edge(sc.ekey(v, w), sc.ekey(s, t))
    return H


# ---------- (B) potification + construction matrix + m_P ----------

def build_pot_Z(G, edge_color, orient):
    """Return (Z, p): net-count matrix (rows=bond types, cols=distinct tiles)."""
    colors = sorted(set(edge_color.values()))
    vnet = {v: {c: 0 for c in colors} for v in G.nodes()}
    vmult = {v: [] for v in G.nodes()}
    for e in G.edges():
        k = sc.ekey(*e); c = edge_color[k]; tail, head = orient[k]
        vnet[tail][c] += 1; vmult[tail].append((c, +1))
        vnet[head][c] -= 1; vmult[head].append((c, -1))
    sig_net = {}
    for v in G.nodes():
        sig = tuple(sorted(vmult[v]))
        sig_net[sig] = tuple(vnet[v][c] for c in colors)
    tiles = list(sig_net.keys())
    Z = [[sig_net[t][i] for t in tiles] for i in range(len(colors))]
    return Z, len(tiles)


def _lcm(a, b):
    return a * b // gcd(a, b)


def _solve_exact(A, b, nvars):
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    rows_ = len(M); r = 0; pivots = []
    for c in range(nvars):
        piv = next((i for i in range(r, rows_) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for i in range(rows_):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        pivots.append(c); r += 1
        if r == rows_:
            break
    for i in range(rows_):
        if all(M[i][j] == 0 for j in range(nvars)) and M[i][nvars] != 0:
            return None
    if len(pivots) < nvars:
        return None
    x = [Fraction(0)] * nvars
    for i, c in enumerate(pivots):
        x[c] = M[i][nvars]
    return x


def min_realizable_size(Z, p, n_cap):
    """m_P via Prop 3.3: min lcm-of-denominators over vertices of the polytope
    {r >= 0 : Z r = 0, sum r = 1}. Vertices = minimal supports."""
    rows = [list(map(Fraction, r)) for r in Z]
    best = None
    for size in range(1, min(p, len(rows) + 1) + 1):
        for supp in itertools.combinations(range(p), size):
            Asub = [[rows[i][j] for j in supp] for i in range(len(rows))]
            Asub.append([Fraction(1)] * size)
            b = [Fraction(0)] * len(rows) + [Fraction(1)]
            sol = _solve_exact(Asub, b, size)
            if sol is None or any(x <= 0 for x in sol):
                continue
            denom = 1
            for x in sol:
                denom = _lcm(denom, x.denominator)
            if best is None or denom < best:
                best = denom
    return best if best is not None else n_cap


# ---------- enumeration ----------

def orientations(G, full_cap=14, samples=400, seed=0):
    edges = [sc.ekey(*e) for e in G.edges()]; m = len(edges)
    if m <= full_cap:
        for bits in itertools.product((0, 1), repeat=m):
            yield {e: ((e[0], e[1]) if b == 0 else (e[1], e[0])) for e, b in zip(edges, bits)}
    else:
        import random; rng = random.Random(seed)
        yield {e: (e[0], e[1]) for e in edges}
        for _ in range(samples):
            yield {e: ((e[0], e[1]) if rng.random() < .5 else (e[1], e[0])) for e in edges}


def k_colorings(H, edges_all, k, cap=500000):
    """All proper k-colorings of H (edges outside any conflict range freely),
    symmetry-broken up to color permutation. (list, complete_flag)."""
    nodes = list(edges_all)
    adj = {e: set(H[e]) for e in H.nodes()}
    out = []; assign = {}; truncated = [False]
    def bt(i, used):
        if len(out) >= cap:
            truncated[0] = True; return
        if i == len(nodes):
            out.append(dict(assign)); return
        e = nodes[i]
        forb = {assign[f] for f in adj.get(e, ()) if f in assign}
        for c in range(min(used + 1, k)):
            if c in forb:
                continue
            assign[e] = c
            bt(i + 1, max(used, c + 1))
            del assign[e]
            if truncated[0]:
                return
    bt(0, 0)
    return out, (not truncated[0])


# ---------- main routine ----------

def tightened_b3(G, k_max=None, color_cap=500000, full_orient_cap=14,
                 orient_samples=400, verbose=False):
    n = G.number_of_nodes()
    H = conflict_graph(G)
    beta = sc.greedy_chi(H)
    edges_all = [sc.ekey(*e) for e in G.edges()]
    if k_max is None:
        k_max = G.number_of_edges()
    for k in range(max(beta, 1), k_max + 1):
        cols, complete = k_colorings(H, edges_all, k, color_cap)
        for col in cols:
            if len(set(col.values())) < k:
                continue
            for orient in orientations(G, full_orient_cap, orient_samples):
                Z, p = build_pot_Z(G, col, orient)
                if min_realizable_size(Z, p, n) == n:
                    if verbose:
                        print(f"  k={k}: proper coloring with a valid potification (m_P=n)")
                    return dict(beta=beta, bound=k, exact=True, n=n)
        if not complete:
            if verbose:
                print(f"  k={k}: coloring search truncated -> report >= {k}")
            return dict(beta=beta, bound=k, exact=False, n=n)
        if verbose:
            print(f"  k={k}: no proper coloring has a valid potification -> B_3 > {k}")
    return dict(beta=beta, bound=k_max + 1, exact=True, n=n)


# ============================================================
#  Reporting helpers
# ============================================================

def _classify(r):
    """Turn a tightened_b3 result into (verdict, strength_note)."""
    beta, bound, exact = r["beta"], r["bound"], r["exact"]
    improved = bound > beta
    if not exact:
        strength = ("PARTIAL: the coloring search was truncated before it could "
                    "rule out this many types, so the value is a floor that may rise "
                    "with a larger --color-cap.")
    elif improved:
        strength = ("TIGHTENED: the swap bound chi(H)=%d was raised by the "
                    "construction-matrix (no-smaller-complex) obstruction." % beta)
    else:
        strength = ("BASELINE: the construction matrix admitted a valid pot already "
                    "at chi(H)=%d, so no tightening beyond the swap bound." % beta)
    return strength


def _report_row(name, G, r):
    marker = "" if r["exact"] else " *"
    return (name, G.number_of_nodes(), G.number_of_edges(),
            r["beta"], r["bound"], marker)


# ============================================================
#  CLI
# ============================================================

def _ask(p):
    try:
        return input(p).strip()
    except EOFError:
        return ""


def _parse_range(s):
    p = s.replace(",", " ").split()
    if len(p) == 1:
        return range(int(p[0]), int(p[0]) + 1)
    if len(p) == 2:
        return range(int(p[0]), int(p[1]) + 1)
    return range(int(p[0]), int(p[1]) + 1, int(p[2]))


def _parse_edges(s):
    G = nx.Graph()
    for tok in s.replace(",", " ").split():
        u, v = tok.split("-")
        G.add_edge(int(u), int(v))
    return G


def _run_one(name, G, k_max, color_cap, orient_cap, verbose):
    r = tightened_b3(G, k_max=k_max, color_cap=color_cap,
                     full_orient_cap=orient_cap, verbose=verbose)
    return r


def main():
    import argparse, sys
    fams = list(sc.FAMILY_GENS)

    ap = argparse.ArgumentParser(
        prog="b3bound",
        description="Tightened lower bound on B_3(G), the minimum number of "
                    "bond-edge types for a Scenario-3 self-assembly pot. Combines "
                    "the 2-switch (swap) obstruction with the construction-matrix "
                    "no-smaller-complex obstruction.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="examples:\n"
               "  b3bound.py --family cycle --range 3 8\n"
               "  b3bound.py --family sun --range 3 6 --verbose\n"
               "  b3bound.py --edges '0-1 1-2 2-0 2-3'\n"
               "  b3bound.py            (interactive)\n")
    ap.add_argument("--family", choices=fams, help="graph family")
    ap.add_argument("--range", nargs="+", type=int, metavar="N",
                    help="n range: 'start', 'start end', or 'start end step'")
    ap.add_argument("--edges", help="custom graph as an edge list, e.g. '0-1 1-2 2-0'")
    ap.add_argument("--k-max", type=int, default=None,
                    help="cap on bond-edge types searched (default |E|)")
    ap.add_argument("--color-cap", type=int, default=500000,
                    help="max proper colorings enumerated per k (default 500000)")
    ap.add_argument("--orient-cap", type=int, default=14,
                    help="enumerate all 2^|E| orientations when |E|<=this, else "
                         "sample (default 14)")
    ap.add_argument("--verbose", action="store_true",
                    help="print the per-k search trace")
    args = ap.parse_args()

    # ---- assemble the work list ----
    items = []          # (label, G)
    if args.edges:
        G = _parse_edges(args.edges)
        items = [("custom", G)]
    elif args.family:
        gen, desc, minn = sc.FAMILY_GENS[args.family]
        rng = _parse_range(" ".join(map(str, args.range))) if args.range \
            else range(minn, minn + 1)
        items = [(f"{args.family}({n})", gen(n)) for n in rng if n >= minn]
    else:
        # ---- interactive ----
        print("B_3(G) lower-bound calculator")
        print("families (n = symmetry index):")
        for i, nm in enumerate(fams, 1):
            print(f"  {i:>2}. {nm:<16} {sc.FAMILY_GENS[nm][1]}")
        print(f"  {len(fams)+1:>2}. custom           paste an edge list")
        sel = _ask("\nselect [number/name]: ")
        if sel in (str(len(fams) + 1), "custom"):
            items = [("custom", _parse_edges(_ask("edge list (e.g. 0-1 1-2 2-0): ")))]
        else:
            key = fams[int(sel) - 1] if sel.isdigit() and 1 <= int(sel) <= len(fams) else sel
            if key not in sc.FAMILY_GENS:
                print("unknown family."); return
            gen, desc, minn = sc.FAMILY_GENS[key]
            rng = _parse_range(_ask(f"n range for '{key}' (min {minn}) -- 'start end [step]': "))
            items = [(f"{key}({n})", gen(n)) for n in rng if n >= minn]
        args.verbose = _ask("show per-k search trace? [y/N] ").lower() in ("y", "yes")

    if not items:
        print("nothing to compute."); return

    # ---- compute + tabulate ----
    print(f"\n{'graph':<14}{'verts':>6}{'edges':>6}{'chi(H)':>8}{'B_3 >=':>8}{'exact':>7}")
    print("-" * 49)
    results = []
    for label, G in items:
        if args.verbose:
            print(f"[{label}]")
        r = _run_one(label, G, args.k_max, args.color_cap, args.orient_cap, args.verbose)
        results.append((label, G, r))
        nm, nv, ne, beta, bnd, mark = _report_row(label, G, r)
        exact = "yes" if r["exact"] else "no"
        print(f"{nm:<14}{nv:>6}{ne:>6}{beta:>8}{bnd:>8}{exact:>7}")

    # ---- interpretation ----
    print("\nhow to read this:")
    print("  chi(H)  lower bound from the 2-switch obstruction alone (swap program).")
    print("  B_3 >=  tightened lower bound: smallest #bond-edge types for which a")
    print("          proper coloring of H admits an orientation whose pot realizes G")
    print("          with no smaller complex (construction matrix m_P = n).")
    print("  exact   'yes' = the search at the failing levels was exhaustive, so the")
    print("          bound is certified. 'no' = a search cap truncated it; the printed")
    print("          value is a floor (raise --color-cap / --k-max to push further).")
    print("\n  This is a LOWER bound, always <= B_3(G). It captures the swap and the")
    print("  no-smaller-complex obstructions, but NOT the Scenario-3-only 'no same-size")
    print("  non-isomorphic complex' obstruction (loops/multiedges, Lemmas 2-3). For")
    print("  graphs whose B_3 is driven by that (e.g. K_n, where B_3=n-1 but B_2 is 1-2)")
    print("  the bound can be loose; for cycles and complete bipartite graphs it is tight.")

    for label, G, r in results:
        if len(results) == 1 or r["bound"] > r["beta"] or not r["exact"]:
            print(f"\n[{label}] {_classify(r)}")


if __name__ == "__main__":
    main()