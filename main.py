import numpy as np
import networkx as nx
from itertools import combinations

# ============================================================
#  Graph families.  Parameter n = the family's symmetry index
#  (e.g. mobius n=3 -> M_6 ; sun n=3 -> 6 vertices).
# ============================================================

def sunlet_graph(n):
    """C_n with one pendant per cycle vertex. 2n verts, 2n edges."""
    G = nx.cycle_graph(n)
    for i in range(n):
        G.add_edge(i, n + i)
    return G

def sun_graph(n):
    """Complete n-sun: clique K_n on 0..n-1; outer vertex n+i adjacent to
    i and (i+1) mod n. 2n verts."""
    G = nx.complete_graph(n)
    for i in range(n):
        G.add_edge(n + i, i)
        G.add_edge(n + i, (i + 1) % n)
    return G

def mobius_ladder_graph(n):
    """Mobius ladder M_2n: cycle C_2n plus n antipodal rungs. n=3 -> M_6. 2n verts."""
    G = nx.cycle_graph(2 * n)
    for i in range(n):
        G.add_edge(i, i + n)
    return G

def web_graph(n):
    """Web graph W_n (Koh/Gallian): circular ladder CL_n with a pendant on each
    outer-cycle vertex. Inner 0..n-1 (cycle), outer n..2n-1 (cycle),
    pendants 2n..3n-1. 3n verts, 4n edges."""
    G = nx.Graph()
    for i in range(n):
        G.add_edge(i, (i + 1) % n)                    # inner cycle
        G.add_edge(n + i, n + (i + 1) % n)            # outer cycle
        G.add_edge(i, n + i)                          # radial rung
        G.add_edge(n + i, 2 * n + i)                  # pendant off outer cycle
    return G

# name -> (generator, vertex-count note, min_n)
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
}

# ============================================================
#  Core: alternating 4-cycles, the triangle criterion, and the
#  exhaustive isomorphism proof.
# ============================================================

def alternating_4cycles(G):
    """Yield (v,w,s,t, delta) for every valid 2-switch target:
    edges vw,st present; vt,sw absent; 4 distinct vertices.
    d_ij = (A^2)[i,j]; delta = (d_vt+d_sw)-(d_vw+d_st)-2*(A[v,s]+A[w,t])."""
    nodes = list(G.nodes())
    idx = {v: k for k, v in enumerate(nodes)}
    A = nx.to_numpy_array(G, nodelist=nodes, dtype=np.int64)
    D = A @ A
    edges = [(idx[u], idx[v]) for u, v in G.edges()]
    for (v, w), (a, b) in combinations(edges, 2):
        for s, t in ((a, b), (b, a)):
            if len({v, w, s, t}) < 4:
                continue
            if A[v, t] or A[s, w]:
                continue
            d = int((D[v, t] + D[s, w]) - (D[v, w] + D[s, t]) - 2 * (A[v, s] + A[w, t]))
            yield nodes[v], nodes[w], nodes[s], nodes[t], d

def switch(G, v, w, s, t):
    """Apply the 2-switch: drop vw,st ; add vt,sw."""
    H = G.copy()
    H.remove_edges_from([(v, w), (s, t)])
    H.add_edges_from([(v, t), (s, w)])
    return H

def certify_triangle(G):
    """Sufficient test: True iff every switch target has delta != 0
    (=> G provably unswappable by the triangle invariant)."""
    return all(d != 0 for *_, d in alternating_4cycles(G))

def verify(G):
    """Ground truth via full isomorphism checks on every switch target.
       total  : number of switch targets
       nz, z  : count with delta != 0 and delta == 0
       ziso   : delta == 0 AND actually isomorphic (real swaps the Δ-test misses)
       viol   : first target with delta != 0 yet isomorphic (would disprove lemma)
       unswap : True iff NO switch on G yields an isomorphic graph."""
    total = nz = z = ziso = 0
    viol = None
    any_iso = False
    for v, w, s, t, d in alternating_4cycles(G):
        total += 1
        iso = nx.is_isomorphic(G, switch(G, v, w, s, t))
        any_iso |= iso
        if d == 0:
            z += 1
            ziso += iso
        else:
            nz += 1
            if iso and viol is None:
                viol = (v, w, s, t)
    return dict(total=total, nz=nz, z=z, ziso=ziso, viol=viol, unswap=not any_iso)

# ============================================================
#  CLI: pick a family and an n-range, run the exhaustive proof.
# ============================================================

def ask(prompt):
    try:
        return input(prompt).strip()
    except EOFError:
        return ""

def parse_range(s):
    p = s.replace(",", " ").split()
    if len(p) == 1:
        return range(int(p[0]), int(p[0]) + 1)
    if len(p) == 2:
        return range(int(p[0]), int(p[1]) + 1)
    return range(int(p[0]), int(p[1]) + 1, int(p[2]))

def main():
    names = list(FAMILY_GENS)
    print("families (n = symmetry index):")
    for i, nm in enumerate(names, 1):
        print(f"  {i:>2}. {nm:<16} {FAMILY_GENS[nm][1]}")

    sel = ask("\nselect family [number or name]: ")
    key = names[int(sel) - 1] if sel.isdigit() and 1 <= int(sel) <= len(names) else sel
    if key not in FAMILY_GENS:
        print("unknown family.")
        return
    gen, desc, minn = FAMILY_GENS[key]

    rng = parse_range(ask(f"n range for '{key}' (min {minn}) -- 'start end [step]': "))

    print(f"\n{key}: {desc}    (exhaustive isomorphism proof)")
    print(f"{'n':>4}{'verts':>7}{'edges':>7}{'targets':>9}{'Δ-cert':>8}{'verdict':>14}")
    print("-" * 57)
    for n in rng:
        if n < minn:
            print(f"{n:>4}   skipped (< min {minn})")
            continue
        G = gen(n)
        r = verify(G)
        cert = "YES" if certify_triangle(G) else "no"
        verdict = "unswappable" if r["unswap"] else "∃ swap"
        print(f"{n:>4}{G.number_of_nodes():>7}{G.number_of_edges():>7}"
              f"{r['total']:>9}{cert:>8}{verdict:>14}")

if __name__ == "__main__":
    main()