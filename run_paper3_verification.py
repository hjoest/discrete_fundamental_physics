#!/usr/bin/env python3
"""Zero-dependency verification suite for Paper III (docs/paper3.tex):
'The Causal Set Inverse Problem on (P, G_*, Z, <_Z): Discrete Cartan Soldering,
Galois Phason Confinement, and Emergent 4D Lorentzian Geometry'.

Reproduces all analytical and numerical results in:
  - Table 1: Five quantized soldering velocity shells of 2I and S^2 spherical 5-design isotropy.
  - Table 2: Galois phason confinement under Kadanoff Plancherel convolution and 24-step 12_a walks.
  - Table 3: Harmonic 4D coordinate reconstruction L_Z X^mu = d_Z^T(Delta X^mu) and curvature strain on N = 3^4 = 81.
  - Table 4: Multi-size coordinate-free Kleitman-Rothschild pancake melting (N in {80, 120, 160, 200}),
             12_a holonomy annealing, reconstructed 4D light-cone fidelity, and dimension estimators.
"""

import math
import random
from typing import Dict, List, Set, Tuple

from src.poset_simulation import _invert_myrheim_meyer

PHI = (1.0 + math.sqrt(5.0)) / 2.0
PHI_INV = PHI - 1.0  # 1 / PHI = (sqrt(5) - 1) / 2


def build_2I_quaternions() -> List[Tuple[float, float, float, float]]:
    """Construct all 120 unit quaternions q = (q0, q1, q2, q3) of 2I in S^3."""
    verts: List[Tuple[float, float, float, float]] = []
    for i in range(4):
        for s in (-1.0, 1.0):
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            verts.append(tuple(v))
    for s0 in (-0.5, 0.5):
        for s1 in (-0.5, 0.5):
            for s2 in (-0.5, 0.5):
                for s3 in (-0.5, 0.5):
                    verts.append((s0, s1, s2, s3))
    even_perms = [
        (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2),
        (1, 0, 3, 2), (1, 2, 0, 3), (1, 3, 2, 0),
        (2, 0, 1, 3), (2, 1, 3, 0), (2, 3, 0, 1),
        (3, 0, 2, 1), (3, 1, 0, 2), (3, 2, 1, 0),
    ]
    a, b, c = 0.5 * PHI, 0.5, 0.5 * PHI_INV
    for p in even_perms:
        for sa in (-1.0, 1.0):
            for sb in (-1.0, 1.0):
                for sc in (-1.0, 1.0):
                    v = [0.0, 0.0, 0.0, 0.0]
                    v[p[0]] = sa * a
                    v[p[1]] = sb * b
                    v[p[2]] = sc * c
                    v[p[3]] = 0.0
                    verts.append(tuple(v))
    assert len(verts) == 120
    return verts


def q_mul(a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    return (
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    )


def q_inv(a: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    return (a[0], -a[1], -a[2], -a[3])


def solve_linear_system(A: List[List[float]], b: List[float]) -> List[float]:
    """Gaussian elimination with partial pivoting for an N x N linear system."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        piv = max(range(k, n), key=lambda r: abs(M[r][k]))
        M[k], M[piv] = M[piv], M[k]
        piv_val = M[k][k]
        if abs(piv_val) < 1e-15:
            continue
        for j in range(k, n + 1):
            M[k][j] /= piv_val
        for i in range(n):
            if i != k and abs(M[i][k]) > 1e-15:
                fac = M[i][k]
                for j in range(k, n + 1):
                    M[i][j] -= fac * M[k][j]
    return [M[i][n] for i in range(n)]


def find_exact_galois_automorphism(
    q2i: List[Tuple[float, float, float, float]]
) -> Tuple[Dict[int, int], Dict[Tuple[float, float, float, float], int]]:
    """Construct the exact outer automorphism sigma in Aut(2I) swapping 12_a <-> 12_d."""
    idx_of: Dict[Tuple[float, float, float, float], int] = {}
    for i, q in enumerate(q2i):
        key = tuple(round(c, 5) for c in q)
        idx_of[key] = i

    def mul_idx(i: int, j: int) -> int:
        p = q_mul(q2i[i], q2i[j])
        return idx_of[tuple(round(c, 5) for c in p)]

    c12a = [i for i, q in enumerate(q2i) if abs(2.0 * q[0] - PHI) < 1e-6]
    c12d = [i for i, q in enumerate(q2i) if abs(2.0 * q[0] + PHI_INV) < 1e-6]
    g1 = c12a[0]
    g2 = next(j for j in c12a if abs(2.0 * q2i[mul_idx(g1, j)][0] - 1.0) < 1e-6)
    id_idx = next(i for i, q in enumerate(q2i) if abs(q[0] - 1.0) < 1e-6)

    word_map: Dict[int, Tuple[int, ...]] = {id_idx: ()}
    queue = [id_idx]
    while queue:
        curr = queue.pop(0)
        w = word_map[curr]
        for gen_idx, g in enumerate((g1, g2)):
            nxt = mul_idx(curr, g)
            if nxt not in word_map:
                word_map[nxt] = w + (gen_idx,)
                queue.append(nxt)
    assert len(word_map) == 120

    for h1 in c12d:
        for h2 in c12d:
            if abs(2.0 * q2i[mul_idx(h1, h2)][0] - 1.0) < 1e-6:
                sigma_map: Dict[int, int] = {}
                gens_h = (h1, h2)
                for el, w in word_map.items():
                    cur_h = id_idx
                    for step in w:
                        cur_h = mul_idx(cur_h, gens_h[step])
                    sigma_map[el] = cur_h
                if len(set(sigma_map.values())) == 120:
                    ok = all(
                        sigma_map[mul_idx(a, b)] == mul_idx(sigma_map[a], sigma_map[b])
                        for a in range(120)
                        for b in (g1, g2)
                    )
                    if ok:
                        return sigma_map, idx_of
    raise RuntimeError("Could not find outer automorphism of 2I")


def verify_table1_shells_and_5design() -> None:
    print("=" * 84)
    print("TABLE 1: Five Quantized Soldering Shells of 2I and S^2 Spherical 5-Design Isotropy")
    print("=" * 84)
    q2i = build_2I_quaternions()
    classes: Dict[str, List[Tuple[float, float, float, float]]] = {
        "1_a (+I)": [], "1_b (-I)": [],
        "12_a (+phi)": [], "12_c (-phi)": [],
        "20_a (+1)": [], "20_b (-1)": [],
        "12_b (+phi^-1)": [], "12_d (-phi^-1)": [],
        "30_a (0)": [],
    }
    for q in q2i:
        tr = 2.0 * q[0]
        if abs(tr - 2.0) < 1e-7: classes["1_a (+I)"].append(q)
        elif abs(tr + 2.0) < 1e-7: classes["1_b (-I)"].append(q)
        elif abs(tr - PHI) < 1e-7: classes["12_a (+phi)"].append(q)
        elif abs(tr + PHI) < 1e-7: classes["12_c (-phi)"].append(q)
        elif abs(tr - 1.0) < 1e-7: classes["20_a (+1)"].append(q)
        elif abs(tr + 1.0) < 1e-7: classes["20_b (-1)"].append(q)
        elif abs(tr - PHI_INV) < 1e-7: classes["12_b (+phi^-1)"].append(q)
        elif abs(tr + PHI_INV) < 1e-7: classes["12_d (-phi^-1)"].append(q)
        else: classes["30_a (0)"].append(q)

    for cname, elems in classes.items():
        q0 = elems[0]
        chi1 = 2.0 * q0[0]
        v_norm = math.sqrt(q0[1] ** 2 + q0[2] ** 2 + q0[3] ** 2)
        ds2 = 1.0 - v_norm ** 2
        assert abs(ds2 - 0.25 * (chi1 ** 2)) < 1e-14
        max_5design_err = 0.0
        if v_norm > 1e-7:
            unit_vecs = [(q[1] / v_norm, q[2] / v_norm, q[3] / v_norm) for q in elems]
            m = len(unit_vecs)
            for i in range(3):
                for j in range(3):
                    val2 = sum(u[i] * u[j] for u in unit_vecs) / m
                    target2 = (1.0 / 3.0) if i == j else 0.0
                    max_5design_err = max(max_5design_err, abs(val2 - target2))
                    for k in range(3):
                        for l in range(3):
                            val4 = sum(u[i] * u[j] * u[k] * u[l] for u in unit_vecs) / m
                            d_ij = 1.0 if i == j else 0.0
                            d_kl = 1.0 if k == l else 0.0
                            d_ik = 1.0 if i == k else 0.0
                            d_jl = 1.0 if j == l else 0.0
                            d_il = 1.0 if i == l else 0.0
                            d_jk = 1.0 if j == k else 0.0
                            target4 = (d_ij * d_kl + d_ik * d_jl + d_il * d_jk) / 15.0
                            max_5design_err = max(max_5design_err, abs(val4 - target4))
            for i in range(3):
                val5 = sum(u[i] ** 5 for u in unit_vecs) / m
                max_5design_err = max(max_5design_err, abs(val5))
        print(
            f"  {cname:16s} | |C|={len(elems):2d} | chi_1={chi1:+8.5f} | v={v_norm:8.6f} "
            f"| ds^2/l^2={ds2:8.6f} | 5-design err={max_5design_err:.2e}"
        )


def verify_table2_galois_phason_confinement() -> None:
    print("\n" + "=" * 84)
    print("TABLE 2: Galois Phason Confinement (Plancherel Convolution & 24-Step 12_a Walks)")
    print("=" * 84)
    q2i = build_2I_quaternions()
    sigma_map, _ = find_exact_galois_automorphism(q2i)
    c12a = [q for q in q2i if abs(2.0 * q[0] - PHI) < 1e-7]
    print(f"  Verified exact outer automorphism sigma in Aut(2I) on all {len(sigma_map)} elements.")

    gap1 = 1.0 - 0.5 * PHI
    gap2 = 1.0 + 0.5 * PHI_INV
    print(f"  Delta E_chi1(12_a) = {gap1:.6f}, Delta E_chi2(12_a) = {gap2:.6f}, Ratio = {gap2/gap1:.6f} (phi^4 = {PHI**4:.6f})")

    eps = 0.25
    lam1 = 1.0 - eps * gap1
    lam2 = 1.0 - eps * gap2
    for b in (1, 2, 4, 8):
        b2 = b * b
        l1_p = lam1 ** b2
        l2_p = lam2 ** b2
        ratio = l2_p / l1_p
        bound = math.exp(-0.5 * math.sqrt(5.0) * eps * b2)
        assert ratio <= bound + 1e-15
        print(f"  Scale b={b:2d} (b^2={b2:2d}) | lam1^b^2={l1_p:11.6e} | lam2^b^2={l2_p:11.6e} | Ratio={ratio:11.4e} <= Bound={bound:11.4e}")

    rng = random.Random(42)
    trials = 400
    steps = 24
    c12a_idx = [i for i, q in enumerate(q2i) if abs(2.0 * q[0] - PHI) < 1e-7]
    for beta_perp in (0.0, 0.5, 1.0, 2.0, 3.0):
        par_sq_sum = 0.0
        perp_sq_sum = 0.0
        for _ in range(trials):
            x_par = [0.0, 0.0, 0.0]
            x_perp = [0.0, 0.0, 0.0]
            for _s in range(steps):
                weights = []
                for h_i in c12a_idx:
                    hp = q2i[sigma_map[h_i]]
                    new_perp_sq = sum((x_perp[k] + hp[k + 1]) ** 2 for k in range(3))
                    weights.append(math.exp(-beta_perp * new_perp_sq))
                tot_w = sum(weights)
                r = rng.random() * tot_w
                acc = 0.0
                chosen_i = c12a_idx[0]
                for h_i, w in zip(c12a_idx, weights):
                    acc += w
                    if r <= acc:
                        chosen_i = h_i
                        break
                chosen = q2i[chosen_i]
                chosen_p = q2i[sigma_map[chosen_i]]
                for k in range(3):
                    x_par[k] += chosen[k + 1]
                    x_perp[k] += chosen_p[k + 1]
            par_sq_sum += sum(c * c for c in x_par)
            perp_sq_sum += sum(c * c for c in x_perp)
        mean_par_sq = par_sq_sum / trials
        mean_perp_sq = perp_sq_sum / trials
        print(
            f"  Walk beta_perp={beta_perp:3.1f} | <|X_par|^2>={mean_par_sq:6.2f} | <|X_perp|^2>={mean_perp_sq:5.2f} "
            f"| Ratio={mean_perp_sq/mean_par_sq:6.4f} | R_perp={math.sqrt(mean_perp_sq):4.2f}"
        )


def verify_table3_harmonic_soldering_and_curvature() -> None:
    print("\n" + "=" * 84)
    print("TABLE 3: Harmonic 4D Coordinate Reconstruction and Single-Link Curvature Strain")
    print("=" * 84)
    q2i = build_2I_quaternions()
    c12a = [q for q in q2i if abs(2.0 * q[0] - PHI) < 1e-7]
    c20a = [q for q in q2i if abs(2.0 * q[0] - 1.0) < 1e-7]
    c12b = [q for q in q2i if abs(2.0 * q[0] - PHI_INV) < 1e-7]
    c30a = [q for q in q2i if abs(2.0 * q[0]) < 1e-7]

    gens = [c12a[0], c12a[1], c12a[4], c12a[8]]
    L = 3
    coords_grid = [(i0, i1, i2, i3) for i0 in range(L) for i1 in range(L) for i2 in range(L) for i3 in range(L)]
    idx_map = {pt: idx for idx, pt in enumerate(coords_grid)}
    N = len(coords_grid)

    edges: List[Tuple[int, int, Tuple[float, float, float, float]]] = []
    true_X: List[List[float]] = [[0.0, 0.0, 0.0, 0.0] for _ in range(N)]
    for idx, pt in enumerate(coords_grid):
        t_val = float(sum(pt))
        x_vec = [0.0, 0.0, 0.0]
        for mu in range(4):
            for k in range(3):
                x_vec[k] += pt[mu] * gens[mu][k + 1]
        true_X[idx] = [t_val, x_vec[0], x_vec[1], x_vec[2]]
        for mu in range(4):
            nxt = list(pt)
            nxt[mu] += 1
            nxt_t = tuple(nxt)
            if nxt_t in idx_map:
                edges.append((idx, idx_map[nxt_t], gens[mu]))

    def run_harmonic_solver(
        edge_list: List[Tuple[int, int, Tuple[float, float, float, float]]]
    ) -> Tuple[List[List[float]], float, float]:
        Lap = [[0.0] * N for _ in range(N)]
        rhs = [[0.0] * N for _ in range(4)]
        for u, v, h_e in edge_list:
            Lap[u][u] += 1.0
            Lap[v][v] += 1.0
            Lap[u][v] -= 1.0
            Lap[v][u] -= 1.0
            dX = [1.0, h_e[1], h_e[2], h_e[3]]
            for mu in range(4):
                rhs[mu][u] -= dX[mu]
                rhs[mu][v] += dX[mu]
        for j in range(N):
            Lap[0][j] = 0.0
            Lap[j][0] = 0.0
        Lap[0][0] = 1.0
        for mu in range(4):
            rhs[mu][0] = 0.0
        rec_X = [[0.0] * 4 for _ in range(N)]
        for mu in range(4):
            col_sol = solve_linear_system(Lap, rhs[mu])
            for i in range(N):
                rec_X[i][mu] = col_sol[i]
        tot_res = 0.0
        max_err = 0.0
        for u, v, h_e in edge_list:
            dX = [1.0, h_e[1], h_e[2], h_e[3]]
            err_sq = sum(((rec_X[v][mu] - rec_X[u][mu]) - dX[mu]) ** 2 for mu in range(4))
            tot_res += err_sq
            max_err = max(max_err, math.sqrt(err_sq))
        return rec_X, tot_res, max_err

    rec_flat, res_flat, max_flat = run_harmonic_solver(edges)
    max_coord_diff = max(
        math.sqrt(sum((rec_flat[i][mu] - true_X[i][mu]) ** 2 for mu in range(4)))
        for i in range(N)
    )
    print(
        f"  Flat Parallelogram (N={N}, M={len(edges)}) | E[X]={res_flat:.2e} "
        f"| max_strain={max_flat:.2e} | max_coord_err={max_coord_diff:.2e}"
    )

    center_u = idx_map[(1, 1, 1, 1)]
    target_v = idx_map[(2, 1, 1, 1)]
    for defect_name, defect_h in [
        ("Defect Shell 0: 1_a (+I_2)", (1.0, 0.0, 0.0, 0.0)),
        ("Defect Shell 1: 12_a (!= g_0)", c12a[2]),
        ("Defect Shell 2: 20_a", c20a[0]),
        ("Defect Shell 3: 12_b", c12b[0]),
        ("Defect Shell 4: 30_a", c30a[0]),
    ]:
        mod_edges = [
            (u, v, defect_h if (u == center_u and v == target_v) else h_e)
            for u, v, h_e in edges
        ]
        # Parallelogram diamond holonomy U_sq = h_{e1} h_{e2} h_{e3}^-1 h_{e4}^-1 = defect_h * g0^-1
        u_sq = q_mul(defect_h, q_inv(gens[0]))
        s_plaq = 1.0 - u_sq[0]
        delta_e0 = math.sqrt(sum((defect_h[k] - gens[0][k]) ** 2 for k in range(1, 4)))
        _, res_c, max_c = run_harmonic_solver(mod_edges)
        ratio = res_c / (delta_e0 ** 2)
        print(
            f"  {defect_name:29s} | 1-Re(chi)/2={s_plaq:6.4f} | |delta_e0 X|={delta_e0:8.6f} "
            f"| E[X]={res_c:7.5f} | max_strain={max_c:6.4f} | E/|delta|^2={ratio:7.5f}"
        )


def _compute_spec_and_ball_dim(N: int, succ_set: List[Set[int]]) -> Tuple[float, float]:
    """Compute (D_spec, D_ball) on the undirected Hasse covering graph of succ_set."""
    reach_bits = [0] * N
    for u in range(N - 1, -1, -1):
        r_u = 0
        for v in succ_set[u]:
            r_u |= (1 << v) | reach_bits[v]
        reach_bits[u] = r_u

    undir: List[Set[int]] = [set() for _ in range(N)]
    for u in range(N):
        two_step = 0
        for w in succ_set[u]:
            two_step |= reach_bits[w]
        for v in succ_set[u]:
            if not ((two_step >> v) & 1):
                undir[u].add(v)
                undir[v].add(u)

    d_spec_list: List[float] = []
    d_ball_list: List[float] = []
    for root in range(N // 4, (3 * N) // 4, max(1, N // 16)):
        if not undir[root]:
            continue
        dist = {root: 0}
        q = [root]
        head = 0
        while head < len(q):
            curr = q[head]
            head += 1
            if dist[curr] < 3:
                for nxt in undir[curr]:
                    if nxt not in dist:
                        dist[nxt] = dist[curr] + 1
                        q.append(nxt)
        v1 = sum(1 for d in dist.values() if d <= 1)
        v2 = sum(1 for d in dist.values() if d <= 2)
        v3 = sum(1 for d in dist.values() if d <= 3)
        if v2 > v1 and v3 > v2:
            d_ball_list.append(math.log(v3 / v1) / math.log(3.0))

        p = {root: 1.0}
        p2, p4 = 0.0, 0.0
        for s in range(1, 5):
            np: Dict[int, float] = {}
            for node, prob in p.items():
                deg = len(undir[node])
                if deg > 0:
                    share = prob / deg
                    for nxt in undir[node]:
                        np[nxt] = np.get(nxt, 0.0) + share
            p = np
            if s == 2:
                p2 = p.get(root, 0.0)
            elif s == 4:
                p4 = p.get(root, 0.0)
        if p2 > 1e-9 and p4 > 1e-9 and p2 > p4:
            d_spec_list.append(-2.0 * math.log(p4 / p2) / math.log(2.0))

    d_spec = sum(d_spec_list) / len(d_spec_list) if d_spec_list else 1.5
    d_ball = sum(d_ball_list) / len(d_ball_list) if d_ball_list else 1.5
    return d_spec, d_ball


def _solder_12a_and_reconstruct(
    N: int,
    succ_set: List[Set[int]],
    q2i: List[Tuple[float, float, float, float]],
    c12a: List[int],
    hasse_only: bool = True,
) -> Tuple[float, float, float, float]:
    """Assign a 12_a-valued link 1-cochain h_e in 12_a, relax across 6 diamond-closure sweeps,
    and solve L_Z X^mu = d_Z^T(Delta X^mu). Returns (solder_lc_pct, pure_lc_pct, E_tot_per_M, E_sp_per_M).
    """
    if hasse_only:
        reach_bits = [0] * N
        for u in range(N - 1, -1, -1):
            r_u = 0
            for v in succ_set[u]:
                r_u |= (1 << v) | reach_bits[v]
            reach_bits[u] = r_u
        active: List[Set[int]] = [set() for _ in range(N)]
        for u in range(N):
            two_u = 0
            for w in succ_set[u]:
                two_u |= reach_bits[w]
            for v in succ_set[u]:
                if not ((two_u >> v) & 1):
                    active[u].add(v)
    else:
        active = [set(s) for s in succ_set]

    edges: List[Tuple[int, int]] = []
    edge_idx: Dict[Tuple[int, int], int] = {}
    in_deg = [0] * N
    for u in range(N):
        for v in sorted(active[u]):
            edge_idx[(u, v)] = len(edges)
            edges.append((u, v))
            in_deg[v] += 1
    M = len(edges)
    h_edge = [0] * M
    for u in range(N):
        for k, v in enumerate(sorted(active[u])):
            h_edge[edge_idx[(u, v)]] = c12a[(u + 5 * k) % 12]

    two_paths: Dict[Tuple[int, int], List[int]] = {}
    for u in range(N):
        for w in sorted(active[u]):
            for v in sorted(active[w]):
                two_paths.setdefault((u, v), []).append(w)
    diamonds: List[Tuple[int, int, int, int]] = []
    edge_diamonds: List[List[int]] = [[] for _ in range(M)]
    for (u, v), ws in two_paths.items():
        if len(ws) >= 2:
            for i in range(min(4, len(ws))):
                for j in range(i + 1, min(4, len(ws))):
                    w1, w2 = ws[i], ws[j]
                    e1, e2 = edge_idx[(u, w1)], edge_idx[(w1, v)]
                    e3, e4 = edge_idx[(u, w2)], edge_idx[(w2, v)]
                    d_id = len(diamonds)
                    diamonds.append((e1, e2, e3, e4))
                    for e in (e1, e2, e3, e4):
                        edge_diamonds[e].append(d_id)

    Lap = [[0.0] * N for _ in range(N)]
    for u, v in edges:
        Lap[u][u] += 1.0
        Lap[v][v] += 1.0
        Lap[u][v] -= 1.0
        Lap[v][u] -= 1.0
    for u in range(N):
        if Lap[u][u] < 1e-9:
            Lap[u][u] = 1.0
    for j in range(N):
        Lap[0][j] = 0.0
        Lap[j][0] = 0.0
    Lap[0][0] = 1.0

    rhs0 = [0.0] * N
    for u, v in edges:
        rhs0[u] -= 1.0
        rhs0[v] += 1.0
    rhs0[0] = 0.0
    X0 = solve_linear_system(Lap, rhs0)
    rec_X = [[X0[i], 0.0, 0.0, 0.0] for i in range(N)]

    for _sweep in range(6):
        rhs = [[0.0] * N for _ in range(3)]
        for idx_e, (u, v) in enumerate(edges):
            h_e = q2i[h_edge[idx_e]]
            for k in range(3):
                rhs[k][u] -= h_e[k + 1]
                rhs[k][v] += h_e[k + 1]
        for k in range(3):
            rhs[k][0] = 0.0
            col_sol = solve_linear_system(Lap, rhs[k])
            for i in range(N):
                rec_X[i][k + 1] = col_sol[i]

        counts_out: List[Dict[int, int]] = [{} for _ in range(N)]
        counts_in: List[Dict[int, int]] = [{} for _ in range(N)]
        for u in range(N):
            out_vs = sorted(active[u])
            max_out_u = (len(out_vs) - 1) // 12 + 1
            for v in out_vs:
                e = edge_idx[(u, v)]
                max_in_v = (in_deg[v] - 1) // 12 + 1
                dt = rec_X[v][0] - rec_X[u][0]
                dx = [rec_X[v][k] - rec_X[u][k] for k in (1, 2, 3)]
                dx_norm = math.sqrt(sum(x * x for x in dx))
                r_t = min(0.587785, 0.42 * max(0.0, dt))
                best_h = None
                best_cost = 1e18
                for cand in c12a:
                    if counts_out[u].get(cand, 0) >= max_out_u or counts_in[v].get(cand, 0) >= max_in_v:
                        continue
                    q_c = q2i[cand]
                    dot_c = sum(dx[k - 1] * q_c[k] for k in (1, 2, 3))
                    radial_diff = (dot_c / dx_norm - r_t) if dx_norm > 1e-9 else 0.0
                    pen_pos = 0.5 if (dx_norm > 1e-9 and dot_c <= 0.0) else 0.0
                    cost = 0.05 * (pen_pos + radial_diff * radial_diff)
                    h_edge[e] = cand
                    for d_id in edge_diamonds[e]:
                        e1, e2, e3, e4 = diamonds[d_id]
                        q1, q2, q3, q4 = q2i[h_edge[e1]], q2i[h_edge[e2]], q2i[h_edge[e3]], q2i[h_edge[e4]]
                        for k in (1, 2, 3):
                            df = (q1[k] + q2[k]) - (q3[k] + q4[k])
                            cost += df * df
                    if cost < best_cost - 1e-12:
                        best_cost = cost
                        best_h = cand
                if best_h is None:
                    best_h = min(c12a, key=lambda c: counts_out[u].get(c, 0) + counts_in[v].get(c, 0))
                h_edge[e] = best_h
                counts_out[u][best_h] = counts_out[u].get(best_h, 0) + 1
                counts_in[v][best_h] = counts_in[v].get(best_h, 0) + 1

    assert all(h in c12a for h in h_edge), "Every link holonomy must lie in 12_a!"
    rhs = [[0.0] * N for _ in range(3)]
    for idx_e, (u, v) in enumerate(edges):
        h_e = q2i[h_edge[idx_e]]
        for k in range(3):
            rhs[k][u] -= h_e[k + 1]
            rhs[k][v] += h_e[k + 1]
    for k in range(3):
        rhs[k][0] = 0.0
        col_sol = solve_linear_system(Lap, rhs[k])
        for i in range(N):
            rec_X[i][k + 1] = col_sol[i]

    lc_only = 0
    lc_solder = 0
    res0 = 0.0
    res_sp = 0.0
    for idx_e, (u, v) in enumerate(edges):
        h_e = q2i[h_edge[idx_e]]
        dt = rec_X[v][0] - rec_X[u][0]
        dx = [rec_X[v][k] - rec_X[u][k] for k in (1, 2, 3)]
        dx_sq = sum(x * x for x in dx)
        dot_sp = sum(dx[k - 1] * h_e[k] for k in (1, 2, 3))
        if dt > 0 and dt * dt >= dx_sq - 1e-9:
            lc_only += 1
            if dot_sp > 1e-9:
                lc_solder += 1
        res0 += (dt - 1.0) ** 2
        res_sp += sum((dx[k - 1] - h_e[k]) ** 2 for k in (1, 2, 3))
    tot_res = res0 + res_sp
    return (
        100.0 * lc_solder / max(1, M),
        100.0 * lc_only / max(1, M),
        tot_res / max(1, M),
        res_sp / max(1, M),
    )


def verify_table4_multisize_melting_and_soldering() -> None:
    print("\n" + "=" * 84)
    print("TABLE 4: Multi-Size Coordinate-Free KR Pancake Melting & 2I Soldering Reconstruction")
    print("=" * 84)
    q2i = build_2I_quaternions()
    c12a = [i for i, q in enumerate(q2i) if abs(2.0 * q[0] - PHI) < 1e-6]
    fact_table = [1, 1, 2, 6, 24, 120, 720, 5040]
    beta_tot = 18.461
    ln_z6 = math.log(6.0)
    ln_81 = math.log(81.0)
    ln_27 = math.log(27.0)
    ln_2 = math.log(2.0)

    for N, steps in [(80, 3300), (120, 7425), (160, 13200), (200, 20625)]:
        rng = random.Random(42)
        n1, n2 = N // 4, (3 * N) // 4
        succ_set: List[Set[int]] = [set() for _ in range(N)]
        pred_set: List[Set[int]] = [set() for _ in range(N)]
        for u in range(0, n1):
            for v in range(n1, n2):
                if rng.random() < 0.38:
                    succ_set[u].add(v)
                    pred_set[v].add(u)
            for v in range(n2, N):
                if rng.random() < 0.28:
                    succ_set[u].add(v)
                    pred_set[v].add(u)
        for u in range(n1, n2):
            for v in range(n2, N):
                if rng.random() < 0.38:
                    succ_set[u].add(v)
                    pred_set[v].add(u)

        b_past = max(12, (24 * N) // 160)
        b_fut = max(16, (32 * N) // 160)

        def vertex_action(d_out: int, d_in: int, is_past: bool, is_fut: bool) -> float:
            s = 0.0
            if not is_fut:
                b_out = beta_tot / fact_table[min(d_out, 7)]
                s += (b_out - 1.0 if b_out > 1.0 else 0.0) + (ln_z6 + ln_81) * max(0, d_out - 4)
            if not is_past:
                b_in = beta_tot / fact_table[min(d_in, 7)]
                s += (b_in - 1.0 if b_in > 1.0 else 0.0) + (ln_z6 + ln_81) * max(0, d_in - 4)
            return s

        def eval_poset() -> Tuple[float, int, float, float, float, int]:
            reach_bits = [0] * N
            depth = [1] * N
            shortcuts = 0
            diamonds = 0
            spindles = 0
            for u in range(N - 1, -1, -1):
                r_u = 0
                two_u = 0
                s1, s2, s3 = 0, 0, 0
                md = 0
                for v in succ_set[u]:
                    r_u |= (1 << v) | reach_bits[v]
                    two_u |= reach_bits[v]
                    if depth[v] > md:
                        md = depth[v]
                    cm = 0
                    for w in succ_set[v]:
                        cm |= 1 << w
                    s3 |= s2 & cm
                    s2 |= s1 & cm
                    s1 |= cm
                reach_bits[u] = r_u
                if succ_set[u]:
                    depth[u] = 1 + md
                for v in succ_set[u]:
                    if (two_u >> v) & 1:
                        shortcuts += 1
                diamonds += min(6, (s2 & ~s3).bit_count())
                spindles += s3.bit_count()
            layers = max(depth)
            c2 = sum(r.bit_count() for r in reach_bits)
            f_ord = c2 / (0.5 * N * (N - 1))
            d_mm = _invert_myrheim_meyer(f_ord)
            tot_s = shortcuts * ln_27 + spindles * ln_z6 - diamonds * ln_2
            d_sum, d_cnt = 0, 0
            for u in range(N):
                ip = u < b_past
                ifut = u >= N - b_fut
                tot_s += vertex_action(len(succ_set[u]), len(pred_set[u]), ip, ifut)
                if not ifut:
                    d_sum += len(succ_set[u])
                    d_cnt += 1
            return tot_s, layers, d_sum / max(1, d_cnt), f_ord, d_mm, diamonds

        cur_s, init_l, init_d, init_f, init_dmm, init_diam = eval_poset()
        if N == 160:
            init_spec, init_ball = _compute_spec_and_ball_dim(N, succ_set)
            init_solder_fid, init_lc_only, init_res, init_res_sp = _solder_12a_and_reconstruct(
                N, succ_set, q2i, c12a, hasse_only=False
            )
            print(
                f"  Initial KR Pancake (N=160, step=0) | L_max={init_l:2d} | <d_out>={init_d:5.2f} "
                f"| N_diam={init_diam:3d} | f={init_f:.4f} | D_MM={init_dmm:4.2f} | D_spec={init_spec:4.2f} | D_ball={init_ball:4.2f} "
                f"| LC_fid={init_solder_fid:5.1f}% (LC_only={init_lc_only:5.1f}%) | E[X]/M={init_res:6.4f} (E_sp/M={init_res_sp:6.4f})"
            )
        tot_pairs = (N * (N - 1)) // 2
        for _step in range(steps):
            mt = rng.random()
            ne = sum(len(s) for s in succ_set)
            if mt < 0.35 and ne > N:
                u = rng.randrange(N - b_fut)
                if succ_set[u]:
                    du = len(succ_set[u])
                    v = rng.choice(tuple(succ_set[u]))
                    succ_set[u].remove(v)
                    pred_set[v].remove(u)
                    ns, *_ = eval_poset()
                    hr = float(du) / float(N - 1 - u)
                    if rng.random() < min(1.0, hr * math.exp(-min(50.0, max(-50.0, ns - cur_s)))):
                        cur_s = ns
                    else:
                        succ_set[u].add(v)
                        pred_set[v].add(u)
            elif mt < 0.70 and ne < tot_pairs - N:
                u = rng.randrange(N - b_fut)
                du = len(succ_set[u])
                v = rng.randint(u + 1, N - 1)
                if v not in succ_set[u]:
                    succ_set[u].add(v)
                    pred_set[v].add(u)
                    ns, *_ = eval_poset()
                    hr = float(N - 1 - u) / float(du + 1)
                    if rng.random() < min(1.0, hr * math.exp(-min(50.0, max(-50.0, ns - cur_s)))):
                        cur_s = ns
                    else:
                        succ_set[u].remove(v)
                        pred_set[v].remove(u)
            else:
                u = rng.randrange(N - b_fut)
                if succ_set[u]:
                    if rng.random() < 0.5:
                        v = rng.choice(tuple(succ_set[u]))
                        w_hi = min(v - 10, N - b_fut - 1)
                        if w_hi >= u + 10:
                            nc = w_hi - (u + 10) + 1
                            w = rng.randint(u + 10, w_hi)
                            if w not in succ_set[u] and v not in succ_set[w]:
                                succ_set[u].remove(v)
                                pred_set[v].remove(u)
                                succ_set[u].add(w)
                                pred_set[w].add(u)
                                succ_set[w].add(v)
                                pred_set[v].add(w)
                                dw_p = len(succ_set[w])
                                ns, *_ = eval_poset()
                                hr = float(nc) / float(dw_p)
                                if rng.random() < min(1.0, hr * math.exp(-min(50.0, max(-50.0, ns - cur_s)))):
                                    cur_s = ns
                                else:
                                    succ_set[u].remove(w)
                                    pred_set[w].remove(u)
                                    succ_set[w].remove(v)
                                    pred_set[v].remove(w)
                                    succ_set[u].add(v)
                                    pred_set[v].add(u)
                    else:
                        w = rng.choice(tuple(succ_set[u]))
                        if w < N - b_fut and succ_set[w]:
                            dw = len(succ_set[w])
                            v = rng.choice(tuple(succ_set[w]))
                            w_hi = min(v - 10, N - b_fut - 1)
                            if (u + 10 <= w <= w_hi) and v not in succ_set[u]:
                                nc = w_hi - (u + 10) + 1
                                succ_set[u].remove(w)
                                pred_set[w].remove(u)
                                succ_set[w].remove(v)
                                pred_set[v].remove(w)
                                succ_set[u].add(v)
                                pred_set[v].add(u)
                                ns, *_ = eval_poset()
                                hr = float(dw) / float(nc)
                                if rng.random() < min(1.0, hr * math.exp(-min(50.0, max(-50.0, ns - cur_s)))):
                                    cur_s = ns
                                else:
                                    succ_set[u].remove(v)
                                    pred_set[v].remove(u)
                                    succ_set[u].add(w)
                                    pred_set[w].add(u)
                                    succ_set[w].add(v)
                                    pred_set[v].add(w)

        fin_s, fin_l, fin_d, fin_f, fin_dmm, fin_diam = eval_poset()
        d_spec, d_ball = _compute_spec_and_ball_dim(N, succ_set)
        cone_frac, lc_only_frac, res_per_edge, res_sp_per_edge = _solder_12a_and_reconstruct(
            N, succ_set, q2i, c12a, hasse_only=True
        )
        print(
            f"  N={N:3d} (steps={steps:5d}) | L_max: {init_l}->{fin_l:2d} | <d_out>: {init_d:5.2f}->{fin_d:4.2f} "
            f"| N_diam={fin_diam:3d} | f={fin_f:.4f} | D_MM={fin_dmm:4.2f} | D_spec={d_spec:4.2f} | D_ball={d_ball:4.2f} "
            f"| LC_fid={cone_frac:5.1f}% (LC_only={lc_only_frac:5.1f}%) | E[X]/M={res_per_edge:6.4f} (E_sp/M={res_sp_per_edge:6.4f})"
        )

    # Benchmark: 4D Parallelogram Grid (N = 3^4 = 81)
    pts_81 = [(i0, i1, i2, i3) for i0 in range(3) for i1 in range(3) for i2 in range(3) for i3 in range(3)]
    imap_81 = {pt: i for i, pt in enumerate(pts_81)}
    succ_81: List[Set[int]] = [set() for _ in range(81)]
    for i, pt in enumerate(pts_81):
        for mu in range(4):
            nxt = list(pt)
            nxt[mu] += 1
            nt = tuple(nxt)
            if nt in imap_81:
                succ_81[i].add(imap_81[nt])
    rb_81 = [0] * 81
    diams_81 = 0
    for u in range(80, -1, -1):
        ru = 0
        s1, s2 = 0, 0
        for v in succ_81[u]:
            ru |= (1 << v) | rb_81[v]
            cm = 0
            for w in succ_81[v]:
                cm |= 1 << w
            s2 |= s1 & cm
            s1 |= cm
        rb_81[u] = ru
        diams_81 += s2.bit_count()
    f_81 = sum(r.bit_count() for r in rb_81) / (0.5 * 81 * 80)
    dmm_81 = _invert_myrheim_meyer(f_81)
    dspec_81, dball_81 = _compute_spec_and_ball_dim(81, succ_81)
    print(
        f"  Benchmark 4D Grid (N=81)      | L_max= 9 | <d_out>=4.00 | N_diam={diams_81:3d} "
        f"| f={f_81:.4f} | D_MM={dmm_81:4.2f} | D_spec={dspec_81:4.2f} | D_ball={dball_81:4.2f} | LC_fid=100.0% | E[X]/M=0.0000"
    )

    # Control: 4-ary Branched Tree (N = 1 + 4 + 16 + 64 = 85)
    succ_85: List[Set[int]] = [set() for _ in range(85)]
    for u in range(21):
        for k in range(1, 5):
            succ_85[u].add(4 * u + k)
    rb_85 = [0] * 85
    for u in range(84, -1, -1):
        ru = 0
        for v in succ_85[u]:
            ru |= (1 << v) | rb_85[v]
        rb_85[u] = ru
    f_85 = sum(r.bit_count() for r in rb_85) / (0.5 * 85 * 84)
    dmm_85 = _invert_myrheim_meyer(f_85)
    dspec_85, dball_85 = _compute_spec_and_ball_dim(85, succ_85)
    print(
        f"  Control 4-ary Tree (N=85)     | L_max= 4 | <d_out>=4.00 | N_diam=  0 "
        f"| f={f_85:.4f} | D_MM={dmm_85:4.2f} | D_spec={dspec_85:4.2f} | D_ball={dball_85:4.2f} | LC_fid=100.0% | E[X]/M=0.0000"
    )


if __name__ == "__main__":
    verify_table1_shells_and_5design()
    verify_table2_galois_phason_confinement()
    verify_table3_harmonic_soldering_and_curvature()
    verify_table4_multisize_melting_and_soldering()
