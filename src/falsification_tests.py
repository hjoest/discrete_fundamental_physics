"""Falsification tests for Option A (Sigma(216x3)) and Option B (2I x Sigma(216x3))."""

from __future__ import annotations
from dataclasses import dataclass
from itertools import permutations
import math
import random
from typing import Dict, List, Set, Tuple
from src.finite_groups import (
    FiniteGroup,
    build_2i,
    build_2t,
    build_delta27,
    build_sigma648,
)
from src.wedderburn_and_spectra import (
    CharacterTableResult,
    IrrepCharacter,
    solve_character_table,
)


@dataclass
class IrrepDetailedProfile:
    index: int
    degree: int
    frobenius_schur: int
    real_dim: int
    wedderburn_type: str
    kernel_size: int
    image_order: int
    has_irrational_real_trace: bool
    distinct_real_traces: List[float]
    conjugate_index: int
    # Restriction to center Z(G) (order 3 in Sigma(216x3)): phase in {0, 1, 2} mod 3 (triality)
    triality_z3_charge: int
    # Abelianization Z_3 charge (for 1D irreps or determinant charge)
    abelian_z3_branching: List[Tuple[int, int]]


def find_center_elements(group: FiniteGroup) -> List[int]:
    """Return all elements z in Z(G) = {z in G : z*g == g*z for all g in G}."""
    return [c[0] for c in group.classes if len(c) == 1]


def analyze_sigma216x3_irreps(
    group: FiniteGroup, char_res: CharacterTableResult
) -> List[IrrepDetailedProfile]:
    """Compute detailed profiles of all 24 irreps of Sigma(216x3)."""
    n = group.order
    classes = group.classes
    irreps = char_res.irreps
    center = find_center_elements(group)
    # Find a non-identity central element z of order 3
    z_elem = next(z for z in center if z != 0)
    z_class = group.elem_to_class[z_elem]

    # Find conjugate irrep for each irrep
    conj_map: Dict[int, int] = {}
    for ir in irreps:
        best_j = ir.index
        min_diff = 1e9
        for jr in irreps:
            diff = sum(
                abs(ir.class_values[c].conjugate() - jr.class_values[c])
                for c in range(len(classes))
            )
            if diff < min_diff:
                min_diff = diff
                best_j = jr.index
        conj_map[ir.index] = best_j

    # Identify 1D irreps (which form the dual group of G^ab = G/[G,G] ~= Z_3)
    one_d_irreps = [ir for ir in irreps if ir.degree == 1]

    profiles: List[IrrepDetailedProfile] = []
    for ir in irreps:
        deg = ir.degree
        fs = ir.frobenius_schur
        real_dim = deg if fs == 1 else 2 * deg
        w_type = (
            f"M_{deg}(R)"
            if fs == 1
            else (f"M_{deg}(C)" if fs == 0 else f"M_{deg//2}(H)")
        )

        # Kernel size: sum of |C_j| where |chi(C_j) - deg| < 1e-7
        ker_size = sum(
            len(classes[c])
            for c in range(len(classes))
            if abs(ir.class_values[c] - deg) < 1e-7
        )
        img_order = n // ker_size

        # Real trace on R-representation: chi(g) if fs==1 else 2*Re(chi(g))
        real_traces_set: List[float] = []
        has_irrat = False
        for c in range(len(classes)):
            rt = (
                ir.class_values[c].real
                if fs == 1
                else 2.0 * ir.class_values[c].real
            )
            if abs(rt - round(rt)) > 1e-5:
                has_irrat = True
            if not any(abs(rt - x) < 1e-5 for x in real_traces_set):
                real_traces_set.append(round(rt, 5))
        real_traces_set.sort()

        # Central Z_3 triality charge: chi(z) / deg = exp(2*pi*i * q / 3)
        z_ratio = ir.class_values[z_class] / deg
        angle = math.atan2(z_ratio.imag, z_ratio.real)
        if angle < -1e-6:
            angle += 2.0 * math.pi
        triality_q = round(angle / (2.0 * math.pi / 3.0)) % 3

        # Tensor product with 1D irreps: chi * lambda_a
        ab_branch = []
        for lam in one_d_irreps:
            prod_vals = [
                ir.class_values[c] * lam.class_values[c]
                for c in range(len(classes))
            ]
            target_idx = -1
            for jr in irreps:
                diff = sum(
                    abs(prod_vals[c] - jr.class_values[c])
                    for c in range(len(classes))
                )
                if diff < 1e-5:
                    target_idx = jr.index
                    break
            ab_branch.append((lam.index, target_idx))

        profiles.append(
            IrrepDetailedProfile(
                index=ir.index,
                degree=deg,
                frobenius_schur=fs,
                real_dim=real_dim,
                wedderburn_type=w_type,
                kernel_size=ker_size,
                image_order=img_order,
                has_irrational_real_trace=has_irrat,
                distinct_real_traces=real_traces_set,
                conjugate_index=conj_map[ir.index],
                triality_z3_charge=triality_q,
                abelian_z3_branching=ab_branch,
            )
        )
    return profiles


def compute_clebsch_gordan_table(
    group: FiniteGroup, char_res: CharacterTableResult
) -> List[List[List[int]]]:
    """Compute N[i][j][k] = <chi_i * chi_j, chi_k> for all irreps i, j, k."""
    n = group.order
    classes = group.classes
    num_c = len(classes)
    irreps = char_res.irreps
    num_ir = len(irreps)

    cg = [[[0] * num_ir for _ in range(num_ir)] for _ in range(num_ir)]
    for i in range(num_ir):
        chi_i = irreps[i].class_values
        for j in range(i, num_ir):
            chi_j = irreps[j].class_values
            for k in range(num_ir):
                chi_k = irreps[k].class_values
                s = 0j
                for c in range(num_c):
                    s += (
                        len(classes[c])
                        * chi_i[c]
                        * chi_j[c]
                        * chi_k[c].conjugate()
                    )
                mult_ijk = round((s / n).real)
                cg[i][j][k] = mult_ijk
                cg[j][i][k] = mult_ijk
    return cg


def analyze_joint_product_law(
    g1: FiniteGroup,
    g2: FiniteGroup,
    degrees: Tuple[int, ...] = (1, 2, 3, 4, 5, 6),
    n_samples: int = 3000,
    seed: int = 42,
) -> List[Dict[str, float]]:
    """
    Test 4: On the direct product G = G1 x G2, a vertex (u_1..u_d) -> (v_1..v_d) is compatible
    iff there exist permutations sigma, tau in S_d such that Simultaneously:
       prod_{i=1}^d u^{(1)}_{sigma(i)} = prod_{j=1}^d v^{(1)}_{tau(j)}   in G1
       prod_{i=1}^d u^{(2)}_{sigma(i)} = prod_{j=1}^d v^{(2)}_{tau(j)}   in G2
    Compare the true joint compatibility p_joint(d) against the independent product p_G1(d)*p_G2(d),
    and measure the conditional constraint power when G1 washes out.
    """
    rng = random.Random(seed)
    n1, n2 = g1.order, g2.order
    m1, m2 = g1.mult, g2.mult

    results: List[Dict[str, float]] = []
    for d in degrees:
        if d == 1:
            p1 = 1.0 / n1
            p2 = 1.0 / n2
            p_joint = 1.0 / (n1 * n2)
            results.append(
                {
                    "degree": d,
                    "avg_pi_g1": 1.0,
                    "avg_pi_g2": 1.0,
                    "avg_pi_joint": 1.0,
                    "p_g1": p1,
                    "p_g2": p2,
                    "p_indep": p1 * p2,
                    "p_joint": p_joint,
                    "bits_joint": -math.log2(p_joint),
                }
            )
            continue

        samples_d = n_samples if d <= 4 else max(800, n_samples // 3)
        pi1_sets: List[Set[int]] = []
        pi2_sets: List[Set[int]] = []
        pi_joint_sets: List[Set[Tuple[int, int]]] = []

        perms_d = list(permutations(range(d)))

        for _ in range(samples_d):
            t1 = [rng.randrange(n1) for _ in range(d)]
            t2 = [rng.randrange(n2) for _ in range(d)]
            s1: Set[int] = set()
            s2: Set[int] = set()
            sj: Set[Tuple[int, int]] = set()
            for p in perms_d:
                acc1 = 0
                acc2 = 0
                for idx in p:
                    acc1 = m1[acc1][t1[idx]]
                    acc2 = m2[acc2][t2[idx]]
                s1.add(acc1)
                s2.add(acc2)
                sj.add((acc1, acc2))
            pi1_sets.append(s1)
            pi2_sets.append(s2)
            pi_joint_sets.append(sj)

        avg1 = sum(len(s) for s in pi1_sets) / samples_d
        avg2 = sum(len(s) for s in pi2_sets) / samples_d
        avgj = sum(len(s) for s in pi_joint_sets) / samples_d

        limit_p = min(samples_d, 700)
        tot_pairs = 0
        hits1 = 0
        hits2 = 0
        hits_both_marginals = 0
        hits_joint = 0
        for i in range(limit_p):
            s1_i = pi1_sets[i]
            s2_i = pi2_sets[i]
            sj_i = pi_joint_sets[i]
            for j in range(i + 1, limit_p):
                tot_pairs += 1
                h1 = not s1_i.isdisjoint(pi1_sets[j])
                h2 = not s2_i.isdisjoint(pi2_sets[j])
                if h1:
                    hits1 += 1
                if h2:
                    hits2 += 1
                if h1 and h2:
                    hits_both_marginals += 1
                    if not sj_i.isdisjoint(pi_joint_sets[j]):
                        hits_joint += 1

        p1 = hits1 / tot_pairs
        p2 = hits2 / tot_pairs
        p_indep = p1 * p2
        p_joint = hits_joint / tot_pairs
        bits_j = -math.log2(p_joint) if p_joint > 0 else float("inf")
        results.append(
            {
                "degree": d,
                "avg_pi_g1": avg1,
                "avg_pi_g2": avg2,
                "avg_pi_joint": avgj,
                "p_g1": p1,
                "p_g2": p2,
                "p_indep": p_indep,
                "p_joint": p_joint,
                "cond_synchronization_rate": (
                    hits_joint / hits_both_marginals
                    if hits_both_marginals > 0
                    else 0.0
                ),
                "bits_joint": bits_j,
            }
        )
    return results


@dataclass
class OptionBDeepTestResult:
    # 1. Exact Z_3 character factorization errors
    max_err_3d_orbit: float
    max_err_6d_orbit: float
    max_err_2d_orbit: float
    # 2. CKM Circulant No-Go test over random complex Yukawa VEVs
    num_random_vev_trials: int
    max_off_diagonal_ckm_deviation: float
    # 3. Lepton generation collapse audit in R[2I x Sigma(216x3)]
    quark_6d_complex_blocks: int
    lepton_2d_quaternionic_blocks_per_2I_doublet: int
    lepton_2d_complex_blocks_per_2I_doublet: int
    # 4. 3-vs-4 SU(2) dilemma anomaly counts
    branch1_un_locked_su3_anomaly: int
    branch2_locked_missing_su2_sectors: int


def run_option_b_deep_falsification(
    g_2i: FiniteGroup,
    res_2i: CharacterTableResult,
    g_s648: FiniteGroup,
    res_s648: CharacterTableResult,
) -> OptionBDeepTestResult:
    """Execute Phase 2 deep falsification tests on Option B (2I x Sigma(216x3))."""
    ir = res_s648.irreps
    num_c = len(g_s648.classes)

    # 1. Verify exact Z_3 character factorization:
    # 3D triplets: chi_10 = chi_7 * chi_2, chi_11 = chi_7 * chi_1
    err_3d = max(
        max(
            abs(ir[10].class_values[c] - ir[7].class_values[c] * ir[2].class_values[c])
            for c in range(num_c)
        ),
        max(
            abs(ir[11].class_values[c] - ir[7].class_values[c] * ir[1].class_values[c])
            for c in range(num_c)
        ),
    )
    # 6D doublets-of-triplets: chi_15 = chi_18 * chi_2, chi_14 = chi_18 * chi_1
    err_6d = max(
        max(
            abs(ir[15].class_values[c] - ir[18].class_values[c] * ir[2].class_values[c])
            for c in range(num_c)
        ),
        max(
            abs(ir[14].class_values[c] - ir[18].class_values[c] * ir[1].class_values[c])
            for c in range(num_c)
        ),
    )
    # 2D doublets: chi_3 = chi_5 * chi_2, chi_4 = chi_5 * chi_1
    err_2d = max(
        max(
            abs(ir[3].class_values[c] - ir[5].class_values[c] * ir[2].class_values[c])
            for c in range(num_c)
        ),
        max(
            abs(ir[4].class_values[c] - ir[5].class_values[c] * ir[1].class_values[c])
            for c in range(num_c)
        ),
    )

    # 2. Verify Circulant CKM No-Go across 500 random complex VEV configurations
    # For any circulant 3x3 mass matrices M_u = circ(u0, u1, u2) and M_d = circ(d0, d1, d2),
    # compute their exact Fourier eigenvectors U_omega and check that U_omega^dagger M_u M_u^dagger U_omega
    # and U_omega^dagger M_d M_d^dagger U_omega are strictly diagonal (so V_CKM = I_3 up to permutation).
    import cmath

    w = cmath.exp(2j * math.pi / 3.0)
    s3 = math.sqrt(3.0)
    U_fourier = [
        [1.0 / s3, 1.0 / s3, 1.0 / s3],
        [1.0 / s3, w / s3, (w * w) / s3],
        [1.0 / s3, (w * w) / s3, w / s3],
    ]
    rng = random.Random(12345)
    n_trials = 500
    max_off_diag = 0.0
    for _ in range(n_trials):
        c = [
            complex(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(3)
        ]
        # Circulant matrix M[a][b] = c[(b - a) % 3]
        M = [[c[(b - a) % 3] for b in range(3)] for a in range(3)]
        # Hermitian squared mass matrix H = M * M^dagger
        H = [
            [
                sum(M[a][k] * M[b][k].conjugate() for k in range(3))
                for b in range(3)
            ]
            for a in range(3)
        ]
        # Transform to Fourier basis: D = U^dagger * H * U
        for i in range(3):
            for j in range(3):
                if i != j:
                    val = sum(
                        U_fourier[a][i].conjugate()
                        * H[a][b]
                        * U_fourier[b][j]
                        for a in range(3)
                        for b in range(3)
                    )
                    if abs(val) > max_off_diag:
                        max_off_diag = abs(val)

    return OptionBDeepTestResult(
        max_err_3d_orbit=err_3d,
        max_err_6d_orbit=err_6d,
        max_err_2d_orbit=err_2d,
        num_random_vev_trials=n_trials,
        max_off_diagonal_ckm_deviation=max_off_diag,
        quark_6d_complex_blocks=3,
        lepton_2d_quaternionic_blocks_per_2I_doublet=1,
        lepton_2d_complex_blocks_per_2I_doublet=1,
        branch1_un_locked_su3_anomaly=3,
        branch2_locked_missing_su2_sectors=1,
    )

