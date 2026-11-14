"""Combinatorial and information-theoretic analysis of the Ordered Product Conservation Law."""

from __future__ import annotations
from dataclasses import dataclass
from itertools import permutations
import math
import random
from typing import Dict, List, Set, Tuple
from src.finite_groups import FiniteGroup


@dataclass
class ProductLawStats:
    group_name: str
    order: int
    num_classes: int
    commutator_order: int
    abelianization_order: int
    # Exact m=2 statistics
    avg_pi_size_m2: float
    p_1_2_exact: float
    p_1_2_formula: float
    p_2_2_exact: float
    # Sampled / exact m=3, 4 statistics
    avg_pi_size_m3: float
    p_3_3_est: float
    avg_pi_size_m4: float
    p_4_4_est: float
    # Degree sweep (d = 1..8): (degree, avg_pi_size, coset_fill_fraction, p_d_d, causal_selectivity_bits)
    degree_sweep: List[Tuple[int, float, float, float, float]]


def permutation_product_set(group: FiniteGroup, tup: Tuple[int, ...]) -> Set[int]:
    """Return Pi(g_1, ..., g_m) = { prod_{i=1}^m g_{sigma(i)} : sigma in S_m }."""
    mult = group.mult
    out: Set[int] = set()
    for perm in set(permutations(tup)):
        acc = 0
        for g in perm:
            acc = mult[acc][g]
        out.add(acc)
    return out


def analyze_product_law(
    group: FiniteGroup, n_samples: int = 8000, seed: int = 42
) -> ProductLawStats:
    """Compute exact (m=2) and Monte Carlo (m>=3) statistics of the ordered product law."""
    n = group.order
    k = group.num_classes
    mult = group.mult

    # 1. Exact O(|G|^2) computation for m=2
    # For each (g1, g2) in G^2, Pi(g1, g2) = {g1*g2, g2*g1}
    cnt1 = [0] * n  # cnt1[x] = number of pairs (h1, h2) with x in Pi(h1, h2)
    cnt2: Dict[Tuple[int, int], int] = {}  # unordered pair {a, b} count for a != b
    total_pi_m2 = 0

    pairs_pi: List[Tuple[int, int]] = []
    for g1 in range(n):
        for g2 in range(n):
            a = mult[g1][g2]
            b = mult[g2][g1]
            if a == b:
                total_pi_m2 += 1
                cnt1[a] += 1
                pairs_pi.append((a, a))
            else:
                total_pi_m2 += 2
                cnt1[a] += 1
                cnt1[b] += 1
                pair = (a, b) if a < b else (b, a)
                cnt2[pair] = cnt2.get(pair, 0) + 1
                pairs_pi.append(pair)

    avg_pi_m2 = total_pi_m2 / (n * n)
    p_1_2_exact = total_pi_m2 / (n * n * n)
    p_1_2_formula = (2 * n - k) / (n * n)
    assert abs(p_1_2_exact - p_1_2_formula) < 1e-15

    # Exact (2 -> 2) compatibility via inclusion-exclusion in O(|G|^2)
    compatible_2_2 = 0
    for a, b in pairs_pi:
        if a == b:
            compatible_2_2 += cnt1[a]
        else:
            compatible_2_2 += cnt1[a] + cnt1[b] - cnt2[(a, b)]
    p_2_2_exact = compatible_2_2 / (n**4)

    # 2. Sampled statistics for m = 3..7 (and exact for m=1, 2)
    rng = random.Random(seed)
    comm_size = len(group.commutator_subgroup)

    degree_sweep: List[Tuple[int, float, float, float, float]] = []
    # d = 1
    p_1_1 = 1.0 / n
    degree_sweep.append((1, 1.0, 1.0 / comm_size, p_1_1, -math.log2(p_1_1)))
    # d = 2
    degree_sweep.append(
        (2, avg_pi_m2, avg_pi_m2 / comm_size, p_2_2_exact, -math.log2(p_2_2_exact))
    )

    avg_pi_m3 = 0.0
    p_3_3 = 0.0
    avg_pi_m4 = 0.0
    p_4_4 = 0.0

    for d in (3, 4, 5, 6):
        samples_d = n_samples if d <= 4 else max(1500, n_samples // 4)
        pi_sets: List[Set[int]] = []
        size_sum = 0
        for _ in range(samples_d):
            tup = tuple(rng.randrange(n) for _ in range(d))
            ps = permutation_product_set(group, tup)
            pi_sets.append(ps)
            size_sum += len(ps)
        avg_sz = size_sum / samples_d

        # Estimate p_{d, d} from paired samples
        half = samples_d // 2
        hits = 0
        for i in range(half):
            if not pi_sets[i].isdisjoint(pi_sets[i + half]):
                hits += 1
        # Refine p_{d, d} using cross-pairs for better precision when small
        total_pairs = 0
        hits_all = 0
        limit_pairs = min(samples_d, 800)
        for i in range(limit_pairs):
            ps_i = pi_sets[i]
            for j in range(i + 1, limit_pairs):
                total_pairs += 1
                if not ps_i.isdisjoint(pi_sets[j]):
                    hits_all += 1
        p_dd = hits_all / total_pairs if total_pairs > 0 else 0.0
        info_bits = -math.log2(p_dd) if p_dd > 0 else float("inf")
        degree_sweep.append((d, avg_sz, avg_sz / comm_size, p_dd, info_bits))

        if d == 3:
            avg_pi_m3 = avg_sz
            p_3_3 = p_dd
        elif d == 4:
            avg_pi_m4 = avg_sz
            p_4_4 = p_dd

    return ProductLawStats(
        group_name=group.name,
        order=n,
        num_classes=k,
        commutator_order=comm_size,
        abelianization_order=group.abelianization_order,
        avg_pi_size_m2=avg_pi_m2,
        p_1_2_exact=p_1_2_exact,
        p_1_2_formula=p_1_2_formula,
        p_2_2_exact=p_2_2_exact,
        avg_pi_size_m3=avg_pi_m3,
        p_3_3_est=p_3_3,
        avg_pi_size_m4=avg_pi_m4,
        p_4_4_est=p_4_4,
        degree_sweep=degree_sweep,
    )
