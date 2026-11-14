"""Character table solver, Real Wedderburn-Artin decomposition, and Cayley graph spectral analysis."""

from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Dict, List, Tuple
from src.finite_groups import FiniteGroup


@dataclass
class IrrepCharacter:
    """An irreducible complex character of a finite group G."""

    index: int
    degree: int
    frobenius_schur: int  # +1 (real R), 0 (complex C), -1 (quaternionic H)
    class_values: List[complex]
    is_faithful: bool

    @property
    def division_algebra(self) -> str:
        if self.frobenius_schur == 1:
            return "R"
        elif self.frobenius_schur == 0:
            return "C"
        elif self.frobenius_schur == -1:
            return "H"
        raise ValueError(f"Invalid Frobenius-Schur indicator: {self.frobenius_schur}")


@dataclass
class WedderburnBlock:
    """A simple matrix algebra component M_m(D) inside the real group algebra R[G]."""

    matrix_size: int
    division_algebra: str  # "R", "C", or "H"
    real_dimension: int
    source_irrep_degrees: List[int]

    def __str__(self) -> str:
        if self.matrix_size == 1:
            return f"{self.division_algebra}"
        return f"M_{self.matrix_size}({self.division_algebra})"


@dataclass
class CharacterTableResult:
    group_name: str
    order: int
    class_sizes: List[int]
    class_orders: List[int]
    irreps: List[IrrepCharacter]
    max_orthogonality_error: float
    wedderburn_blocks: List[WedderburnBlock]

    def format_wedderburn_summary(self) -> str:
        counts: Dict[str, int] = {}
        for b in self.wedderburn_blocks:
            s = str(b)
            counts[s] = counts.get(s, 0) + 1
        parts = []
        for s, c in counts.items():
            parts.append(f"{s}^{c}" if c > 1 else s)
        return " (+) ".join(parts)

    def has_connes_algebra_components(self) -> Dict[str, int]:
        """Count blocks matching C (or R), H, M_3(C), and quark-like M_6(C) in R[G]."""
        c_blocks = sum(
            1
            for b in self.wedderburn_blocks
            if b.division_algebra in ("C", "R") and b.matrix_size == 1
        )
        h_blocks = sum(
            1
            for b in self.wedderburn_blocks
            if b.division_algebra == "H" and b.matrix_size == 1
        )
        m3c_blocks = sum(
            1
            for b in self.wedderburn_blocks
            if b.division_algebra == "C" and b.matrix_size == 3
        )
        m6c_blocks = sum(
            1
            for b in self.wedderburn_blocks
            if b.division_algebra == "C" and b.matrix_size == 6
        )
        return {
            "scalar_1d": c_blocks,
            "H_weak": h_blocks,
            "M3C_color": m3c_blocks,
            "M6C_quark_bifund": m6c_blocks,
        }


def solve_character_table(group: FiniteGroup) -> CharacterTableResult:
    """Compute the exact complex character table, Frobenius-Schur indicators, and R[G] decomposition."""
    n = group.order
    classes = group.classes
    k = len(classes)
    elem_to_class = group.elem_to_class
    mult = group.mult

    # Class multiplication coefficients N[r][j][m]
    N = [[[0] * k for _ in range(k)] for _ in range(k)]
    for r in range(k):
        for j in range(k):
            for x in classes[r]:
                row_x = mult[x]
                for y in classes[j]:
                    m = elem_to_class[row_x[y]]
                    N[r][j][m] += 1
            for m in range(k):
                N[r][j][m] //= len(classes[m])

    # Build generic Hermitian linear combination M of commuting normal matrices
    # T^(r)_{jm} = sqrt(|C_m| / |C_j|) N[r][j][m]
    M = [[0j] * k for _ in range(k)]
    for r in range(k):
        cr = math.sin(r + 1.123456789) * 3.14159265 + 0.718281828 * ((r + 1) ** 1.5)
        dr = math.cos(r + 2.987654321) * 2.71828182 + 0.618033988 * ((r + 1) ** 1.7)
        for j in range(k):
            len_j = len(classes[j])
            for m in range(k):
                len_m = len(classes[m])
                t_jm = math.sqrt(len_m / len_j) * N[r][j][m]
                t_mj = math.sqrt(len_j / len_m) * N[r][m][j]
                sym = 0.5 * (t_jm + t_mj)
                skew = 0.5 * (t_jm - t_mj)
                M[j][m] += cr * sym + 1j * dr * skew

    # Complex Jacobi unitary diagonalization of Hermitian matrix M
    V = [[0j] * k for _ in range(k)]
    for i in range(k):
        V[i][i] = 1.0 + 0j
    A = [row[:] for row in M]

    for _sweep in range(80):
        max_off = 0.0
        for p in range(k):
            for q in range(p + 1, k):
                apq = abs(A[p][q])
                if apq > max_off:
                    max_off = apq
                if apq > 1e-15:
                    phase = A[p][q] / apq
                    app = A[p][p].real
                    aqq = A[q][q].real
                    tau = (aqq - app) / (2.0 * apq)
                    t = (
                        math.copysign(
                            1.0 / (abs(tau) + math.sqrt(1.0 + tau * tau)), tau
                        )
                        if tau != 0.0
                        else 1.0
                    )
                    c = 1.0 / math.sqrt(1.0 + t * t)
                    s = t * c
                    phase_conj = phase.conjugate()

                    for idx in range(k):
                        aip = A[idx][p]
                        aiq = A[idx][q]
                        A[idx][p] = c * aip - s * phase_conj * aiq
                        A[idx][q] = s * phase * aip + c * aiq
                    for idx in range(k):
                        api = A[p][idx]
                        aqi = A[q][idx]
                        A[p][idx] = c * api - s * phase * aqi
                        A[q][idx] = s * phase_conj * api + c * aqi
                    for idx in range(k):
                        vip = V[idx][p]
                        viq = V[idx][q]
                        V[idx][p] = c * vip - s * phase_conj * viq
                        V[idx][q] = s * phase * vip + c * viq
        if max_off < 1e-14:
            break

    # Extract irreducible characters and Frobenius-Schur indicators
    raw_irreps = []
    for col in range(k):
        v0 = V[0][col]
        u0 = v0 / abs(v0)
        chi: List[complex] = []
        for j in range(k):
            val = math.sqrt(n / len(classes[j])) * (V[j][col] / u0).conjugate()
            # Clean tiny numerical noise
            re_v = 0.0 if abs(val.real) < 1e-11 else val.real
            im_v = 0.0 if abs(val.imag) < 1e-11 else val.imag
            chi.append(complex(re_v, im_v))

        deg = round(chi[0].real)
        # Frobenius-Schur indicator: iota(chi) = (1/|G|) sum_{g in G} chi(g^2)
        fs_sum = 0j
        for g in range(n):
            g2 = mult[g][g]
            fs_sum += chi[elem_to_class[g2]]
        fs_ind = round((fs_sum / n).real)

        # Check faithfulness: kernel is {g : chi(g) == deg}
        is_faithful = all(
            abs(chi[j] - deg) > 1e-7 for j in range(1, k)
        )
        raw_irreps.append((deg, fs_ind, chi, is_faithful))

    # Sort irreps by degree, then by Frobenius-Schur indicator (+1, 0, -1), then by character values
    raw_irreps.sort(
        key=lambda x: (
            x[0],
            -x[1],
            tuple(-round(c.real, 5) for c in x[2]),
            tuple(-round(c.imag, 5) for c in x[2]),
        )
    )

    irreps = [
        IrrepCharacter(
            index=i,
            degree=deg,
            frobenius_schur=fs_ind,
            class_values=chi,
            is_faithful=is_f,
        )
        for i, (deg, fs_ind, chi, is_f) in enumerate(raw_irreps)
    ]

    # Verify character orthogonality
    max_err = 0.0
    for a in range(k):
        for b in range(k):
            ip = (
                sum(
                    len(classes[j])
                    * irreps[a].class_values[j]
                    * irreps[b].class_values[j].conjugate()
                    for j in range(k)
                )
                / n
            )
            target = 1.0 if a == b else 0.0
            err = abs(ip - target)
            if err > max_err:
                max_err = err

    assert sum(ir.degree**2 for ir in irreps) == n, "Degree sum-of-squares mismatch!"
    assert max_err < 1e-8, f"Character orthogonality error too high: {max_err}"

    wedderburn_blocks = compute_wedderburn_blocks(irreps)
    class_sizes = [len(c) for c in classes]
    class_orders = [group.elem_orders[c[0]] for c in classes]

    return CharacterTableResult(
        group_name=group.name,
        order=n,
        class_sizes=class_sizes,
        class_orders=class_orders,
        irreps=irreps,
        max_orthogonality_error=max_err,
        wedderburn_blocks=wedderburn_blocks,
    )


def compute_wedderburn_blocks(irreps: List[IrrepCharacter]) -> List[WedderburnBlock]:
    """Compute the Real Wedderburn-Artin decomposition of R[G] from its IrrepCharacters."""
    blocks: List[WedderburnBlock] = []
    complex_irreps: List[IrrepCharacter] = []

    for ir in irreps:
        d = ir.degree
        if ir.frobenius_schur == 1:
            blocks.append(
                WedderburnBlock(
                    matrix_size=d,
                    division_algebra="R",
                    real_dimension=d * d,
                    source_irrep_degrees=[d],
                )
            )
        elif ir.frobenius_schur == -1:
            assert d % 2 == 0, "Quaternionic irrep must have even complex degree"
            m = d // 2
            blocks.append(
                WedderburnBlock(
                    matrix_size=m,
                    division_algebra="H",
                    real_dimension=4 * m * m,
                    source_irrep_degrees=[d],
                )
            )
        elif ir.frobenius_schur == 0:
            complex_irreps.append(ir)

    # Complex irreps (iota = 0) pair up with their complex conjugates {chi, bar(chi)}
    assert len(complex_irreps) % 2 == 0, "Complex irreps with iota=0 must come in pairs"
    by_deg: Dict[int, int] = {}
    for ir in complex_irreps:
        by_deg[ir.degree] = by_deg.get(ir.degree, 0) + 1
    for d, cnt in sorted(by_deg.items()):
        assert cnt % 2 == 0
        for _ in range(cnt // 2):
            blocks.append(
                WedderburnBlock(
                    matrix_size=d,
                    division_algebra="C",
                    real_dimension=2 * d * d,
                    source_irrep_degrees=[d, d],
                )
            )

    alg_Order = {"R": 0, "C": 1, "H": 2}
    blocks.sort(key=lambda b: (alg_Order[b.division_algebra], b.matrix_size))
    return blocks


def solve_product_wedderburn(
    res1: CharacterTableResult, res2: CharacterTableResult
) -> CharacterTableResult:
    """Compute the exact character table and Real Wedderburn decomposition of G1 x G2.

    Uses the exact tensor-product identities:
      deg(chi_1 (x) chi_2) = deg(chi_1) * deg(chi_2)
      iota(chi_1 (x) chi_2) = iota(chi_1) * iota(chi_2)
    """
    prod_irreps: List[IrrepCharacter] = []
    idx = 0
    for ir1 in res1.irreps:
        for ir2 in res2.irreps:
            d = ir1.degree * ir2.degree
            fs = ir1.frobenius_schur * ir2.frobenius_schur
            prod_irreps.append(
                IrrepCharacter(
                    index=idx,
                    degree=d,
                    frobenius_schur=fs,
                    class_values=[],
                    is_faithful=ir1.is_faithful and ir2.is_faithful,
                )
            )
            idx += 1

    prod_irreps.sort(key=lambda x: (x.degree, -x.frobenius_schur))
    blocks = compute_wedderburn_blocks(prod_irreps)
    total_order = res1.order * res2.order
    assert sum(b.real_dimension for b in blocks) == total_order

    return CharacterTableResult(
        group_name=f"{res1.group_name} x {res2.group_name}",
        order=total_order,
        class_sizes=[s1 * s2 for s1 in res1.class_sizes for s2 in res2.class_sizes],
        class_orders=[
            math.lcm(o1, o2) for o1 in res1.class_orders for o2 in res2.class_orders
        ],
        irreps=prod_irreps,
        max_orthogonality_error=max(
            res1.max_orthogonality_error, res2.max_orthogonality_error
        ),
        wedderburn_blocks=blocks,
    )


@dataclass
class CayleySpectrumResult:
    group_name: str
    gen_name: str
    order: int
    degree: int
    triangles_w3: int
    plaquettes_w4: int
    kappa_adj: float
    kappa_lap: float
    bare_mh_over_mw_adj: float
    bare_mh_over_mw_lap: float
    distinct_adj_eigenvalues: List[Tuple[float, int]]  # (eigenvalue, multiplicity)


def analyze_cayley_spectrum(
    group: FiniteGroup,
    char_res: CharacterTableResult,
    S: List[int],
    gen_name: str,
) -> CayleySpectrumResult:
    """Compute exact closed-walk counts, spectral moment ratios kappa, and eigenvalues for Cay(G, S)."""
    n = group.order
    mult = group.mult
    inv = group.inv
    S_set = set(S)
    assert 0 not in S_set, "Identity should not be in generating set S"
    for s in list(S_set):
        assert inv[s] in S_set, "Generating set S must be symmetric (inverse-closed)"
    S_list = sorted(S_set)
    k_s = len(S_list)

    # Compute exact closed walk counts from e=0 of lengths 2, 3, 4
    # distribution after 1 step: dist1[s] = 1 for s in S
    dist1 = {s: 1 for s in S_list}
    dist2: Dict[int, int] = {}
    for x, c1 in dist1.items():
        row_x = mult[x]
        for s in S_list:
            y = row_x[s]
            dist2[y] = dist2.get(y, 0) + c1

    w2 = dist2.get(0, 0)
    assert w2 == k_s

    # W_3 = sum_{g} dist2[g] * 1(g^{-1} in S)
    w3 = sum(cnt for g, cnt in dist2.items() if inv[g] in S_set)
    # W_4 = sum_{g} dist2[g] * dist2[g^{-1}] (since S = S^{-1}, dist2[g^{-1}] == dist2[g])
    w4 = sum(cnt * dist2.get(inv[g], 0) for g, cnt in dist2.items())

    # Adjacency moments Tr(A^m) = |G| * W_m
    # kappa(A_S) = |G| Tr(A^4) / (Tr(A^2))^2 = W_4 / k_S^2
    kappa_adj = w4 / (k_s * k_s)

    # Laplacian Delta_S = k_S I - A_S:
    # (1/|G|) Tr(Delta^2) = k_S^2 + W_2 = k_S(k_S + 1)
    # (1/|G|) Tr(Delta^4) = k_S^4 + 6 k_S^2 W_2 - 4 k_S W_3 + W_4
    lap_m2 = k_s * (k_s + 1)
    lap_m4 = (k_s**4) + 6 * (k_s**3) - 4 * k_s * w3 + w4
    kappa_lap = lap_m4 / (lap_m2 * lap_m2)

    # Check if S is class-symmetric (union of conjugacy classes)
    is_class_symmetric = all(
        all(c_el in S_set for c_el in group.classes[group.elem_to_class[s]])
        for s in S_list
    )

    eig_Summary: List[Tuple[float, int]] = []
    if is_class_symmetric:
        # Exact eigenvalues via Schur's lemma: lambda_chi = (1/d_chi) sum_{s in S} chi(s)
        eig_counts: Dict[float, int] = {}
        for ir in char_res.irreps:
            s_sum = sum(
                ir.class_values[group.elem_to_class[s]].real for s in S_list
            )
            lam = round(s_sum / ir.degree, 6)
            if abs(lam) < 1e-9:
                lam = 0.0
            eig_counts[lam] = eig_counts.get(lam, 0) + (ir.degree * ir.degree)
        eig_Summary = sorted(eig_counts.items(), key=lambda x: -x[0])

    return CayleySpectrumResult(
        group_name=group.name,
        gen_name=gen_name,
        order=n,
        degree=k_s,
        triangles_w3=w3,
        plaquettes_w4=w4,
        kappa_adj=kappa_adj,
        kappa_lap=kappa_lap,
        bare_mh_over_mw_adj=2.0 * math.sqrt(kappa_adj),
        bare_mh_over_mw_lap=2.0 * math.sqrt(kappa_lap),
        distinct_adj_eigenvalues=eig_Summary,
    )
