# Deductive Dependency Graph: From Axioms 1.1–1.3 and 2.1–2.2 to the Four-Factor Design Theorem (Theorem 2.15)

This document traces the complete deductive dependency graph from the five foundational axioms (**Axioms 1.1–1.3** in Section 1 and **Axioms 2.1–2.2** in Section 2) through the intermediate propositions and no-go theorems to the **Four-Factor Design Theorem** (**Theorem 2.15**, `\label{thm:design}`) and **Proposition 2.16** (`\label{prop:anomaly_cancellation}`) in [`paper.tex`](paper.tex#L85-L389).

---

## Mermaid Deductive Dependency Graph

```mermaid
flowchart TD
    subgraph S1["Section 1: Foundational Kinematic & Conservation Axioms"]
        A11["Axiom 1.1: Causal Spacetime Poset (P, <=_P)"]
        A12["Axiom 1.2: Finite Group Fiber G & Real Group Algebra R[G]"]
        A13["Axiom 1.3: Total Event Poset Z & Local Ordered Product Law"]

        P14["Prop 1.4: Division-Algebra Decomposition of R[G] into M_n(R), M_n(C), M_n(H)"]
        P15["Prop 1.5: Frobenius-Schur Indicator Classification iota in {+1, 0, -1}"]
        P16["Prop 1.6: Resolution of Covering-Relation Triviality Trap"]
        P17["Prop 1.7: Abelianization Kirchhoff Conservation Law on G^ab = G/[G,G]"]
    end

    subgraph S2A["Section 2: Minimal Holodivisionality & Merolytic Independence Axioms"]
        A21["Axiom 2.1 (Defs 2.1-2.2): G is minimally holodivisional"]
        A22["Axiom 2.2 (Def 2.3): G is merolytic"]
    end

    subgraph S21["Section 2.1: Minimal Faithful Subgroups in U(1), SU(2), and SU(3)"]
        P23["Prop 2.3: Unitary Embedding H in U(d), SU(d)"]
        P24["Prop 2.4: Minimal Unitary Degrees d in {1, 2, 3} -> C + H + M_3(C)"]
        P25["Prop 2.5: Minimal Candidate Classification: Z_3; T_7, Delta(27), M_27; ADE in SU(2)"]
        P26["Prop 2.6: M_3(C) Projective Unitary Operator Basis & M_8(R) Normalizer Closure"]
        T27["Thm 2.7: Elimination of T_7 & M_27; Selection of Z_3 < Delta(27) < Sigma(216x3)"]
        P28["Prop 2.8: Equivalence of Perfection [H,H]=H & Chiral Non-Crystallographic SU(2) Pair"]
        T29["Thm 2.9: Elimination of Solvable C_n, 2D_2n, 2T, 2O; Selection of Perfect 2I"]
    end

    subgraph S22["Section 2.2: Normal-Series Unfolding No-Go Theorems & Four-Factor Design Theorem"]
        P210["Prop 2.10: Frobenius-Schur Multiplicativity & Conjugate-Component Reality Constraint"]
        T211["Thm 2.11 (No-Go 1): Omitting Sigma(216x3) -> No M_8(R) Octet & No Solvable Quaternionic Factor"]
        T212["Thm 2.12 (No-Go 2): Omitting 2I -> Solvable Group, Delta(27) Kernel & 4D Crystal"]
        T213["Thm 2.13 (No-Go 3): Omitting Delta(27)_f -> Quotient-Kernel Collapse 3->2 & Circulant Triviality"]
        T214["Thm 2.14 (No-Go 4): Omitting Z_3^PQ -> Conjugate-Component Reality Lock & Single-Abelianization Obstruction"]
        T215["Thm 2.15 (Four-Factor Design Theorem): Unique Minimal Group G* = 2I x Sigma(216x3) x Delta(27)_f x Z_3^PQ"]
        P216["Prop 2.16: 3-Generation Discrete Ibanez-Ross Anomaly Cancellation (3 | N_gen)"]
    end

    %% Section 1 edges
    A12 --> P14
    P14 --> P15
    A11 --> P16
    A13 --> P16
    A11 --> P17
    A12 --> P17
    A13 --> P17

    %% Section 2.1 edges
    A12 --> P23
    P14 --> P24
    P15 --> P24
    P17 --> P24
    P23 --> P24
    A21 --> P24

    P24 --> P25
    A21 --> P25

    P14 --> P26
    P15 --> P26
    P25 --> P26

    A21 --> T27
    P25 --> T27
    P26 --> T27

    P15 --> P28
    P17 --> P28
    P25 --> P28

    A21 --> T29
    P17 --> T29
    P28 --> T29

    %% Section 2.2 edges
    P15 --> P210

    A21 --> T211
    P26 --> T211
    T27 --> T211

    A21 --> T212
    P28 --> T212
    T29 --> T212

    A22 --> T213
    P210 --> T213
    T27 --> T213

    A22 --> T214
    P210 --> T214
    P25 --> T214

    A21 --> T215
    A22 --> T215
    P25 --> T215
    T27 --> T215
    T29 --> T215
    P210 --> T215
    T211 --> T215
    T212 --> T215
    T213 --> T215
    T214 --> T215

    T215 --> P216
```

---

## Step-by-Step Deductive Ledger

| Node | Statement in [`paper.tex`](paper.tex) | Direct Inputs | Mathematical Output Fed Forward |
| :--- | :--- | :--- | :--- |
| **[Axiom 1.1](paper.tex#L85-L88)** | Causal Spacetime Poset $(P, \le_P)$ | Foundational postulate | Directed covering links $x \prec_P y$ |
| **[Axiom 1.2](paper.tex#L90-L93)** | Finite Group Fiber $G$ and Real Group Algebra $\R[G]$ | Foundational postulate | Finite group $G$ and semisimple real algebra $\R[G]$ |
| **[Axiom 1.3](paper.tex#L95-L107)** | Total Event Poset $Z \subseteq P \times G$ & Ordered Product Law | Foundational postulate | Induced strict order $\prec_Z$ and vertex equation $\prod u_{\sigma(i)} = \prod v_{\tau(j)}$ |
| **[Prop 1.4](paper.tex#L109-L120)** | Division-Algebra Decomposition of $\R[G]$ | Axiom 1.2 (Maschke, Artin–Wedderburn, Frobenius) | $\R[G] \cong \bigoplus M_{n_\alpha}(\mathbb{D}_\alpha)$ with $\mathbb{D}_\alpha \in \{\R, \C, \quat\}$ |
| **[Prop 1.5](paper.tex#L122-L134)** | Frobenius–Schur Classification of Simple Components of $\R[G]$ | Axiom 1.2, Prop 1.4 | $\iota(\chi) \in \{+1, 0, -1\} \leftrightarrow \{\R, \C, \quat\}$ and $d_\chi = 2n_\alpha$ for $\quat$ |
| **[Prop 1.6](paper.tex#L136-L143)** | Resolution of Covering-Relation Triviality Trap | Axioms 1.1, 1.3 | Antichain fibers $Z_x$ and covering equivalence $\prec_Z \iff \prec_P$ |
| **[Prop 1.7](paper.tex#L145-L157)** | Abelianization Kirchhoff Conservation Law | Axioms 1.1, 1.2, 1.3 | Additive charge conservation in $G^{\mathrm{ab}} = G/[G,G]$ ($G^{\mathrm{ab}} \neq \{1\}$ vs. $[K,K]=K$) |
| **[Axiom 2.1 / Defs 2.1–2.2](paper.tex#L165-L192)** | Minimal Holodivisionality ($G$ is *minimally holodivisional*) | Algebraic minimality postulate | Requires $\{\R,\C,\quat\}$ in perfect/solvable & abelian/non-abelian sectors at minimal $(d_\chi, |H|, \text{exp})$ with irreducible $N_{\SU(d_\chi)}(H)$ action on traceless $H/Z(H)$ |
| **[Axiom 2.2 / Def 2.3](paper.tex#L180-L197)** | Merolytic Independence ($G$ is *merolytic*) | Algebraic independence postulate | Terms of $1 \triangleleft Z(H) \triangleleft H \triangleleft N_{\SU(d_\chi)}(H)$ are independent direct factors: (1) $N_{\SU(d_\chi)}(H)/H$ representations act faithfully on an independent factor $\cong H$; (2) conjugate pairs of $N_{\SU(d_\chi)}(H)/H$ are separated by an independent factor $\cong Z(H)$ into disjoint simple Wedderburn components |
| **[Prop 2.3](paper.tex#L208-L215)** | Unitary Embedding of Finite Group Representations | Axiom 1.2 (Weyl's theorem) | Every finite subgroup embeds into $\U(d_\chi)$ or $\SU(d_\chi)$ |
| **[Prop 2.4](paper.tex#L217-L230)** | Minimal Unitary Degrees of Non-Real Division Algebras | Props 1.4, 1.5, 1.7, 2.3, Axiom 2.1 | $d_\chi = 1$ ($\U(1) \to \C$), $d_\chi = 2$ ($\SU(2) \to \quat$), $d_\chi = 3$ ($\SU(3) \to M_3(\C)$) |
| **[Prop 2.5](paper.tex#L232-L245)** | Classification of Minimal Finite Group Candidates | Prop 2.4, Axiom 2.1 | Selects $\Z_3 \subset \U(1)$; candidates $\{T_7, \Delta(27), M_{27}\} \subset \SU(3)$; ADE candidates $\{C_n, 2D_{2n}, 2T, 2O, 2I\} \subset \SU(2)$ |
| **[Prop 2.6](paper.tex#L249-L261)** | $M_3(\C)$ Operator Basis & Irreducible $M_8(\R)$ Normalizer Closure | Props 1.4, 1.5, 2.5 | Equivalence: $9$-element basis $H/Z(H)$ + transitive $N_{\SU(3)}(H)$ action $\iff |H|=27, \text{exp } 3 \iff H \cong \Delta(27), N_{\SU(3)}(H) \cong \Sigma(216\times 3)$ |
| **[Thm 2.7](paper.tex#L263-L270)** | Elimination of $T_7$ and $M_{27}$; Selection of $\Delta(27)$ & $\Sigma(216\times 3)$ | Axiom 2.1, Props 2.5, 2.6 | Rules out $T_7$ and $M_{27}$; yields normal series $\Z_3 \triangleleft \Delta(27) \triangleleft \Sigma(216 \times 3)$ |
| **[Prop 2.8](paper.tex#L274-L286)** | Equivalence of Perfection & Chiral Non-Crystallographic $\SU(2)$ Pair | Props 1.5, 1.7, 2.5 | Equivalence in $\SU(2)$: $[H,H]=H \iff \mathbf{2}_{-1} \oplus \mathbf{2}_{-1}'$ ($5$-design) $\iff H \cong 2I$ |
| **[Thm 2.9](paper.tex#L288-L295)** | Elimination of Solvable $C_n, 2D_{2n}, 2T, 2O$; Selection of $2I$ | Axiom 2.1, Props 1.7, 2.8 | Rules out solvable ADE subgroups; selects unique perfect faithful subgroup $H = 2I \subset \SU(2)$ |
| **[Prop 2.10](paper.tex#L301-L313)** | Frobenius–Schur Multiplicativity & Conjugate-Component Reality | Prop 1.5 | $\iota(\bigboxtimes \chi_j) = \prod \iota(\chi_j)$ and Hermitian conjugation $\rho_{\bar{\chi}}(x) = \rho_\chi(x)^*$ on shared complex components $M_{d_\chi}(\C)$ |
| **[Thm 2.11](paper.tex#L317-L324)** | No-Go 1: Eliminating All Subsets Omitting $\Sigma(216 \times 3)$ | Axiom 2.1, Prop 2.6, Thm 2.7 | Violates Def 2.1(2)–(3): without $\Sigma(216\times 3) = \Delta(27) \rtimes 2T$, no $M_8(\R)$ octet and no solvable quaternionic subgroup independent of $2I$ exist |
| **[Thm 2.12](paper.tex#L326-L333)** | No-Go 2: Eliminating All Subsets Omitting $2I$ | Axiom 2.1, Prop 2.8, Thm 2.9 | Violates Def 2.1(2): without $2I$, $G$ is solvable and all $d_\R \le 4$ irreps have $\Delta(27)$ in their kernel with crystallographic traces |
| **[Thm 2.13](paper.tex#L335-L342)** | No-Go 3: Eliminating All Subsets Omitting $\Delta(27)_f$ | Axiom 2.2, Prop 2.10, Thm 2.7 | Violates Def 2.3(1): on $\Delta(27)$-trivial quotient representations, abelian $\Z_3$ triplets collapse $3 \to 2$ simple Wedderburn components and force $\Z_3$-circulant matrices ($|U_1^\dagger U_2| = I_3$) |
| **[Thm 2.14](paper.tex#L344-L356)** | No-Go 4: Eliminating Ungraded $G_3$ (Omitting $\Z_3^{\mathrm{PQ}}$) | Axiom 2.2, Prop 2.10, Prop 2.5(1) | Violates Def 2.3(2): conjugate quotient doublets $\Sigma[4],\Sigma[3]$ fuse into one $M_2(\C)$ component locking $\Tr(M^{(1)}M^{(1)\dagger}) = \Tr(M^{(2)}M^{(2)\dagger})$ unless separated by $\Z_3^{\mathrm{PQ}}$ |
| **[Thm 2.15](paper.tex#L358-L374)** | **Four-Factor Design Theorem for $\Gstar$** | Axioms 2.1, 2.2, Prop 2.5(1), Thms 2.7, 2.9, Prop 2.10, Thms 2.11–2.14 | Uniquely selects $\Gstar = 2I \times \Sigma(216 \times 3) \times \Delta(27)_f \times \Z_3^{\mathrm{PQ}}$ (order $6{,}298{,}560$) |
| **[Prop 2.16](paper.tex#L378-L389)** | 3-Generation Discrete Ibáñez–Ross Anomaly Cancellation | Thm 2.15 | Single-generation discrete anomalies $\equiv (1,1,1) \pmod 3$ cancel iff $3 \mid N_{\mathrm{gen}}$ ($N_{\mathrm{gen}} = \dim\mathbf{3}_f = 3$) |
