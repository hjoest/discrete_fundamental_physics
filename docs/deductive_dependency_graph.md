# Deductive Dependency Graph: From Axioms 1.1–1.3 and Axiom 2.1 to the Four-Factor Design Theorem (Theorem 2.12) and Companion Papers II–III

This document traces the complete deductive dependency graph of Paper I ([`paper.tex`](paper.tex#L85-L404)) from the four foundational axioms (**Axioms 1.1–1.3** in Section 1 and **Axiom 2.1**, *Optimal Unitary Design Approximation*, in Section 2) through the intermediate propositions and no-go theorems to the **Four-Factor Design Theorem** (**Theorem 2.12**, `\label{thm:design}`) and **Proposition 2.13** (`\label{prop:anomaly_cancellation}`), which in turn provide the algebraic input $\Gstar = 2I \times \Sigma(216 \times 3) \times \Delta(27)_f \times \Z_3^{\mathrm{PQ}}$ for the 4D Poset Collider simulations in Paper II ([`paper2.tex`](paper2.tex)) and the Quaternionic Cartan Soldering / Causal Set Inverse Problem theorems (Propositions 1.2–1.3 and Theorems 2.1, 2.3–2.4) in Paper III ([`paper3.tex`](paper3.tex)).

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
        P16["Prop 1.6: Frobenius-Schur Multiplicativity, Real Collapse & Conjugate-Pair Fusion"]
        P17["Prop 1.7: Unitary Embedding in U(d) & Projective Reduction to SU(d)"]
        P18["Prop 1.8: Resolution of Covering-Relation Triviality Trap"]
        P19["Prop 1.9: Abelianization Kirchhoff Conservation & Perfect-vs.-Solvable Dichotomy"]
    end

    subgraph S2A["Section 2: Optimal Unitary Design Approximation Axiom"]
        A21["Axiom 2.1 (Defs 2.1-2.2): Optimal Unitary Design Approximation"]
    end

    subgraph S21["Section 2.1: Classification of Maximal Unitary t-Designs & Clifford Design Hierarchies"]
        P22["Prop 2.2: Minimal Unitary Degrees d in {1, 2, 3} -> C + H + M_3(C)"]
        P23["Prop 2.3: Minimal Candidate Classification: Z_3; T_7, Delta(27), M_27; ADE in SU(2)"]
        P24["Prop 2.4: Clifford Unitary Design Hierarchy & Irreducible Adjoint Closure in SU(3)"]
        T25["Thm 2.5: Elimination of T_7 & M_27; Unique Solvable Complex Clifford Hierarchy Z_3 < Delta(27) < Sigma(216x3)"]
        P26["Prop 2.6: Equivalence of Perfection [H,H]=H & Maximal Perfect Unitary 5-Design in SU(2)"]
        T27["Thm 2.7: Elimination of Solvable C_n, 2D_2n, 2T, 2O; Selection of Unique Maximal Perfect 5-Design 2I"]
    end

    subgraph S22["Section 2.2: Design-Level Unfolding No-Go Theorems & Four-Factor Design Theorem"]
        T28["Thm 2.8 (No-Go 1): Omitting Sigma(216x3) -> No 2-Design M_8(R) Octet & No Solvable Quaternionic Factor"]
        T29["Thm 2.9 (No-Go 2): Omitting 2I -> Solvable Group, No Perfect 5-Design & 4D Crystal"]
        T210["Thm 2.10 (No-Go 3): Omitting H = Delta(27) -> Quotient-Kernel Collapse 3->2 & Circulant Triviality"]
        T211["Thm 2.11 (No-Go 4): Omitting Z(H) = Z_3 -> Conjugate-Pair Fusion Lock & Single-Abelianization Obstruction"]
        T212["Thm 2.12 (Four-Factor Design Theorem): Unique Group G* = 2I x Sigma(216x3) x Delta(27) x Z_3"]
        P213["Prop 2.13: 3-Generation Discrete Ibanez-Ross Anomaly Cancellation (3 | N_gen)"]
    end

    %% Section 1 edges
    A12 --> P14
    P14 --> P15
    P15 --> P16
    A12 --> P17
    A11 --> P18
    A13 --> P18
    A11 --> P19
    A12 --> P19
    A13 --> P19

    %% Section 2.1 edges
    P15 --> A21
    P16 --> A21
    P17 --> A21
    P19 --> A21

    P14 --> P22
    P15 --> P22
    P17 --> P22
    A21 --> P22

    P22 --> P23
    A21 --> P23

    P14 --> P24
    P15 --> P24
    P23 --> P24

    A21 --> T25
    P23 --> T25
    P24 --> T25

    P15 --> P26
    P19 --> P26
    P23 --> P26

    A21 --> T27
    P19 --> T27
    P26 --> T27

    %% Section 2.2 edges
    A21 --> T28
    P24 --> T28
    T25 --> T28

    A21 --> T29
    P26 --> T29
    T27 --> T29

    A21 --> T210
    P16 --> T210
    T25 --> T210

    A21 --> T211
    P16 --> T211
    P24 --> T211

    A21 --> T212
    P23 --> T212
    T25 --> T212
    T27 --> T212
    P16 --> T212
    T28 --> T212
    T29 --> T212
    T210 --> T212
    T211 --> T212

    T212 --> P213
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
| **[Prop 1.6](paper.tex#L136-L152)** | Frobenius–Schur Multiplicativity, Real Collapse, and Conjugate-Pair Fusion | Prop 1.5 | $\iota(\bigboxtimes \chi_j) = \prod \iota(\chi_j)$, real/abelian-triplet collapse ($\iota = (-1)(-1) = +1$ and $3 \to 2$), and conjugate-pair fusion $\rho_{\bar{\chi}}(x) = \rho_\chi(x)^*$ in $M_{d_\chi}(\C)$ |
| **[Prop 1.7](paper.tex#L154-L161)** | Unitary and Special Unitary Embedding of Finite Group Representations | Axiom 1.2 (Weyl's theorem) | Every representation embeds into $\U(d_\chi)$, and projective tensor conjugation ($d_\chi \ge 2$) factors through $\PU(d_\chi) \cong \SU(d_\chi)/\Z_{d_\chi}$ and embeds at minimal order into $\SU(d_\chi)$ |
| **[Prop 1.8](paper.tex#L163-L170)** | Resolution of Covering-Relation Triviality Trap | Axioms 1.1, 1.3 | Antichain fibers $Z_x$ and covering equivalence $\prec_Z \iff \prec_P$ |
| **[Prop 1.9](paper.tex#L172-L184)** | Abelianization Kirchhoff Conservation Law & Perfect-vs.-Solvable Dichotomy | Axioms 1.1, 1.2, 1.3 | Additive charge conservation in $G^{\mathrm{ab}} = G/[G,G]$, partitioning subgroups into solvable ($K^{\mathrm{ab}} \neq \{1\}$) and perfect ($[K,K]=K \iff K^{\mathrm{ab}} = \{1\}$) |
| **[Axiom 2.1 / Defs 2.1–2.2](paper.tex#L192-L209)** | Optimal Unitary Design Approximation | Props 1.5–1.7, 1.9 & algebraic design optimality postulate | $\Gstar$ is the finite group of minimal order and minimal exponent whose real group algebra $\R[\Gstar]$ simultaneously and faithfully realizes, without real collapse or conjugate-pair fusion, the maximal-degree perfect projective unitary $t$-design and every level of the minimal solvable complex ($\iota = 0$) Clifford unitary design hierarchy |
| **[Prop 2.2](paper.tex#L223-L236)** | Minimal Unitary Degrees of Non-Real Division Algebras | Props 1.4, 1.5, 1.7, Axiom 2.1 | $d_\chi = 1$ ($\U(1) \to \C$), $d_\chi = 2$ ($\SU(2) \to \quat$), $d_\chi = 3$ ($\SU(3) \to M_3(\C)$) |
| **[Prop 2.3](paper.tex#L238-L250)** | Classification of Minimal Finite Group Candidates | Prop 2.2, Axiom 2.1 | Selects $\Z_3 \subset \U(1)$; candidates $\{T_7, \Delta(27), M_{27}\} \subset \SU(3)$; ADE candidates $\{C_n, 2D_{2n}, 2T, 2O, 2I\} \subset \SU(2)$ |
| **[Prop 2.4](paper.tex#L254-L266)** | Clifford Unitary Design Hierarchy & Irreducible Adjoint Closure in $\SU(3)$ | Props 1.4, 1.5, 2.3, Defs 2.1–2.2 | Equivalence: $9$-element $1$-design basis $H/Z(H)$ + $2$-design irreducible $N_{\SU(3)}(H)$ action $\iff |H|=27, \text{exp } 3 \iff H \cong \Delta(27), N_{\SU(3)}(H) \cong \Sigma(216\times 3)$ and irrep classification $\Delta[0..10], \Sigma[0..23]$ |
| **[Thm 2.5](paper.tex#L268-L275)** | Elimination of $T_7$ and $M_{27}$; Selection of $\Delta(27)$ & $\Sigma(216\times 3)$ | Axiom 2.1, Props 2.3, 2.4 | Rules out $T_7$ and $M_{27}$; proves $1 \triangleleft \Z_3 \triangleleft \Delta(27) \triangleleft \Sigma(216 \times 3)$ is the unique minimal solvable complex ($\iota = 0$) Clifford unitary design hierarchy |
| **[Prop 2.6](paper.tex#L279-L291)** | Equivalence of Perfection & Maximal Perfect Unitary $5$-Design in $\SU(2)$ | Props 1.5, 1.9, 2.3, Def 2.1 | Equivalence in $\SU(2)$: $[H,H]=H \iff \mathbf{2}_{-1} \oplus \mathbf{2}_{-1}'$ (unique maximal perfect $5$-design) $\iff H \cong 2I$ with $5$ real ($A_5, z_{\mathrm{spin}}=+1$) and $4$ quaternionic ($z_{\mathrm{spin}}=-1$) irreps ($\sum d_\chi^2 = 60 + 60$) |
| **[Thm 2.7](paper.tex#L293-L300)** | Elimination of Solvable $C_n, 2D_{2n}, 2T, 2O$; Selection of $2I$ | Axiom 2.1, Props 1.9, 2.6 | Rules out solvable ADE subgroups; selects unique maximal-degree perfect unitary $5$-design $2I \subset \SU(2)$ |
| **[Thm 2.8](paper.tex#L306-L313)** | No-Go 1: Eliminating All Subsets Omitting $\Sigma(216 \times 3)$ | Axiom 2.1, Prop 2.4, Thm 2.5 | Violates Axiom 2.1: without $\Sigma(216\times 3) = \Delta(27) \rtimes 2T$, no $2$-design $M_8(\R)$ octet and no solvable quaternionic subgroup independent of $2I$ exist |
| **[Thm 2.9](paper.tex#L315-L322)** | No-Go 2: Eliminating All Subsets Omitting $2I$ | Axiom 2.1, Prop 2.6, Thm 2.7 | Violates Axiom 2.1: without $2I$, $G$ is solvable and all $d_\R \le 4$ irreps have $\Delta(27)$ in their kernel with crystallographic traces |
| **[Thm 2.10](paper.tex#L324-L331)** | No-Go 3: Eliminating All Subsets Omitting $H \cong \Delta(27)$ | Axiom 2.1, Prop 1.6, Thm 2.5 | Violates Axiom 2.1: on $\Delta(27)$-trivial quotient representations, abelian $\Z_3$ triplets collapse $3 \to 2$ simple Wedderburn components and force $\Z_3$-circulant matrices ($|U_1^\dagger U_2| = I_3$) |
| **[Thm 2.11](paper.tex#L333-L345)** | No-Go 4: Eliminating Ungraded $G_3$ (Omitting $Z(H) \cong \Z_3$) | Axiom 2.1, Prop 1.6, Prop 2.4 | Violates Axiom 2.1: conjugate quotient doublets $\Sigma[4],\Sigma[3]$ fuse into one $M_2(\C)$ component locking $\Tr(M^{(1)}M^{(1)\dagger}) = \Tr(M^{(2)}M^{(2)\dagger})$ unless separated by $Z(H) \cong \Z_3$ |
| **[Thm 2.12](paper.tex#L347-L363)** | **Four-Factor Design Theorem for $\Gstar$** | Axiom 2.1, Prop 1.6, Prop 2.3(1), Thms 2.5, 2.7, Thms 2.8–2.11 | Uniquely selects $\Gstar = 2I \times \Sigma(216 \times 3) \times \Delta(27) \times \Z_3$ (order $6{,}298{,}560$) from the $1 + 3 = 4$ design levels |
| **[Prop 2.13](paper.tex#L393-L404)** | 3-Generation Discrete Ibáñez–Ross Anomaly Cancellation | Thm 2.12, Table 1 | Single-generation discrete anomalies $\equiv (1,1,1) \pmod 3$ cancel iff $3 \mid N_{\mathrm{gen}}$ ($N_{\mathrm{gen}} = \dim\mathbf{3}_f = 3$) |
