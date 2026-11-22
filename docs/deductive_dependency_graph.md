# Deductive Dependency Graph: From Axioms 1.1–1.3 and Axiom 2.1 to the Four-Factor Design Theorem (Theorem 2.14) and Companion Papers II–III

This document traces the complete deductive dependency graph of Paper I ([`paper.tex`](paper.tex#L85-L377)) from the four foundational axioms (**Axioms 1.1–1.3** in Section 1 and **Axiom 2.1**, *Optimal Unitary Design Approximation*, in Section 2) through the intermediate propositions and no-go theorems to the **Four-Factor Design Theorem** (**Theorem 2.14**, `\label{thm:design}`) and **Proposition 2.15** (`\label{prop:anomaly_cancellation}`), which in turn provide the algebraic input $\Gstar = 2I \times \Sigma(216 \times 3) \times \Delta(27)_f \times \Z_3^{\mathrm{PQ}}$ for the 4D Poset Collider simulations in Paper II ([`paper2.tex`](paper2.tex)) and the Quaternionic Cartan Soldering / Causal Set Inverse Problem theorems (Propositions 1.2–1.3 and Theorems 2.1, 2.3–2.4) in Paper III ([`paper3.tex`](paper3.tex)).

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

    subgraph S2A["Section 2: Optimal Unitary Design Approximation Axiom"]
        A21["Axiom 2.1 (Defs 2.1-2.2): Optimal Unitary Design Approximation"]
    end

    subgraph S21["Section 2.1: Classification of Maximal Unitary t-Designs & Clifford Design Hierarchies"]
        P22["Prop 2.2: Unitary Embedding H in U(d), SU(d)"]
        P23["Prop 2.3: Minimal Unitary Degrees d in {1, 2, 3} -> C + H + M_3(C)"]
        P24["Prop 2.4: Minimal Candidate Classification: Z_3; T_7, Delta(27), M_27; ADE in SU(2)"]
        P25["Prop 2.5: M_3(C) Projective Unitary 1-Design Basis & 2-Design M_8(R) Normalizer Closure"]
        T26["Thm 2.6: Elimination of T_7 & M_27; Unique Solvable Complex Clifford Hierarchy Z_3 < Delta(27) < Sigma(216x3)"]
        P27["Prop 2.7: Equivalence of Perfection [H,H]=H & Maximal Perfect Unitary 5-Design in SU(2)"]
        T28["Thm 2.8: Elimination of Solvable C_n, 2D_2n, 2T, 2O; Selection of Unique Maximal Perfect 5-Design 2I"]
    end

    subgraph S22["Section 2.2: Design-Level Unfolding No-Go Theorems & Four-Factor Design Theorem"]
        P29["Prop 2.9: Frobenius-Schur Multiplicativity & Conjugate-Component Reality Constraint"]
        T210["Thm 2.10 (No-Go 1): Omitting Sigma(216x3) -> No 2-Design M_8(R) Octet & No Solvable Quaternionic Factor"]
        T211["Thm 2.11 (No-Go 2): Omitting 2I -> Solvable Group, No Perfect 5-Design & 4D Crystal"]
        T212["Thm 2.12 (No-Go 3): Omitting Delta(27)_f -> Quotient-Kernel Collapse 3->2 & Circulant Triviality"]
        T213["Thm 2.13 (No-Go 4): Omitting Z_3^PQ -> Conjugate-Component Reality Lock & Single-Abelianization Obstruction"]
        T214["Thm 2.14 (Four-Factor Design Theorem): Unique Group G* = 2I x Sigma(216x3) x Delta(27)_f x Z_3^PQ"]
        P215["Prop 2.15: 3-Generation Discrete Ibanez-Ross Anomaly Cancellation (3 | N_gen)"]
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
    A12 --> P22
    P14 --> P23
    P15 --> P23
    P17 --> P23
    P22 --> P23
    A21 --> P23

    P23 --> P24
    A21 --> P24

    P14 --> P25
    P15 --> P25
    P24 --> P25

    A21 --> T26
    P24 --> T26
    P25 --> T26

    P15 --> P27
    P17 --> P27
    P24 --> P27

    A21 --> T28
    P17 --> T28
    P27 --> T28

    %% Section 2.2 edges
    P15 --> P29

    A21 --> T210
    P25 --> T210
    T26 --> T210

    A21 --> T211
    P27 --> T211
    T28 --> T211

    A21 --> T212
    P29 --> T212
    T26 --> T212

    A21 --> T213
    P29 --> T213
    P24 --> T213

    A21 --> T214
    P24 --> T214
    T26 --> T214
    T28 --> T214
    P29 --> T214
    T210 --> T214
    T211 --> T214
    T212 --> T214
    T213 --> T214

    T214 --> P215
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
| **[Axiom 2.1 / Defs 2.1–2.2](paper.tex#L165-L182)** | Optimal Unitary Design Approximation | Algebraic design optimality postulate | $\Gstar$ is the finite group of minimal order and minimal exponent whose real group algebra $\R[\Gstar]$ simultaneously and faithfully realizes, without real or conjugate-pair collapse, every level of the maximal-degree perfect and solvable complex ($\iota = 0$) projective unitary $t$-design hierarchies |
| **[Prop 2.2](paper.tex#L196-L203)** | Unitary Embedding of Finite Group Representations | Axiom 1.2 (Weyl's theorem) | Every finite subgroup embeds into $\U(d_\chi)$ or $\SU(d_\chi)$ |
| **[Prop 2.3](paper.tex#L205-L218)** | Minimal Unitary Degrees of Non-Real Division Algebras | Props 1.4, 1.5, 1.7, 2.2, Axiom 2.1 | $d_\chi = 1$ ($\U(1) \to \C$), $d_\chi = 2$ ($\SU(2) \to \quat$), $d_\chi = 3$ ($\SU(3) \to M_3(\C)$) |
| **[Prop 2.4](paper.tex#L220-L233)** | Classification of Minimal Finite Group Candidates | Prop 2.3, Axiom 2.1 | Selects $\Z_3 \subset \U(1)$; candidates $\{T_7, \Delta(27), M_{27}\} \subset \SU(3)$; ADE candidates $\{C_n, 2D_{2n}, 2T, 2O, 2I\} \subset \SU(2)$ |
| **[Prop 2.5](paper.tex#L237-L249)** | $M_3(\C)$ Projective Unitary $1$-Design Basis & $2$-Design $M_8(\R)$ Normalizer Closure | Props 1.4, 1.5, 2.4 | Equivalence: $9$-element $1$-design basis $H/Z(H)$ + $2$-design transitive $N_{\SU(3)}(H)$ action $\iff |H|=27, \text{exp } 3 \iff H \cong \Delta(27), N_{\SU(3)}(H) \cong \Sigma(216\times 3)$ |
| **[Thm 2.6](paper.tex#L251-L258)** | Elimination of $T_7$ and $M_{27}$; Selection of $\Delta(27)$ & $\Sigma(216\times 3)$ | Axiom 2.1, Props 2.4, 2.5 | Rules out $T_7$ and $M_{27}$; proves $1 \triangleleft \Z_3 \triangleleft \Delta(27) \triangleleft \Sigma(216 \times 3)$ is the unique minimal solvable complex ($\iota = 0$) Clifford unitary design hierarchy |
| **[Prop 2.7](paper.tex#L262-L274)** | Equivalence of Perfection & Maximal Perfect Unitary $5$-Design in $\SU(2)$ | Props 1.5, 1.7, 2.4 | Equivalence in $\SU(2)$: $[H,H]=H \iff \mathbf{2}_{-1} \oplus \mathbf{2}_{-1}'$ (unique maximal perfect $5$-design) $\iff H \cong 2I$ |
| **[Thm 2.8](paper.tex#L276-L283)** | Elimination of Solvable $C_n, 2D_{2n}, 2T, 2O$; Selection of $2I$ | Axiom 2.1, Props 1.7, 2.7 | Rules out solvable ADE subgroups; selects unique maximal-degree perfect unitary $5$-design $2I \subset \SU(2)$ |
| **[Prop 2.9](paper.tex#L289-L301)** | Frobenius–Schur Multiplicativity & Conjugate-Component Reality | Prop 1.5 | $\iota(\bigboxtimes \chi_j) = \prod \iota(\chi_j)$ and Hermitian conjugation $\rho_{\bar{\chi}}(x) = \rho_\chi(x)^*$ on shared complex components $M_{d_\chi}(\C)$ |
| **[Thm 2.10](paper.tex#L305-L312)** | No-Go 1: Eliminating All Subsets Omitting $\Sigma(216 \times 3)$ | Axiom 2.1, Prop 2.5, Thm 2.6 | Violates Axiom 2.1: without $\Sigma(216\times 3) = \Delta(27) \rtimes 2T$, no $2$-design $M_8(\R)$ octet and no solvable quaternionic subgroup independent of $2I$ exist |
| **[Thm 2.11](paper.tex#L314-L321)** | No-Go 2: Eliminating All Subsets Omitting $2I$ | Axiom 2.1, Prop 2.7, Thm 2.8 | Violates Axiom 2.1: without $2I$, $G$ is solvable and all $d_\R \le 4$ irreps have $\Delta(27)$ in their kernel with crystallographic traces |
| **[Thm 2.12](paper.tex#L323-L330)** | No-Go 3: Eliminating All Subsets Omitting $\Delta(27)_f$ | Axiom 2.1, Prop 2.9, Thm 2.6 | Violates Axiom 2.1: on $\Delta(27)$-trivial quotient representations, abelian $\Z_3$ triplets collapse $3 \to 2$ simple Wedderburn components and force $\Z_3$-circulant matrices ($|U_1^\dagger U_2| = I_3$) |
| **[Thm 2.13](paper.tex#L332-L344)** | No-Go 4: Eliminating Ungraded $G_3$ (Omitting $\Z_3^{\mathrm{PQ}}$) | Axiom 2.1, Prop 2.9, Prop 2.4(1) | Violates Axiom 2.1: conjugate quotient doublets $\Sigma[4],\Sigma[3]$ fuse into one $M_2(\C)$ component locking $\Tr(M^{(1)}M^{(1)\dagger}) = \Tr(M^{(2)}M^{(2)\dagger})$ unless separated by $\Z_3^{\mathrm{PQ}}$ |
| **[Thm 2.14](paper.tex#L346-L362)** | **Four-Factor Design Theorem for $\Gstar$** | Axiom 2.1, Prop 2.4(1), Thms 2.6, 2.8, Prop 2.9, Thms 2.10–2.13 | Uniquely selects $\Gstar = 2I \times \Sigma(216 \times 3) \times \Delta(27)_f \times \Z_3^{\mathrm{PQ}}$ (order $6{,}298{,}560$) from the $1 + 3 = 4$ design levels |
| **[Prop 2.15](paper.tex#L366-L377)** | 3-Generation Discrete Ibáñez–Ross Anomaly Cancellation | Thm 2.14 | Single-generation discrete anomalies $\equiv (1,1,1) \pmod 3$ cancel iff $3 \mid N_{\mathrm{gen}}$ ($N_{\mathrm{gen}} = \dim\mathbf{3}_f = 3$) |
