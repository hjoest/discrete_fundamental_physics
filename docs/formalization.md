# Mathematical Formalization: Causal Posets, Finite Group Fibers, and the Ordered Product Law

## 1. Foundational Primitives and the Two Structural Traps

We seek a strictly finitist formulation of spacetime and internal gauge structure using only:
1. A countable, locally finite partially ordered set $(P, \le_P)$ (the **base causal poset**),
2. A finite group $G$ (the **internal symmetry group**), and
3. An **event poset** $(Z, \le_Z)$ with projections $\pi_P: Z \to P$ and $\pi_G: Z \to G$.

In constructing $(Z, \le_Z)$ from $P$ and $G$, two mathematical traps must be avoided:

### Trap 1: The Antisymmetry / Closed Timelike Curve (CTC) Trap
If one defines $Z = P \times G$ with the naive pullback preorder
$$(p_1, g_1) \le_Z (p_2, g_2) \iff p_1 \le_P p_2,$$
then for any $g_1 \neq g_2 \in G$, both $(p, g_1) \le_Z (p, g_2)$ and $(p, g_2) \le_Z (p, g_1)$ hold. Antisymmetry fails, turning each fiber $F_p = \{p\} \times G$ into a clique of mutually causal events at zero spacetime separation (instantaneous closed timelike loops) and making the transitive reduction (covering relation $\prec_Z$) ill-defined.

### Trap 2: The Complete Bipartite Fiber Trap
Suppose one repairs antisymmetry by setting $Z = P \times G$ with the strict pullback order
$$(p_1, g_1) <_Z (p_2, g_2) \iff p_1 <_P p_2.$$
Then every fiber $F_p = \{p\} \times G$ is an antichain (spacelike), so $(Z, \le_Z)$ is a valid poset. However, its covering relation $\prec_Z$ is then the complete bipartite relation between adjacent fibers:
$$(p_1, g_1) \prec_Z (p_2, g_2) \iff p_1 \prec_P p_2 \quad \forall\, g_1, g_2 \in G.$$
Consequently, for any event $z = (p, g) \in Z$, its immediate predecessors are $\mathrm{Pred}_Z(p, g) = \mathrm{Pred}_P(p) \times G$. An ordered product over $\mathrm{Pred}_Z(p, g)$ is therefore an ordered product of $\lvert\mathrm{Pred}_P(p)\rvert$ complete copies of the entire group $G$.

By a theorem of finite group theory (Rhemtulla; Dénes–Hermann), for any perfect group $G = [G, G]$ (such as $2I \cong \mathrm{SL}(2, \mathbb{F}_5)$ or $\mathrm{PSL}(2, \mathbb{F}_7)$), *every* element of $G$ can be expressed as an ordered product of all elements of $G$. Hence the conservation law $\prod_{u \prec_Z z} \pi_G(u_{\sigma(i)}) = \prod_{z \prec_Z v} \pi_G(v_{\tau(j)})$ would hold trivially at every interior vertex of every poset $P$, carrying zero dynamical information.

---

## 2. Well-Posed Axioms for $(P, G, Z, \prec_Z)$

Both traps are simultaneously resolved by requiring $\pi_P: Z \to P$ to be a **strictly order-preserving projection** while allowing the causal links $\prec_Z$ (or occupied subset $Z \subseteq P \times G$) to be non-complete across fibers.

### Axiom 1 (Base Causal Set)
$(P, <_P)$ is a countable, locally finite strict partially ordered set:
- **Irreflexive & Asymmetric:** $p \not<_P p$, and $p_1 <_P p_2 \implies p_2 \not<_P p_1$.
- **Transitive:** $p_1 <_P p_2 \text{ and } p_2 <_P p_3 \implies p_1 <_P p_3$.
- **Locally Finite:** For all $p_1, p_2 \in P$, the causal interval $[p_1, p_2]_P := \{q \in P \mid p_1 \le_P q \le_P p_2\}$ is finite.
The covering relation (link relation) is $p_1 \prec_P p_2 \iff p_1 <_P p_2 \text{ and } [p_1, p_2]_P = \{p_1, p_2\}$.

### Axiom 2 (Internal Finite Group)
$G$ is a finite group of order $\lvert G \rvert < \infty$ with identity $e$.

### Axiom 3 (Causal Bundle Poset $Z \subseteq P \times G$)
The event space is a subset $Z \subseteq P \times G$ equipped with a locally finite strict partial order $<_Z$ (with covering relation $\prec_Z$) satisfying:
1. **Strict Base Projection (Spacelike Fibers):**
   $$z_1 <_Z z_2 \implies \pi_P(z_1) <_P \pi_P(z_2).$$
   *Immediate consequence:* Every fiber $F_p := Z \cap (\{p\} \times G)$ is an **antichain** in $Z$ (distinct internal states at the same spacetime point $p$ are mutually incomparable / spacelike), and $<_Z$ is automatically antisymmetric and acyclic.
2. **Finite Degree:** For every $z \in Z$, the immediate predecessor set and immediate successor set,
   $$\mathrm{Pred}_Z(z) := \{u \in Z \mid u \prec_Z z\}, \qquad \mathrm{Succ}_Z(z) := \{v \in Z \mid z \prec_Z v\},$$
   are finite.

This axiom encompasses two natural physical realizations:
- **Realization A (Populated Sub-Poset):** $Z \subset P \times G$ is a sparse subset of occupied fiber states with $z_1 <_Z z_2 \iff \pi_P(z_1) <_P \pi_P(z_2)$ restricted to $Z$. When $\lvert F_p \rvert = 1$ for all $p$, $Z$ is the graph of a field section $\phi: P \to G$; when $1 \le \lvert F_p \rvert \le k$, multiple internal multiplet states are co-occupied at $p$.
- **Realization B (Discrete Connection on $P \times G$):** $Z = P \times G$, and covering links $(p_1, g_1) \prec_Z (p_2, g_2)$ exist only when $p_1 \prec_P p_2$ and $g_1^{-1} g_2 \in S_{p_1 p_2} \subset G$, where $S_{p_1 p_2}$ is a local connection / generator step set.

### Axiom 4 (Local Ordered Product Conservation Law)
For any $m$-tuple $\mathbf{g} = (g_1, \dots, g_m) \in G^m$, define its **permutation product set**:
$$\Pi(\mathbf{g}) := \left\{ \prod_{i=1}^m g_{\sigma(i)} \;\middle|\; \sigma \in S_m \right\} \subseteq G.$$
At every interior event $z \in Z$ (where $\mathrm{Pred}_Z(z) = \{u_1, \dots, u_m\} \neq \emptyset$ and $\mathrm{Succ}_Z(z) = \{v_1, \dots, v_n\} \neq \emptyset$), the incoming tuple $\mathbf{g}_{\text{in}}(z) = (\pi_G(u_1), \dots, \pi_G(u_m))$ and outgoing tuple $\mathbf{g}_{\text{out}}(z) = (\pi_G(v_1), \dots, \pi_G(v_n))$ satisfy:
$$\Pi\big(\mathbf{g}_{\text{in}}(z)\big) \cap \Pi\big(\mathbf{g}_{\text{out}}(z)\big) \neq \emptyset,$$
i.e., there exist permutations $\sigma \in S_m$ and $\tau \in S_n$ such that
$$\prod_{i=1}^m \pi_G\big(u_{\sigma(i)}\big) = \prod_{j=1}^n \pi_G\big(v_{\tau(j)}\big).$$

---

## 3. Exact Algebraic Properties of the Conservation Law

### Proposition 1 (Abelianization & Exact Kirchhoff Current Conservation)
Let $[G, G]$ be the commutator subgroup of $G$, and let $\alpha: G \to G^{\text{ab}} := G/[G, G]$ be the canonical abelianization homomorphism. For any $\mathbf{g} = (g_1, \dots, g_m) \in G^m$:
$$\Pi(\mathbf{g}) \subseteq g_1 g_2 \cdots g_m [G, G].$$
Consequently, $\alpha(\Pi(\mathbf{g}))$ is a singleton, and Axiom 4 implies the ordering-independent abelian conservation law at every interior event $z \in Z$:
$$\sum_{u \prec_Z z} \alpha\big(\pi_G(u)\big) = \sum_{z \prec_Z v} \alpha\big(\pi_G(v)\big) \quad \text{in } G^{\text{ab}}.$$

*Proof:* In $G/[G, G]$ the group operation is commutative, so $\alpha(g_{\sigma(1)} \cdots g_{\sigma(m)}) = \sum_{i=1}^m \alpha(g_i)$ for every $\sigma \in S_m$. $\square$

### Proposition 2 (Automorphism Covariance & Time-Reversal Duality)
For any tuple $\mathbf{g} = (g_1, \dots, g_m) \in G^m$:
1. **Automorphism / Global Gauge Covariance:** For every $\psi \in \mathrm{Aut}(G)$ (and in particular every inner automorphism $\psi_h(g) = h g h^{-1}$),
   $$\Pi\big(\psi(g_1), \dots, \psi(g_m)\big) = \psi\big(\Pi(g_1, \dots, g_m)\big).$$
   Thus $\mathrm{Aut}(G)$ acts as a global symmetry group on the space of valid configurations $\pi_G: Z \to G$.
2. **Inversion / Order-Reversal Duality:**
   $$\Pi\big(g_1^{-1}, \dots, g_m^{-1}\big) = \big(\Pi(g_1, \dots, g_m)\big)^{-1}.$$
   Reversing the causal order of $Z$ ($<_Z \mapsto >_Z$) while inverting group elements ($g \mapsto g^{-1}$) preserves all solutions of Axiom 4 (discrete $CPT$ / time-reversal invariance).

### Proposition 3 (Closed-Form $(1 \to 2)$ Branching Probability and Commuting Probability)
Let $k(G)$ denote the number of conjugacy classes of $G$, and let $\mathrm{cp}(G) := \lvert\{(g_1, g_2) \in G^2 \mid g_1 g_2 = g_2 g_1\}\rvert / \lvert G \rvert^2 = k(G) / \lvert G \rvert$ be the commuting probability of $G$.
Under uniform random sampling of group labels in $G$, the exact fraction of $(1 \to n)$ transitions ($g_1 \to (h_1, \dots, h_n)$) satisfying Axiom 4 is:
$$p_{1, n}(G) = \frac{\mathbb{E}_{\mathbf{h} \in G^n}\big[\lvert\Pi(\mathbf{h})\rvert\big]}{\lvert G \rvert}.$$
In particular, for elementary binary branching $(1 \to 2)$ or binary merging $(2 \to 1)$,
$$\mathbb{E}_{(h_1, h_2) \in G^2}\big[\lvert\Pi(h_1, h_2)\rvert\big] = 2 - \mathrm{cp}(G) = 2 - \frac{k(G)}{\lvert G \rvert},$$
which yields the exact identity:
$$p_{1, 2}(G) = p_{2, 1}(G) = \frac{2\lvert G \rvert - k(G)}{\lvert G \rvert^2}.$$

*Proof:* For a pair $(h_1, h_2) \in G^2$, $\Pi(h_1, h_2) = \{h_1 h_2, h_2 h_1\}$ has cardinality $1$ if $h_1$ and $h_2$ commute (which occurs for $k(G)\lvert G \rvert$ pairs by Burnside's lemma) and cardinality $2$ otherwise. For each fixed $(h_1, h_2)$, there are $\lvert\Pi(h_1, h_2)\rvert$ values of $g_1 \in G$ such that $g_1 \in \Pi(h_1, h_2)$. Dividing by $\lvert G \rvert^3$ gives the result. $\square$

---

## 4. Real Group Algebra $\mathbb{R}[G]$ and Connes' Standard Model Algebra $\mathcal{A}_F$

In the Connes–Chamseddine Spectral Standard Model, the internal algebra is the **real** finite-dimensional $^*$-algebra
$$\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C}).$$
Notice that over $\mathbb{C}$, the quaternions $\mathbb{H}$ cannot be distinguished from $M_2(\mathbb{C})$ because $\mathbb{H} \otimes_{\mathbb{R}} \mathbb{C} \cong M_2(\mathbb{C})$. Therefore, to derive $\mathcal{A}_F$ from a finite group $G$, one must analyze the **real group algebra $\mathbb{R}[G]$** rather than $\mathbb{C}[G]$.

### Theorem (Real Wedderburn–Artin Decomposition via Frobenius–Schur Indicators)
Let $\mathrm{Irr}_{\mathbb{C}}(G)$ be the set of complex irreducible characters of $G$. For each $\chi \in \mathrm{Irr}_{\mathbb{C}}(G)$ of degree $d_\chi = \chi(e)$, the **Frobenius–Schur indicator** is
$$\iota(\chi) := \frac{1}{\lvert G \rvert} \sum_{g \in G} \chi(g^2) \in \{+1, 0, -1\}.$$
The real group algebra $\mathbb{R}[G]$ decomposes into simple matrix algebras over the three real division algebras $\mathbb{R}, \mathbb{C}, \mathbb{H}$ as:
$$\mathbb{R}[G] \;\cong\; \bigoplus_{\substack{\chi \in \mathrm{Irr}_{\mathbb{C}}(G) \\ \iota(\chi) = +1}} M_{d_\chi}(\mathbb{R}) \;\oplus\; \bigoplus_{\substack{\{\chi, \bar{\chi}\} \\ \iota(\chi) = 0}} M_{d_\chi}(\mathbb{C}) \;\oplus\; \bigoplus_{\substack{\chi \in \mathrm{Irr}_{\mathbb{C}}(G) \\ \iota(\chi) = -1}} M_{d_\chi / 2}(\mathbb{H}).$$

This yields an exact algebraic dictionary for recovering $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$:
1. **Quaternionic Electroweak Sector ($\mathbb{H}$):** Requires at least one 2-dimensional complex irrep with **$\iota(\chi) = -1$** (pseudoreal / quaternionic spinor), contributing $M_{2/2}(\mathbb{H}) = \mathbb{H}$.
2. **Complex Color Sector ($M_3(\mathbb{C})$):** Requires at least one conjugate pair $\{\chi, \bar{\chi}\}$ of 3-dimensional complex irreps with **$\iota(\chi) = 0$**, contributing $M_3(\mathbb{C})$.
3. **Abelian / Scalar Sector ($\mathbb{C}$ or $\mathbb{R} \oplus \mathbb{R}$):** Requires a 1-dimensional complex conjugate pair ($d_\chi = 1, \iota(\chi) = 0 \implies \mathbb{C}$) or real/complex subalgebra projection.
4. **Generation Multiplicity ($N_g$):** The number of distinct copies of the fundamental $(\mathbb{H}, M_3(\mathbb{C}))$ blocks (or the dimension of the flavor multiplicity space) determines the number of fermion generations $N_g$.

---

## 5. Discrete Spectral Action and Bare Physical Observables

Let $S = S^{-1} \subset G \setminus \{e\}$ be a symmetric generating set of $G$ of degree $k_S = \lvert S \rvert$. The adjacency operator $A_S$ and discrete Laplacian $\Delta_S$ on the Cayley graph $\mathrm{Cay}(G, S)$ act on $\psi \in \ell^2(G)$ by:
$$(A_S \psi)(g) = \sum_{s \in S} \psi(g s), \qquad \Delta_S = k_S I - A_S.$$

### 1. Class-Symmetric Generators and Character Eigenvalues
When $S$ is a union of conjugacy classes of $G$ ($g S g^{-1} = S$), $A_S$ lies in the center $Z(\mathbb{R}[G])$ of the group algebra. By Schur's lemma, on each irreducible representation space $V_\chi$ (which appears with multiplicity $d_\chi$ in the regular representation $\ell^2(G)$, contributing $d_\chi^2$ states), $A_S$ acts as the scalar eigenvalue:
$$\lambda_\chi(A_S) = \frac{1}{d_\chi} \sum_{s \in S} \chi(s), \qquad \lambda_\chi(\Delta_S) = k_S - \frac{1}{d_\chi} \sum_{s \in S} \chi(s).$$

### 2. Spectral Quartic Ratio $\kappa(G, S)$ and Closed 4-Walks (Plaquettes)
In the Connes–Chamseddine spectral action, the bare ratio of the quartic Higgs self-coupling $\lambda(\Lambda)$ to the squared Yukawa trace at the cutoff scale $\Lambda$ is controlled by the dimensionless spectral moment ratio of the internal Dirac/Yukawa operator $D_G$:
$$\kappa(D_G) := \lvert G \rvert \frac{\mathrm{Tr}(D_G^4)}{\big(\mathrm{Tr}(D_G^2)\big)^2} = \frac{\frac{1}{\lvert G \rvert} \sum_{i=1}^{\lvert G \rvert} \mu_i^4}{\left(\frac{1}{\lvert G \rvert} \sum_{i=1}^{\lvert G \rvert} \mu_i^2\right)^2},$$
where $\{\mu_i\}_{i=1}^{\lvert G \rvert}$ are the eigenvalues of $D_G$.
- By the Cauchy–Schwarz inequality, $\kappa(D_G) \ge 1$, with equality if and only if all non-zero eigenvalues have identical magnitude (flat spectrum) and no zero modes exist.
- For $D_G = A_S$, $\mathrm{Tr}(A_S^2) = \lvert G \rvert k_S$ and $\mathrm{Tr}(A_S^4) = \lvert G \rvert W_4(S)$, where $W_4(S) = \lvert\{(s_1, s_2, s_3, s_4) \in S^4 \mid s_1 s_2 s_3 s_4 = e\}\rvert$ is the number of 4-step closed walks (plaquettes) from $e$. Thus:
  $$\kappa(A_S) = \frac{W_4(S)}{k_S^2}.$$
- Under the standard spectral action normalization (Chamseddine–Connes), the bare quartic coupling at the unification scale is $\lambda(\Lambda) = \frac{\pi^2}{2 f_0} \kappa(D_G) = g^2(\Lambda) \kappa(D_G)$, giving a tree-level bare Higgs-to-$W$ mass ratio:
  $$\frac{m_H(\Lambda)}{m_W(\Lambda)} = 2 \sqrt{\kappa(D_G)} \quad \text{(before RG running from $\Lambda$ to $M_Z$)}.$$

### 3. Gauge Coupling Normalization and the Bare Weinberg Angle
When the weak $SU(2)_L$ and strong $SU(3)_c$ sectors arise from distinct factor groups $G_{\text{weak}} \times G_{\text{color}}$ or distinct Wedderburn blocks $\mathcal{B}_2 \cong M_{n_2}(\mathbb{H})$ and $\mathcal{B}_3 \cong M_{n_3}(\mathbb{C})$ inside $\mathbb{R}[G]$, the canonical trace on the regular representation $\ell^2(G)$ weights each block by its Plancherel measure $d_\chi^2 / \lvert G \rvert$. This provides a rigid, parameter-free prediction for the bare coupling ratios at the unification scale $\Lambda$.
