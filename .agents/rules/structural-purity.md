---
trigger: always_on
---

# Workspace Rules: Mathematical-Physics Structural Purity

When writing or editing papers, documentation, and theoretical derivations in this repository, enforce **strict structural purity and separation**:

1. **Pure Axioms (`\begin{axiom}...\end{axiom}`):**
   - An axiom must state **only the bare postulate or definition** being assumed.
   - Never embed mathematical consequences, known theorems (e.g., Artin–Wedderburn, Frobenius–Schur), proofs, motivations (*"To avoid..."*), or outcome-driven physical jargon (e.g., *"Lorentz spinor, gauge, Higgs, flavor, Yukawa, Connes' Standard Model"*) inside an axiom.
   - Any immediate mathematical consequence of an axiom must be placed *outside* the axiom environment (as an explicit theorem, proposition, or explanatory paragraph).

2. **Strict Statement–Proof Separation in Mathematical Sections (`\begin{theorem}...\end{theorem}` + `\begin{proof}...\end{proof}` in Sections 1–2):**
   - In the mathematical foundation sections (Sections 1 and 2), every theorem, proposition, or lemma must contain **only the hypotheses and the formal claim** (in italics), with no interwoven derivations (*"Because X, we have Y..."*) or conversational remarks, and must be followed immediately by an upright `\begin{proof}...\end{proof}` block containing the derivation, kept as concise as possible without belaboring trivial algebraic steps.
   - From Section 3 onward (where the focus shifts from pure mathematics to physics), use an integrated physics exposition style without separate `\begin{proof}...\end{proof}` blocks.

3. **Strict Deductive Direction (No Circularity or Outcome-Driven Assumptions):**
   - Never assume or preview physical outcomes (such as Connes' Standard Model algebra $\mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$, $d = 4$ spacetime, chiral spinors, or $N_{\mathrm{gen}} = 3$ generations) inside foundational axioms or early setup paragraphs before they have been derived.
   - Physical identifications must appear strictly as **outcomes** of the mathematical theorems, never as hidden premises.

4. **Clear Separation of Classical Theorems (No Theorem Mishmash):**
   - Never lump multiple distinct theorems into a single conflated clause (e.g., *"By Maschke's, Artin–Wedderburn's, and Frobenius's theorems..."*) or jam conclusions from different theorems into one unexplained equation.
   - State each theorem, its specific hypothesis, and its exact mathematical conclusion in clear, sequential steps.
