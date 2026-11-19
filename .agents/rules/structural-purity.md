---
trigger: always_on
---

# Workspace Rules: Mathematical-Physics Structural Purity

When writing or editing papers, documentation, and theoretical derivations in this repository, enforce **strict structural purity and separation**:

1. **Pure Axioms (`\begin{axiom}...\end{axiom}`):**
   - An axiom must state **only the bare postulate or definition** being assumed.
   - Never embed mathematical consequences, known theorems (e.g., Artin–Wedderburn, Frobenius–Schur), proofs, motivations (*"To avoid..."*), or outcome-driven physical jargon (e.g., *"Lorentz spinor, gauge, Higgs, flavor, Yukawa, Connes' Standard Model"*) inside an axiom.
   - Any immediate mathematical consequence of an axiom must be placed *outside* the axiom environment (as an explicit theorem, proposition, or explanatory paragraph).

2. **Strict Statement–Proof Separation (`\begin{theorem}...\end{theorem}` + `\begin{proof}...\end{proof}`):**
   - Every theorem, proposition, or lemma must contain **only the hypotheses and the formal claim** (in italics), with no interwoven derivations (*"Because X, we have Y..."*) or conversational remarks.
   - Every theorem, proposition, or lemma must be followed immediately by an upright `\begin{proof}...\end{proof}` block containing the derivation, kept as concise as possible without belaboring trivial algebraic steps.

3. **Strict Deductive Direction (No Circularity or Outcome-Driven Assumptions):**
   - Never assume or preview physical outcomes (such as Connes' Standard Model algebra $\mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$, $d = 4$ spacetime, chiral spinors, or $N_{\mathrm{gen}} = 3$ generations) inside foundational axioms or early setup paragraphs before they have been derived.
   - Physical identifications must appear strictly as **outcomes** of the mathematical theorems, never as hidden premises.
