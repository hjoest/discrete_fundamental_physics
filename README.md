# Discrete Fundamental Physics: Posets, Aperiodic Tilings, and Finite Groups

This repository explores a strictly finitist approach to fundamental physics in which spacetime and matter emerge from:
1. **A locally finite causal poset $(P, \le_P)$** representing discrete spacetime events and causal precedence,
2. **A finite group $G$** representing internal gauge, spin, and flavor degrees of freedom, and
3. **An event poset $Z \subseteq P \times G$** governed by a **local ordered product conservation law**:
   $$\exists\, \sigma \in S_m,\, \tau \in S_n : \quad \prod_{i=1}^m \pi_G(u_{\sigma(i)}) = \prod_{j=1}^n \pi_G(v_{\tau(j)})$$
   balancing the projected group elements of immediate predecessors $u \prec_Z z$ and immediate successors $z \prec_Z v$.

## Repository Structure

- [`docs/idea_outline.md`](docs/idea_outline.md) — Structured synthesis of the core theoretical framework, physical dictionary, and candidate finite groups.
- [`docs/gemini_chat_transcript.md`](docs/gemini_chat_transcript.md) — Complete, clean Markdown transcript of the 18-exchange exploratory dialogue (`Fundamental Physics Theory Based on Aperiodic 4D Spacetime Tiling, Partial Order, and a Finite Group`).
- [`docs/formalization.md`](docs/formalization.md) — Mathematical formalization of the axioms for $(P, G, Z, \prec_Z)$, resolution of the covering-relation triviality trap, real group algebra ($\mathbb{R}[G]$) Wedderburn decomposition, and spectral action observables.
- [`docs/computational_results.md`](docs/computational_results.md) — Full quantitative results from the computational exploration of candidate groups, Real Wedderburn decompositions, Cayley spectra, and the ordered product conservation law.
- [`docs/first_outcome.md`](docs/first_outcome.md) — Synthesis of the first outcome: why $\Sigma(216 \times 3)$ (order 648) emerges as the standout Standard Model candidate, and comparison of Architectural Options A and B.
- [`docs/falsification_report.md`](docs/falsification_report.md) — Execution of the 4-part falsification protocol on Option A ($\Sigma(216 \times 3)$) and Option B ($2I \times \Sigma(216 \times 3)$), ruling out Option A via three exact no-go theorems and evaluating Option B's survival and remaining vulnerabilities.
- [`src/`](src/), [`run_analysis.py`](run_analysis.py), & [`run_falsification.py`](run_falsification.py) — Zero-dependency Python computational toolkit for analyzing candidate finite groups, their Frobenius–Schur / Wedderburn decompositions, Clebsch–Gordan Yukawa tensors, Cayley Laplacian spectra, and combinatorial solutions to the ordered product conservation law on causal posets.
