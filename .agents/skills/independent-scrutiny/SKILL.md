---
name: independent-scrutiny
description: >-
  Runs an independent mathematical and physical scrutiny of the repository using
  Claude Opus 5.5 (claude-opus-5-5-high), writes or reduces the findings in
  docs/independent_scrutiny_report.md without committing, returns control to the
  paper-writing model, and proposes a plan to the user for addressing the remaining
  findings. Use when the user asks to scrutinize, audit, or review the paper/repository
  with Opus or run the independent scrutiny workflow.
---

# Independent Scrutiny Workflow (Claude Opus 5.5 -> Paper-Writing Model)

## Overview

This skill runs a two-stage cross-model scrutiny and proposal workflow using the `swarm` tool:

1. **Stage 1 (Dispatch to Claude Opus 5.5 Scrutinizer):** The paper-writing model (`coordinator`) dispatches the scrutiny task to a `scrutinizer` worker pinned to `claude-opus-5-5-high` ("Claude Opus 5.5 (High)") with `workspace:"inherit"`.
2. **Stage 2 (Scrutiny on Claude Opus 5.5):** Running as Claude Opus 5.5, `scrutinizer` reads any preexisting `docs/independent_scrutiny_report.md` and recent commits, removes topics that have been successfully resolved or openly acknowledged in the paper, audits the paper/code/reports for any remaining or newly introduced issues, writes the reduced report to `docs/independent_scrutiny_report.md` **without committing**, and posts `@coordinator [done] <summary>` to wake the paper-writing model.
3. **Stage 3 (Proposal on the Paper-Writing Model):** Woken by `@coordinator [done]`, the paper-writing model reads `docs/independent_scrutiny_report.md` and presents a structured proposal to the user on how to handle each remaining point before making edits.

---

## Workflow

### Stage 1: Dispatch Scrutiny to Claude Opus 5.5 (Executed on the Paper-Writing Model)

1. **Check swarm status:**
   - Call `swarm{action:"status"}`.
2. **If no swarm exists yet (first run in this conversation):**
   - Call `swarm{action:"title", title:"Independent Scrutiny"}`.
   - Call `swarm{action:"add", name:"scrutinizer", model:"claude-opus-5-5-high", workspace:"inherit", brief:"..."}` with the Stage 2 instructions below.
   - If `swarm_monitor` is in your available skills list, also add the monitor (`swarm{action:"add", name:"monitor", brief:"You are this swarm's stall monitor. Load the swarm_monitor skill and follow it. Post only to break a stall."}`).
   - Include the `open_board` value from the `add`/`status` result on its own line at the end of your reply and **end your turn** so `scrutinizer` runs on Claude Opus 5.5.
3. **If `scrutinizer` already exists in the swarm (`waiting` or `stopped`):**
   - Wake it with `post_message{text:"@scrutinizer Re-audit the repository against the latest commits, update/reduce docs/independent_scrutiny_report.md without committing, and post @coordinator [done] with a summary."}` and **end your turn**.

---

### Stage 2: Perform Scrutiny & Write/Reduce Report (Executed by `scrutinizer` on Claude Opus 5.5)

1. **Read preexisting findings and recent diffs:**
   - Read `docs/independent_scrutiny_report.md` if it exists (or `docs/honest_summary.md` if migrating).
   - Run `git log --oneline -n 5` and inspect recent commits (`git show --stat HEAD`, `git diff HEAD~1 HEAD`) to see what changes the paper-writing model recently applied.
2. **Audit and reduce preexisting topics:**
   - Verify each preexisting critique against `docs/paper.tex`, `run_falsification.py`, `run_scrutiny.py`, `src/*.py`, and `docs/*.md`.
   - **Remove** any topic that has been successfully fixed or is now openly and accurately stated in the paper (e.g., in Section 8 *Open Problems* or Section 9 *Critiques*) without contradiction elsewhere in the paper or generated reports.
3. **Scrutinize for remaining or newly introduced issues:**
   - Check mathematical claims, group-theoretic theorems, renormalization-group factors, and numerical predictions in `docs/paper.tex` against what `src/falsification_tests.py`, `run_falsification.py`, and `run_scrutiny.py` actually compute.
   - Check that generated reports (`docs/*phase*.md`, `docs/scrutiny_audit.md`, `docs/viable_group_G.md`, `docs/falsification_report.md`) and summary files (`README.md`, `docs/final_synthesis.md`, `docs/assessment.md`, `docs/design_theorem.md`) do not contradict the revised paper.
4. **Write `docs/independent_scrutiny_report.md`:**
   - Write (or overwrite) `docs/independent_scrutiny_report.md` containing **only** the issues that remain incorrect, unsupported, or internally contradictory (if all issues are resolved, state clearly that no correctness or honesty issues remain).
   - Include the commit hash audited (`git rev-parse --short HEAD`), the date, a **Scope** paragraph listing the openly acknowledged limitations that are excluded, and how to reproduce any numerical checks.
5. **Do NOT commit anything:**
   - Verify with `git status -s -uall` that `docs/independent_scrutiny_report.md` is left untracked/uncommitted and no git commit has been created.
6. **Notify `@coordinator`:**
   - Call `post_message{text:"@coordinator [done] <summary of removed and remaining items in docs/independent_scrutiny_report.md>"}` and end your turn.

---

### Stage 3: Read Report & Propose Resolution Plan (Executed on the Paper-Writing Model)

When woken by `@coordinator [done]` from `scrutinizer`:

1. **Read the board and `docs/independent_scrutiny_report.md`:**
   - Call `read_board` and `read_file` on `/home/hjoest/Development/src/bitbucket.org/hjoest/physics/docs/independent_scrutiny_report.md`.
   - If a `monitor` worker is active, post `[done] Independent scrutiny complete and proposal ready` (with no `@`) on the board so the monitor goes quiet.
2. **Analyze each finding:**
   - Inspect the exact lines referenced in `docs/paper.tex`, `run_falsification.py`, `run_scrutiny.py`, `src/*.py`, or `docs/*.md`.
3. **Present a structured proposal to the user:**
   - Summarize what Claude Opus 5.5 removed as resolved and what remains in `docs/independent_scrutiny_report.md`.
   - For each remaining point, propose a concrete resolution (e.g., code/report fix, theorem refinement in Section 9, or open-problem addition in Section 8) and wait for the user's approval before editing or committing.
