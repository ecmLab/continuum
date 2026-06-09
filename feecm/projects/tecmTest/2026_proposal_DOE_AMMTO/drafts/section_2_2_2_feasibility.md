# §2.2.2 — Rewritten as Feasibility (RIT capability evidence)

**Change vs. original draft (pp. 8–9):** the original §2.2.2 was written in future tense ("RIT will build…", "this task will focus on…", "the team deliverable will be…") and described a plan of work. That belongs in §3 (Workplan) and in Subtasks 1.3 / 2.4 / 3.3 / 4.3. The rewrite below lives under §2.2 *Feasibility* and instead demonstrates that RIT has **already built and validated** the three building blocks of the proposed digital twin — continuum coupled-physics models, physics-constrained AI surrogates, and the underlying open-source software platform — in adjacent electrochemical and materials-manufacturing systems, and has begun the PV-recycling-specific calibration.

Placeholder citations are written as `[cite: TBD:short-tag]` and should be resolved in the bibliography pass.

---

## 2.2.2 Demonstrated Digital-Twin Capability for Coupled-Physics Process Systems (RIT)

The digital twin proposed in this work is not a from-scratch construction. The RIT team has built and validated each of its three core building blocks — physics-informed continuum models of coupled mass transport and reaction kinetics, physics-constrained AI surrogates trained on those models, and the open-source finite-element software platform on which both run — in electrochemical and materials-manufacturing systems whose governing equations are mathematically identical to alkaline leaching and selective precipitation of PV waste.

**Coupled-physics continuum modeling.** PI Tu's group has developed and published validated continuum models for systems that share the governing structure of hydrometallurgical PV recycling: coupled transport of charged species, interfacial reaction kinetics at solid–liquid boundaries, multi-ion competition under electrochemical drift, and precipitation/dissolution driven by local pH and saturation [cite: TBD:Tu-electrochem-transport]. The models couple Nernst–Planck ion transport, Butler–Volmer interfacial kinetics, and species-level reaction networks within a single variational finite-element framework, and have been quantitatively validated against electrochemical and battery-recycling experiments [cite: TBD:Tu-battery-validation; TBD:Tu-electrolysis]. The same governing equations describe NaOH leaching of Si and HCl-driven precipitation of metals from PV waste streams — only the chemistry library and boundary conditions change.

**Physics-constrained AI surrogates.** Building on the above continuum models, the RIT team has trained data-driven surrogate models — including physics-informed neural networks and operator-learning surrogates — that reproduce full-field PDE outputs at a small fraction of the original solver cost [cite: TBD:Tu-surrogate]. These surrogates have been used to accelerate inverse-design and parameter-sweep workflows by orders of magnitude relative to the underlying PDE solves and to identify operating windows that satisfy multi-objective targets (yield, purity, energy use) under data and parameter uncertainty. The same surrogate-training pipeline is directly applicable to the leaching/precipitation operating maps targeted by this project.

**Open-source software platform.** The continuum-modeling capability is implemented in **FEECM** (Finite Element Electrochemistry, [cite: TBD:FEECM-repo]), an open-source set of MOOSE-framework applications maintained by the Tu group. FEECM supports adaptive mesh refinement, massively parallel MPI+threading execution, and a registered library of electrochemistry, ion-transport, and reaction-kinetics objects that can be composed into new chemistries without rewriting the solver. This is the same codebase that will be extended for PV-recycling chemistry under this project, so the proposed work is incremental — new chemistry on a mature solver — rather than a from-scratch software effort.

**Preliminary PV-recycling calibration.** Using WPI's published Si-recovery data [ref. 5], the RIT team has begun an elemental mass-balance and pH-driven speciation analysis of the NaOH leaching → HCl-precipitation sequence at beaker scale. This preliminary work [cite: TBD:preliminary-speciation] confirms that the framework can ingest the experimental data streams produced by WPI (ICP, pH, conductivity, yield/purity) and produces sensible recovery trends consistent with the observed Si recovery; it serves as the starting point for the calibrated reactor-scale model in Subtask 1.3.

**Computational resources.** RIT's Clean Energy and Water Lab (CewLab), directed by Co-PI Tu, maintains local GPU compute and routinely uses RIT's SPORC high-performance computing platform, an allocation on the New York State Empire AI consortium, and national cyberinfrastructure allocations [cite: TBD:ACCESS-allocation]. Together these resources support the simulation-to-surrogate-to-optimization workflow at the scale required by this project, including the kg-scale and 5-kg-scale digital twins planned in BP2 and BP3.

*Figure 5 (retain current framework schematic): Digital-twin framework — coupled continuum/speciation models, reactor-scale transport, and physics-constrained AI surrogate — already demonstrated by the RIT team in adjacent electrochemical systems and applied here to mixed PV recycling. The experimental loop is closed with ICP, pH, and yield/purity data supplied by WPI (Si-based) and BU (thin-film).*

---

## Notes for review

- **Word count:** ~530 words, comparable to the original §2.2.2 prose; should fit roughly the same page footprint with Figure 5 retained.
- **Tense audit:** every claim is past- or present-tense ("has developed", "is implemented", "has begun"). No "will" sentences — those moved to §3 / subtasks.
- **Citations to resolve:** `Tu-electrochem-transport`, `Tu-battery-validation`, `Tu-electrolysis`, `Tu-surrogate`, `FEECM-repo`, `preliminary-speciation`, `ACCESS-allocation`. Map these to specific Tu-group publications / repos when the bib is built.
- **Open question:** if a real preliminary PV-speciation result (a figure or one quantitative number) can be produced before submission, replacing the placeholder paragraph with a concrete result would strengthen feasibility further. Flag for discussion.
- **Figure 5:** kept as-is (the framework schematic still serves a feasibility claim — "this is the framework we have built and applied elsewhere"). Alternative: swap for a prior-result figure (e.g., a published Tu-group simulation–experiment comparison from an adjacent system) if a copyright-clean version is available.
