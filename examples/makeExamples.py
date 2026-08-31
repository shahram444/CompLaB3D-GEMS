#!/usr/bin/env python3
"""
Build the CompLB3D example suite.

Fifteen tiny cases, one per capability, each in its own folder with everything
it needs: CompLaB.xml, the geometry, the kinetics headers it compiles against,
any metabolic model file, and a README saying what to expect.

Why a generator rather than fifteen hand-written files: the cases share a
domain, a set of material numbers and a house style, and a correction to a tag
name has to land in all of them at once. Regenerating is also how the suite
stays honest after validateExamples.py finds something.

    python3 makeGeometry.py      # first, writes input_shared/
    python3 makeExamples.py      # then this
    python3 validateExamples.py  # then check

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
"""

import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "_build"))
import blocks as K                      # noqa: E402
import kinetics as KIN                  # noqa: E402

ROOT = os.path.dirname(os.path.abspath(__file__))
SHARED = os.path.join(ROOT, "input_shared")


# ===========================================================================
# The toy metabolic model, used by examples 09, 10 and 12.
#
# Small enough to solve by hand, which is the whole point: the expected growth
# rate is a formula, not a number this code produced.
#
#   metabolites (rows)   S_c   O_c   P_c
#   reactions (columns)  0 EX_S   "S_c ->"    negative flux imports S
#                        1 EX_O   "O_c ->"    negative flux imports O
#                        2 EX_P   "P_c ->"    positive flux exports P
#                        3 BIO    S_c + 0.5 O_c -> P_c, the objective
#
# Steady state gives v0 = -v3 and v1 = -0.5 v3, so with uptake bounds Vs and Vo
#
#       growth  =  min( Vs,  2 * Vo )
#
# and CompLaB sets each bound by Michaelis-Menten, V = Vmax * C / (Kc + C).
# ===========================================================================
TOY_MODEL = """<?xml version="1.0" ?>
<!-- Toy metabolic model for the CompLB3D examples.

     Three metabolites, four reactions, solvable by hand:

       rows      S_c   O_c   P_c
       columns   0 EX_S    S_c ->            negative flux = import
                 1 EX_O    O_c ->            negative flux = import
                 2 EX_P    P_c ->            positive flux = export
                 3 BIO     S_c + 0.5 O_c -> P_c        (the objective)

     Steady state forces v0 = -v3 and v1 = -0.5 v3, so with uptake bounds
     Vs and Vo the optimum is

           growth = min( Vs, 2 * Vo )

     Generated in the same format extractMM.py writes: S is flattened
     row major, metabolite major, so S[i][j] is entry i * nrxn + j.
-->
<Metabolic_Model>
    <nmet>3</nmet>
    <nrxn>4</nrxn>
    <objLoc>3</objLoc>

    <!-- S, 3 rows of 4 -->
    <S>
        -1  0  0  -1
         0 -1  0  -0.5
         0  0 -1   1
    </S>

    <b>0 0 0</b>
    <c>0 0 0 1</c>

    <lb>-1000 -1000     0    0</lb>
    <ub> 1000  1000  1000 1000</ub>
</Metabolic_Model>
"""


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)




RATES_SYM = '''# A rate law CompLB3D reads at run time, instead of one compiled into
# defineKinetics.hh.  Iron reduction by Geobacter:
#
#     acetate + 4.8 Fe(III)  =  4.8 Fe(II) + 1.845 HCO3-
#
# The growth line is the sort of expression tools/fit_symbolic.py returns from a
# sweep of a metabolic model.  The two substrate lines are NOT fitted: they are
# that reaction, written out, as multiples of growth.  Deriving them keeps the
# 1 : 4.8 ratio exact, which four independent fits would not.

provenance example 16, iron reduction by Geobacter, written by hand

units  per_hour
vars   acetate Fe3 Geobacter

rate   growth  = 0.0527 * acetate/(0.05 + acetate) * Fe3/(0.01 + Fe3)
rate   acetate = -1.0 * growth * Geobacter / 0.42
rate   Fe3     = -4.8 * growth * Geobacter / 0.42

# The box the law is valid over.  Enforced, not advisory: every evaluation is
# clamped to it and the clamping is counted, and the run reports at the end how
# often it happened.  A fitted rate law does not fail outside its range, it
# returns a confident number, so that count is the thing to read.
range  acetate 0.01  2.0
range  Fe3     0.005 1.0
'''


def case(number, name, xml, readme, geometry,
         kinetics=KIN.KINETICS_INERT, abiotic=KIN.ABIOTIC_INERT,
         model=None, extra_files=None):
    """Lay down one complete, self-contained example folder."""
    folder = os.path.join(ROOT, "%02d_%s" % (number, name))
    write(os.path.join(folder, "CompLaB.xml"), xml)
    write(os.path.join(folder, "README.md"), readme)
    write(os.path.join(folder, "defineKinetics.hh"), kinetics)
    write(os.path.join(folder, "defineAbioticKinetics.hh"), abiotic)
    os.makedirs(os.path.join(folder, "input"), exist_ok=True)
    shutil.copy(os.path.join(SHARED, geometry),
                os.path.join(folder, "input", "geometry.dat"))
    if model:
        write(os.path.join(folder, "input", "toy_model.xml"), TOY_MODEL)
    for fn, txt in (extra_files or {}).items():
        write(os.path.join(folder, fn), txt)
    print("  %02d_%s" % (number, name))
    return folder


def readme(number, title, what, expect, switches, notes=""):
    return """# Example %02d: %s

%s

## What to expect

%s

## The switches this case exercises

%s
%s
## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake %s.. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
""" % (number, title, what, expect, switches,
       ("\n## Notes\n\n%s\n" % notes) if notes else "",
       "-DENABLE_GLPK=ON " if number in (9, 12) else
       ("-DENABLE_COBRAPY=ON " if number == 10 else ""))


# ===========================================================================
def build():
    print("Building examples in %s" % ROOT)

    # -----------------------------------------------------------------------
    # 01  flow, and a conservative tracer riding on it
    # -----------------------------------------------------------------------
    xml = (K.header(1, "FLOW AND A CONSERVATIVE TRACER", [
        "The smallest thing that exercises both solvers: Navier Stokes for the",
        "flow, advection diffusion for a tracer that rides on it and reacts",
        "with nothing.",
        "",
        "biotic_mode is false, so the microbiology section is skipped entirely.",
        "Both kinetics hooks are off. If this case does not run, nothing else",
        "will, so it is the first thing runAllExamples.sh checks.",
        "",
        "Peclet 1 means advection and diffusion are comparable. The tracer",
        "should sweep from the left face to the right and level out.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>false</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=1, ade_max_iT=400) + """
    <chemistry>
        <number_of_substrates>1</number_of_substrates>
""" + K.substrate(0, "tracer", "0.", left=("Dirichlet", "1."),
                  right=("Dirichlet", "0."),
                  comment="held at 1 on the left, 0 on the right") + """
    </chemistry>

    <microbiology>
        <number_of_microbes>0</number_of_microbes>
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=100))

    case(1, "flow_only", xml, readme(
        1, "flow and a conservative tracer",
        "Navier Stokes plus advection diffusion, no chemistry and no biology.\n"
        "The tracer is held at 1 on the inlet face and 0 on the outlet, and\n"
        "reacts with nothing.",
        "- the flow solver converging, then the tracer sweeping left to right\n"
        "- a monotone profile: concentration falls with x, never rises\n"
        "- no `[NEG!]` warnings\n"
        "- every z plane identical, because z is periodic and the geometry is\n"
        "  extruded",
        "`biotic_mode` false, `<Peclet>` non zero, Dirichlet solute boundaries.",
        "This is the smoke test. If it fails, the build or the geometry is\n"
        "wrong and nothing further is worth debugging."),
        "channel.dat")

    # -----------------------------------------------------------------------
    # 02  pure diffusion, no flow
    # -----------------------------------------------------------------------
    xml = (K.header(2, "DIFFUSION ONLY, NO FLOW", [
        "The same tracer with the flow switched off, which is what",
        "Peclet 0 does. Nothing advects; the profile is pure diffusion.",
        "",
        "This one is worth having separately because it has an analytic answer.",
        "In a channel with one face held at 1 and the other at 0, and no",
        "reaction, the steady state is a straight line from 1 to 0 along x.",
        "Anything else means the diffusion coefficient, the boundary condition",
        "or the relaxation time is not doing what you think.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>false</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=0, ade_max_iT=2000) + """
    <chemistry>
        <number_of_substrates>1</number_of_substrates>
""" + K.substrate(0, "tracer", "0.", left=("Dirichlet", "1."),
                  right=("Dirichlet", "0.")) + """
    </chemistry>

    <microbiology>
        <number_of_microbes>0</number_of_microbes>
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=200))

    case(2, "diffusion_only", xml, readme(
        2, "diffusion only, no flow",
        "The example 01 tracer with `<Peclet>0</Peclet>`, so the flow solver\n"
        "never runs and transport is pure diffusion.",
        "- a **straight line** from 1 at the inlet to 0 at the outlet once it\n"
        "  converges. That is the analytic steady state for diffusion with no\n"
        "  reaction, and it is the check worth making\n"
        "- convergence reported before `ade_max_iT`\n"
        "- no z dependence",
        "`<Peclet>0</Peclet>`, which switches the flow solver off entirely.",
        "Curvature in the profile means a reaction is firing somewhere it\n"
        "should not, or a boundary is not being held. Both are worth chasing\n"
        "here, where there is nothing else going on to hide them."),
        "channel.dat")

    # -----------------------------------------------------------------------
    # 03  abiotic kinetics
    # -----------------------------------------------------------------------
    xml = (K.header(3, "ABIOTIC KINETICS:  A + B  ->  C", [
        "Chemistry with no biology. Two reactants enter from opposite faces,",
        "meet in the middle, and make a product.",
        "",
        "The rate law lives in defineAbioticKinetics.hh in this folder, which",
        "is COMPILED IN, not read at run time. That is why every example",
        "carries its own copy and why runAllExamples.sh rebuilds between cases.",
        "",
        "One A plus one B makes exactly one C, so the totals must satisfy",
        "dA + dC = 0 and dB + dC = 0 at every step. The runner checks it.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>false</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>true</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=0, ade_max_iT=1500) + """
    <chemistry>
        <number_of_substrates>3</number_of_substrates>
""" + K.substrate(0, "A", "0.", left=("Dirichlet", "1."), right=("Neumann", "0."),
                  comment="enters from the left")
       + K.substrate(1, "B", "0.", left=("Neumann", "0."), right=("Dirichlet", "1."),
                     comment="enters from the right")
       + K.substrate(2, "C", "0.", left=("Neumann", "0."), right=("Neumann", "0."),
                     comment="the product; no source, no sink at the walls") + """
    </chemistry>

    <microbiology>
        <number_of_microbes>0</number_of_microbes>
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=200))

    case(3, "abiotic_kinetics", xml, readme(
        3, "abiotic kinetics, A + B to C",
        "A second order irreversible reaction with no biology. A enters from\n"
        "the left, B from the right, and C forms where they overlap.",
        "- a peak in C near the middle of the domain, with A and B drawn down\n"
        "  around it\n"
        "- the reaction front stationary, because there is no flow to push it\n"
        "- **mass conserved**: total A + total C constant, total B + total C\n"
        "  constant",
        "`enable_abiotic_kinetics`, and a hand written rate law in\n"
        "`defineAbioticKinetics.hh`.",
        "Turn `<Peclet>` up to 1 and the front moves downstream. That single\n"
        "change is the cheapest way to see transport and reaction competing."),
        "channel.dat", abiotic=KIN.ABIOTIC_ABC)

    # -----------------------------------------------------------------------
    # 04  equilibrium speciation
    # -----------------------------------------------------------------------
    xml = (K.header(4, "AQUEOUS SPECIATION:  THE CARBONATE SYSTEM", [
        "The equilibrium solver, on the smallest tableau that is still real",
        "chemistry: bicarbonate, carbonate and carbonic acid over two",
        "components.",
        "",
        "    components        HCO3   H",
        "    HCO3   (component)   1    0     logK  0",
        "    H      (component)   0    1     logK  0",
        "    CO3                  1   -1     logK  -10.33",
        "    H2CO3                1    1     logK  6.35",
        "",
        "A component's own row is the identity row with logK 0. Species that",
        "do not take part get a row of zeros, which is how minerals and",
        "biomass are kept out of the tableau.",
        "",
        "KNOWN LIMITS, worth stating plainly: activities are ideal, there is",
        "no temperature dependence, no gas phase, no redox couple and no",
        "mineral saturation index. The logK values are conditional constants",
        "at your own ionic strength and 25 C.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>false</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=0, ade_max_iT=600) + """
    <chemistry>
        <number_of_substrates>4</number_of_substrates>
""" + K.substrate(0, "HCO3", "1e-3", left=("Dirichlet", "2e-3"), right=("Neumann", "0."),
                  comment="component: total carbonate")
       + K.substrate(1, "H", "1e-7", left=("Dirichlet", "1e-6"), right=("Neumann", "0."),
                     comment="component: protons, acid enters on the left")
       + K.substrate(2, "CO3", "0.", left=("Neumann", "0."), right=("Neumann", "0."))
       + K.substrate(3, "H2CO3", "0.", left=("Neumann", "0."), right=("Neumann", "0.")) + """
    </chemistry>

    <microbiology>
        <number_of_microbes>0</number_of_microbes>
    </microbiology>

    <equilibrium>
        <enabled>true</enabled>

        <components>HCO3 H</components>

        <stoichiometry>
            <species0>1  0</species0>
            <species1>0  1</species1>
            <species2>1 -1</species2>
            <species3>1  1</species3>
        </stoichiometry>

        <logK>
            <species0>0.</species0>
            <species1>0.</species1>
            <species2>-10.33</species2>
            <species3>6.35</species3>
        </logK>
    </equilibrium>
""" + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=100))

    case(4, "equilibrium", xml, readme(
        4, "aqueous speciation, the carbonate system",
        "Four species over two components, re-equilibrated in every pore voxel\n"
        "at every step. Acid enters from the left, so the speciation should\n"
        "shift along the domain.",
        "- the solver reporting the tableau it parsed: 2 components, 4 species\n"
        "- H2CO3 dominant where the acid is, CO3 dominant where it is not\n"
        "- **totals conserved**: total carbonate, summed over HCO3, CO3 and\n"
        "  H2CO3, changes only through transport, never through speciation\n"
        "- convergence of the Anderson accelerated iteration reported each step",
        "`<equilibrium><enabled>true`, plus the `<components>`,\n"
        "`<stoichiometry>` and `<logK>` tableau.",
        "This is the most expensive part of the code by a wide margin, which\n"
        "is why every other example leaves it off. Even on this 3456 voxel\n"
        "domain you will see it in the step time."),
        "channel.dat")

    # -----------------------------------------------------------------------
    # 05, 06, 07  the three biomass solvers, same chemistry
    # -----------------------------------------------------------------------
    for n, solver, tag, blurb, extra_note in [
        (5, "CA", "cellular_automaton",
         "the cellular automaton: biomass that exceeds the maximum density\n"
         "     spills into neighbouring voxels, which is how a biofilm spreads",
         "`CA_method` picks how the excess is shared. `fraction` splits it in\n"
         "proportion to the space available in each neighbour."),
        (6, "FD", "finite_difference",
         "finite difference: biomass diffuses with its own coefficient,\n"
         "     solved on the same grid",
         "FD is the only solver that requires `<biomass_diffusion_coefficients>`,\n"
         "and it is rejected outright for planktonic microbes: a microbe with\n"
         "no `<microbeN>` entry under `material_numbers` cannot use it."),
        (7, "LBM", "lattice_boltzmann",
         "lattice Boltzmann: biomass gets its own advection diffusion\n"
         "     lattice, the same machinery the solutes use",
         "LBM costs a full extra lattice per microbe. It is the right choice\n"
         "when biomass genuinely moves with the flow, and overkill when it\n"
         "does not."),
    ]:
        visc = "0" if solver == "CA" else "0"
        dbio = "3e-13" if solver in ("FD", "LBM") else None
        xml = (K.header(n, "BIOTIC KINETICS, %s BIOMASS SOLVER" % solver, [
            "One organism growing on one donor by Monod kinetics, with the",
            "biomass transported by %s." % solver,
            "",
            "Examples 05, 06 and 07 are THE SAME CASE with only",
            "<solver_type> changed. Run all three and compare: the biomass",
            "totals should be close, the spatial pattern different. That is",
            "the honest way to see what the solver choice actually costs you.",
            "",
            "The rate law is in defineKinetics.hh in this folder:",
            "",
            "    growth  = mu_max * S/(Ks + S) * B",
            "    donor   consumed at growth / Y",
            "    product made at growth * fP",
            "    decay   kd * B, first order",
            "",
            "Note that <half_saturation_constants> in the XML is NOT what",
            "drives this. A hand written kinetics file takes its parameters",
            "from its own table; the XML tag feeds the FBA and surrogate",
            "paths. The tag is present here because the parser expects it.",
        ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>true</biotic_mode>
        <enable_kinetics>true</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=0, ade_max_iT=1000,
               materials="\n                <microbe0>3</microbe0>") + """
    <chemistry>
        <number_of_substrates>2</number_of_substrates>
""" + K.substrate(0, "donor", "0.1", left=("Dirichlet", "0.5"), right=("Neumann", "0."),
                  comment="fed from the left so the colony does not starve")
           + K.substrate(1, "product", "0.", left=("Neumann", "0."), right=("Neumann", "0.")) + """
    </chemistry>

    <microbiology>
        <number_of_microbes>1</number_of_microbes>
        <maximum_biomass_density>50</maximum_biomass_density>
        <thrd_biofilm_fraction>0.01</thrd_biofilm_fraction>
        <CA_method>fraction</CA_method>
""" + K.microbe(0, "Bug", solver, "kinetics", "1.0", decay="0.0",
                visc=visc, Dbio=dbio, Kc="0.05 0", Vmax="0 0",
                comment="seeded on material number 3") + """
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=200))

        case(n, "biotic_%s" % tag, xml, readme(
            n, "biotic kinetics, %s biomass solver" % solver,
            "One organism, Monod growth on one donor, biomass transported by %s.\n"
            "Examples 05, 06 and 07 differ only in `<solver_type>`." % solver,
            "- biomass rising from the seeded patch, then levelling off as\n"
            "  growth and decay balance\n"
            "- the donor drawn down around the colony and replenished from the\n"
            "  left boundary\n"
            "- a product plume spreading from the colony\n"
            "- no `[NEG!]` warnings: the donor must never go negative",
            "`<solver_type>%s</solver_type>`, %s" % (solver, blurb),
            extra_note),
            "one_microbe.dat", kinetics=KIN.KINETICS_MONOD)

    # -----------------------------------------------------------------------
    # 08  two microbes, two different solvers, competing
    # -----------------------------------------------------------------------
    xml = (K.header(8, "TWO MICROBES COMPETING, TWO DIFFERENT SOLVERS", [
        "Two organisms on opposite walls, both eating the same donor, each",
        "using a different biomass solver.",
        "",
        "This is the case that proves <solver_type> and <reaction_type> are",
        "PER MICROBE, not global. Bug1 uses the cellular automaton, Bug2 uses",
        "lattice Boltzmann, and they run in the same simulation.",
        "",
        "They differ in kinetics too, in defineKinetics.hh: Bug1 has the",
        "higher maximum growth rate, Bug2 the lower half saturation constant.",
        "So Bug1 wins where the donor is plentiful and Bug2 wins where it is",
        "scarce, which is the classic r versus K trade off. Bug1 sits nearer",
        "the inlet, so it should end up ahead.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>true</biotic_mode>
        <enable_kinetics>true</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=1, ade_max_iT=1000,
               materials="\n                <microbe0>3</microbe0>"
                         "\n                <microbe1>4</microbe1>") + """
    <chemistry>
        <number_of_substrates>2</number_of_substrates>
""" + K.substrate(0, "donor", "0.1", left=("Dirichlet", "0.5"), right=("Neumann", "0."))
       + K.substrate(1, "product", "0.", left=("Neumann", "0."), right=("Neumann", "0.")) + """
    </chemistry>

    <microbiology>
        <number_of_microbes>2</number_of_microbes>
        <maximum_biomass_density>50</maximum_biomass_density>
        <thrd_biofilm_fraction>0.01</thrd_biofilm_fraction>
        <CA_method>fraction</CA_method>
""" + K.microbe(0, "Bug1_fast", "CA", "kinetics", "1.0", visc="0",
                Kc="0.05 0", Vmax="0 0",
                comment="high mu_max, high Ks: wins where food is plentiful")
       + K.microbe(1, "Bug2_efficient", "LBM", "kinetics", "1.0", visc="0",
                   Dbio="3e-13", Kc="0.005 0", Vmax="0 0",
                   comment="low mu_max, low Ks: wins where food is scarce") + """
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=200))

    case(8, "two_microbes", xml, readme(
        8, "two microbes competing, two different solvers",
        "Two organisms on opposite walls sharing one donor. Bug1 uses the\n"
        "cellular automaton and grows fast; Bug2 uses lattice Boltzmann and is\n"
        "more efficient at low substrate.",
        "- both populations growing, then their **ratio going flat**\n"
        "- Bug1 ahead, because it sits nearer the inlet where the donor is\n"
        "  plentiful and its higher mu_max pays off\n"
        "- a donor shadow between the two colonies",
        "Per microbe `<solver_type>` and `<reaction_type>`, two `<microbeN>`\n"
        "material numbers, `bfilm_count` of 2.",
        "Swap the two `<solver_type>` lines and rerun. The ratio should barely\n"
        "move: the solver decides how biomass spreads, not how much of it\n"
        "there is. If it moves a lot, that is worth understanding before you\n"
        "trust either solver on a real problem."),
        "two_microbes.dat", kinetics=KIN.KINETICS_MONOD)

    # -----------------------------------------------------------------------
    # 09 / 10  flux balance analysis, GLPK and COBRApy
    # -----------------------------------------------------------------------
    for n, tag, solver_name, rxn, switch in [
        (9, "fba_glpk", "GLPK", "glpk", "enable_fba_glpk"),
        (10, "fba_cobrapy", "COBRApy", "cobrapy", "enable_fba_cobrapy"),
    ]:
        others = {"enable_fba_glpk": "false", "enable_fba_cobrapy": "false"}
        others[switch] = "true"
        xml = (K.header(n, "FLUX BALANCE ANALYSIS THROUGH %s" % solver_name, [
            "One organism whose growth comes from a linear program solved in",
            "every voxel, at every step, instead of from a hand written rate",
            "law.",
            "",
            "THE MODEL IS SMALL ENOUGH TO SOLVE BY HAND. input/toy_model.xml",
            "has three metabolites and four reactions:",
            "",
            "    0 EX_S   S_c ->             negative flux imports the donor",
            "    1 EX_O   O_c ->             negative flux imports the acceptor",
            "    2 EX_P   P_c ->             positive flux exports the product",
            "    3 BIO    S_c + 0.5 O_c -> P_c        the objective",
            "",
            "Steady state forces v0 = -v3 and v1 = -0.5 v3, so",
            "",
            "    growth = min( Vs, 2 * Vo )",
            "",
            "where each V is the Michaelis Menten bound CompLaB builds from",
            "the local concentration:  V = Vmax * C / (Kc + C).",
            "",
            "With Vmax 10 and 4, Kc 0.05 and 0.05, and both concentrations",
            "well above their half saturation constants, Vs approaches 10 and",
            "2 * Vo approaches 8, so the acceptor limits and growth settles",
            "near 8 mmol/gDW/h. THAT NUMBER COMES FROM THE STOICHIOMETRY, NOT",
            "FROM THIS CODE, which is what makes it a test.",
            "",
            "Examples 09 and 10 are the same case on the two different",
            "solvers. They receive an identical linear program, so they must",
            "return the same answer. If they disagree, one of the two",
            "couplings is wrong.",
        ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>true</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>%s</enable_fba_glpk>
        <enable_fba_cobrapy>%s</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <fba_concentration_basis>free</fba_concentration_basis>
        <glpk_lpsolver>1</glpk_lpsolver>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" % (others["enable_fba_glpk"], others["enable_fba_cobrapy"])
            + K.domain("geometry.dat", peclet=0, ade_max_iT=500,
                       materials="\n                <microbe0>3</microbe0>") + """
    <chemistry>
        <number_of_substrates>3</number_of_substrates>
        <fix_concentration>no no no</fix_concentration>
        <fix_lower_bounds>no no no</fix_lower_bounds>
""" + K.substrate(0, "S", "1.0", left=("Dirichlet", "1.0"), right=("Neumann", "0."),
                  comment="donor, held up from the left")
            + K.substrate(1, "O", "1.0", left=("Dirichlet", "1.0"), right=("Neumann", "0."),
                          comment="acceptor, held up from the left")
            + K.substrate(2, "P", "0.", left=("Neumann", "0."), right=("Neumann", "0."),
                          comment="product, excreted") + """
    </chemistry>

    <microbiology>
        <number_of_microbes>1</number_of_microbes>
        <maximum_biomass_density>50</maximum_biomass_density>
        <thrd_biofilm_fraction>0.01</thrd_biofilm_fraction>
        <CA_method>fraction</CA_method>
""" + K.microbe(0, "Toybug", "FD", rxn, "1.0", decay="0.01",
                visc="0", Dbio="3e-13", Kc="0.05 0.05 0", extra="""            <model_filename>toy_model</model_filename>

            <!-- which reaction exchanges each substrate, in substrate order.
                 Use -1 for a substrate this organism does not exchange. -->
            <exchange_reaction_indices>0 1 2</exchange_reaction_indices>

            <!-- Michaelis Menten Vmax, mmol/gDW/h, one per substrate.
                 Vs approaches 10 and 2 * Vo approaches 8, so the acceptor
                 limits and growth settles near 8. -->
            <fba_maximum_uptake_flux>10. 4. 0.</fba_maximum_uptake_flux>

            <!-- hard floor and ceiling on the exchange fluxes.  Exchange
                 reactions are written "metabolite ->", so eating is negative. -->
            <substrate_lower_bounds>-10. -4. 0.</substrate_lower_bounds>
            <substrate_upper_bounds>100. 100. 100.</substrate_upper_bounds>

            <objective_direction>maximize</objective_direction>

            <!-- grams dry weight per mole of cells: converts CompLB3D's molar
                 biomass into the gDW a metabolic model expects. -->
            <biomass_molar_mass>24.6</biomass_molar_mass>
""") + """
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=100))

        case(n, tag, xml, readme(
            n, "flux balance analysis through %s" % solver_name,
            "One organism on a four reaction toy model, small enough that the\n"
            "answer is a formula:\n\n"
            "```\n"
            "growth = min( Vs, 2 * Vo ),   V = Vmax * C / (Kc + C)\n"
            "```\n\n"
            "With Vmax 10 and 4 and both concentrations well above their half\n"
            "saturation constants, the acceptor limits and growth should settle\n"
            "near **8 mmol/gDW/h**.",
            "- the model loading as 3 metabolites, 4 reactions\n"
            "- the metabolic summary showing `%s` with 1 microbe\n"
            "- the growth rate approaching 8, from the stoichiometry rather\n"
            "  than from this code\n"
            "- product P accumulating at the same rate biomass is made, since\n"
            "  the biomass reaction makes exactly one P per unit growth\n"
            "- no `[NEG!]` warnings" % solver_name,
            "`<%s>true`, `<reaction_type>%s</reaction_type>`, and the\n"
            "`<model_filename>` / `<exchange_reaction_indices>` /\n"
            "`<fba_maximum_uptake_flux>` group." % (switch, rxn),
            "Examples 09 and 10 must agree. That comparison is the strongest\n"
            "check available on the FBA coupling, because the two solvers share\n"
            "nothing but the linear program.\n\n"
            + ("COBRApy needs `src/complab3d_cobrapy.py` present at `<src_path>`\n"
               "at run time; it is imported by name, not compiled in. Expect one\n"
               "to two orders of magnitude slower than GLPK."
               if n == 10 else
               "Drop `<fba_maximum_uptake_flux>` to `4. 10. 0.` and the donor\n"
               "becomes limiting instead, so growth should settle near 4. That\n"
               "is a second independent check for one line of editing.")),
            "one_microbe.dat", model=True)

    # -----------------------------------------------------------------------
    # 11  surrogate network
    # -----------------------------------------------------------------------
    xml = (K.header(11, "SURROGATE MODEL: A TRAINED NETWORK INSTEAD OF FBA", [
        "The same idea as examples 09 and 10, with the linear program replaced",
        "by a small neural network fitted to it offline. A few hundred",
        "multiplications instead of a solve, so 100 to 1000 times faster.",
        "",
        "This case uses the network SHIPPED in surrogateModel.hh: Geobacter",
        "metallireducens on acetate and Fe(III), 2 inputs, four hidden layers",
        "of ten tansig neurons, one linear output.",
        "",
        "THE VMAX VALUES BELOW ARE NOT ARBITRARY. A network is only valid",
        "inside the range it was trained on, and outside it a network does not",
        "fail, it lies. The shipped one was trained on",
        "",
        "    acetate   0.00089 to 9.998   mmol/gDW/h",
        "    Fe(III)   0.000035 to 0.4997 mmol/gDW/h",
        "",
        "so Vmax is set to 8 and 0.4, safely inside. Those bounds were",
        "recovered by inverting the mapminmax scaling; run",
        "surrogate_training/inspectSurrogate.py on any network to get them.",
        "",
        "Raise <fba_maximum_uptake_flux> above 0.5 on the second substrate and",
        "the simulation leaves the training box. It will still produce",
        "numbers. They will not mean anything.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>true</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>true</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=0, ade_max_iT=500,
               materials="\n                <microbe0>3</microbe0>") + """
    <chemistry>
        <number_of_substrates>2</number_of_substrates>
""" + K.substrate(0, "acetate", "1.0", left=("Dirichlet", "1.0"), right=("Neumann", "0."))
       + K.substrate(1, "Fe3", "0.5", left=("Dirichlet", "0.5"), right=("Neumann", "0.")) + """
    </chemistry>

    <microbiology>
        <number_of_microbes>1</number_of_microbes>
        <maximum_biomass_density>50</maximum_biomass_density>
        <thrd_biofilm_fraction>0.01</thrd_biofilm_fraction>
        <CA_method>fraction</CA_method>
""" + K.microbe(0, "Geobacter", "FD", "surrogate", "1.0", decay="0.01",
                visc="0", Dbio="3e-13", Kc="0.05 0.01", extra="""
            <!-- Vmax, mmol/gDW/h.  Both are INSIDE the shipped network's
                 training range; see the banner above before changing them. -->
            <fba_maximum_uptake_flux>8. 0.4</fba_maximum_uptake_flux>

            <biomass_molar_mass>24.6</biomass_molar_mass>
""") + """
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=100))

    case(11, "surrogate", xml, readme(
        11, "surrogate model, a trained network instead of FBA",
        "The shipped Geobacter network standing in for a flux balance solve.\n"
        "Two inputs, acetate and Fe(III) uptake; one output, growth rate.",
        "- growth rates in the low 1e-2 1/h range: the shipped network's whole\n"
        "  output range is 0 to 0.0527 1/h\n"
        "- a step time far below examples 09 and 10 on the same domain, which\n"
        "  is the entire reason to use a surrogate\n"
        "- both substrates drawn down together",
        "`<enable_surrogate>true`, `<reaction_type>surrogate</reaction_type>`,\n"
        "and `surrogateModel.hh` compiled in.",
        "**The Vmax values matter here in a way they do not elsewhere.** The\n"
        "shipped network has never seen a Fe(III) uptake above 0.5 mmol/gDW/h.\n"
        "`<fba_maximum_uptake_flux>` is set to 0.4 to stay inside that. Raise\n"
        "it and the run leaves the training box, where a network extrapolates\n"
        "badly and returns confident nonsense.\n\n"
        "To fit your own network, see `surrogate_training/` in the parent\n"
        "folder. `inspectSurrogate.py` prints any network's valid range."),
        "one_microbe.dat")

    # -----------------------------------------------------------------------
    # 12  combined reaction type
    # -----------------------------------------------------------------------
    xml = (K.header(12, "COMBINED REACTION TYPES: FBA PLUS HAND WRITTEN KINETICS", [
        "One organism running BOTH a metabolic solver and a hand written rate",
        "law, with the two rates added. That is what the combined reaction",
        "types are for:",
        "",
        "    glpk_and_kinetics        surrogate_and_kinetics",
        "    cobrapy_and_kinetics",
        "",
        "Flux balance analysis supplies the growth; defineKinetics.hh supplies",
        "what the metabolic model does not describe. Here that is MAINTENANCE:",
        "a constant per biomass drain on the donor that yields no growth. Real",
        "cells pay it and a plain FBA objective does not include it, so the",
        "model over predicts growth at low substrate without it.",
        "",
        "THE SIGN DISCIPLINE MATTERS. defineKinetics.hh in this folder writes",
        "subsR but deliberately leaves bioR at zero, because growth belongs to",
        "the FBA path and writing it in both places would double count. That",
        "is the classic mistake with combined types.",
        "",
        "Note also the two separate Vmax tags. <maximum_uptake_flux> is read",
        "by the kinetics path in mol per mol cells per second;",
        "<fba_maximum_uptake_flux> is read by the metabolic path in",
        "mmol/gDW/h. One number cannot satisfy both meanings, so a combined",
        "microbe needs both tags.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>true</biotic_mode>
        <enable_kinetics>true</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>true</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <fba_concentration_basis>free</fba_concentration_basis>
        <glpk_lpsolver>1</glpk_lpsolver>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=0, ade_max_iT=500,
               materials="\n                <microbe0>3</microbe0>") + """
    <chemistry>
        <number_of_substrates>3</number_of_substrates>
        <fix_concentration>no no no</fix_concentration>
        <fix_lower_bounds>no no no</fix_lower_bounds>
""" + K.substrate(0, "S", "1.0", left=("Dirichlet", "1.0"), right=("Neumann", "0."))
       + K.substrate(1, "O", "1.0", left=("Dirichlet", "1.0"), right=("Neumann", "0."))
       + K.substrate(2, "P", "0.", left=("Neumann", "0."), right=("Neumann", "0.")) + """
    </chemistry>

    <microbiology>
        <number_of_microbes>1</number_of_microbes>
        <maximum_biomass_density>50</maximum_biomass_density>
        <thrd_biofilm_fraction>0.01</thrd_biofilm_fraction>
        <CA_method>fraction</CA_method>
""" + K.microbe(0, "Toybug", "FD", "glpk_and_kinetics", "1.0", decay="0.01",
                visc="0", Dbio="3e-13", Kc="0.05 0.05 0",
                Vmax="0. 0. 0.", extra="""            <model_filename>toy_model</model_filename>
            <exchange_reaction_indices>0 1 2</exchange_reaction_indices>

            <!-- the metabolic path's own Vmax, mmol/gDW/h.  Separate from
                 <maximum_uptake_flux> above, which the kinetics path reads in
                 different units. -->
            <fba_maximum_uptake_flux>10. 4. 0.</fba_maximum_uptake_flux>

            <substrate_lower_bounds>-10. -4. 0.</substrate_lower_bounds>
            <substrate_upper_bounds>100. 100. 100.</substrate_upper_bounds>
            <objective_direction>maximize</objective_direction>
            <biomass_molar_mass>24.6</biomass_molar_mass>
""") + """
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + K.io(vtk=100))

    case(12, "mixed_reaction_types", xml, readme(
        12, "combined reaction types, FBA plus hand written kinetics",
        "The example 09 organism with `<reaction_type>glpk_and_kinetics`, so\n"
        "the linear program and `defineKinetics.hh` both run and their rates\n"
        "add. The kinetics file contributes maintenance: a donor drain that\n"
        "yields no growth.",
        "- the same growth rate as example 09, near 8 mmol/gDW/h, because\n"
        "  maintenance does not add biomass\n"
        "- the donor drawn down **faster** than in example 09, by exactly the\n"
        "  maintenance term\n"
        "- the metabolic summary showing 1 GLPK microbe, and the kinetics\n"
        "  path also active",
        "`<reaction_type>glpk_and_kinetics</reaction_type>`, both\n"
        "`<enable_kinetics>` and `<enable_fba_glpk>` true, and both\n"
        "`<maximum_uptake_flux>` and `<fba_maximum_uptake_flux>` present.",
        "Diff this case's donor profile against example 09's. The difference\n"
        "IS the maintenance term, and it is the cheapest way to confirm the\n"
        "two paths are genuinely adding rather than one silently winning."),
        "one_microbe.dat", kinetics=KIN.KINETICS_MAINTENANCE, model=True)

    # -----------------------------------------------------------------------
    # 13  precipitation and pore clogging
    # -----------------------------------------------------------------------
    xml = (K.header(13, "PRECIPITATION AND PORE CLOGGING", [
        "A mineral forms where two solutes meet, fills its voxels, and blocks",
        "the flow.",
        "",
        "    Fe2+ + HS-  ->  FeS(s)",
        "",
        "FeS is declared <immobile>true</immobile>, so it never advects or",
        "diffuses: it accumulates exactly where it forms. Once a voxel reaches",
        "max_precipRho it is converted to solid and the flow is re-solved",
        "around it.",
        "",
        "max_precipRho comes from the mineral, not from the model: it is",
        "1000 divided by the molar volume in cm3/mol. Mackinawite is",
        "87.91 g/mol over 4.30 g/cm3, so 20.44 cm3/mol, so 48.9 mol/L. That",
        "number is used below.",
        "",
        "surface_only 1 means growth only on wall adjacent voxels, which is",
        "heterogeneous nucleation and is the physically right choice. Set it",
        "to 0 and the mineral grows anywhere, which seals the pore throat far",
        "too readily.",
        "",
        "KNOWN LIMIT: the conversion is one way. Nothing here reopens a sealed",
        "voxel. Example 14 does that, and example 15 does both.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>false</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>true</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=1, ade_max_iT=2000) + """
    <chemistry>
        <number_of_substrates>3</number_of_substrates>
""" + K.substrate(0, "Fe2", "0.", left=("Dirichlet", "5."), right=("Neumann", "0."),
                  comment="enters from the left")
       + K.substrate(1, "HS", "0.", left=("Neumann", "0."), right=("Dirichlet", "5."),
                     comment="enters from the right")
       + K.substrate(2, "FeS", "0.", D="0.", left=("Neumann", "0."), right=("Neumann", "0."),
                     immobile=True, comment="the mineral: immobile, never transported") + """
    </chemistry>

    <microbiology>
        <number_of_microbes>0</number_of_microbes>
    </microbiology>
""" + K.OFF_EQ + """
    <precipitation>
        <enabled>true</enabled>
        <solid_substrate>2</solid_substrate>     <!-- index of the mineral substrate -->
        <max_precipRho>48.9</max_precipRho>      <!-- mol/L in a completely full voxel -->
        <surface_only>1</surface_only>           <!-- 1 = heterogeneous nucleation -->
        <perm_ratio>0</perm_ratio>               <!-- 0 = a filled voxel is an impermeable wall -->
        <update_interval>100</update_interval>   <!-- steps between geometry re-checks -->
    </precipitation>
""" + K.OFF_DISSOL + K.io(vtk=200))

    case(13, "precipitation", xml, readme(
        13, "precipitation and pore clogging",
        "FeS forms where an iron front and a sulfide front meet, accumulates\n"
        "in place because it is declared immobile, and seals its voxels once\n"
        "it reaches the full mineral density.",
        "- FeS building up in a band across the middle\n"
        "- the porosity reported in the log **falling in steps**, one per\n"
        "  update interval\n"
        "- the flow re-solved after each conversion, and the velocity rising\n"
        "  in the remaining open path\n"
        "- **mass conserved**: every mole of FeS formed removes one mole of\n"
        "  Fe2+ and one of HS-",
        "`<precipitation><enabled>true`, `<immobile>true</immobile>` on the\n"
        "mineral substrate, and the rate law in `defineAbioticKinetics.hh`.",
        "`max_precipRho` is 48.9 mol/L, which is 1000 divided by mackinawite's\n"
        "molar volume of 20.44 cm3/mol. Derive it from your own mineral rather\n"
        "than copying this number.\n\n"
        "If the domain seals completely the run stops and reports a percolation\n"
        "limit. That is the code working, not failing."),
        "grains.dat", abiotic=KIN.ABIOTIC_PRECIP)

    # -----------------------------------------------------------------------
    # 14  dissolution and pore reopening
    # -----------------------------------------------------------------------
    xml = (K.header(14, "DISSOLUTION AND PORE REOPENING", [
        "The mirror image of example 13: acid eats the grains away and their",
        "voxels open back up to flow.",
        "",
        "    CaCO3(s) + H+  ->  Ca2+ + HCO3-",
        "",
        "WHAT CAN DISSOLVE is only what is declared in a <phaseN> block, and",
        "that is deliberate. If the rule were simply that a solid voxel with",
        "no mineral left becomes pore, then an inert quartz grain, which never",
        "held any mineral, is already below any threshold and the whole grain",
        "pack would dissolve on the first step. Declaring phases makes quartz",
        "safe by construction.",
        "",
        "initial_fill is the whole point of this case: the four grains START",
        "solid, at full density, and are eaten away. Dissolution is not",
        "limited to material that precipitated first.",
        "",
        "The rate law is defineDissolutionRate() in defineAbioticKinetics.hh,",
        "a separate hook from the precipitation one. It is handed the water",
        "TOUCHING the mineral, averaged over the open neighbours, because a",
        "solid voxel has no water of its own that moves. A voxel with no open",
        "neighbour is unreachable and does not dissolve at all.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>false</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>true</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=1, ade_max_iT=2000) + """
    <chemistry>
        <number_of_substrates>3</number_of_substrates>
""" + K.substrate(0, "H", "1e-4", left=("Dirichlet", "1e-2"), right=("Neumann", "0."),
                  comment="acid, entering from the left")
       + K.substrate(1, "Ca", "0.", left=("Neumann", "0."), right=("Neumann", "0."),
                     comment="released by dissolution")
       + K.substrate(2, "calcite", "0.", D="0.", left=("Neumann", "0."), right=("Neumann", "0."),
                     immobile=True, comment="the mineral inventory of each grain") + """
    </chemistry>

    <microbiology>
        <number_of_microbes>0</number_of_microbes>
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + """
    <dissolution>
        <enabled>true</enabled>
        <reopen_fraction>0.9</reopen_fraction>   <!-- reopen below 0.9 x full -->
        <surface_only>1</surface_only>           <!-- only wetted voxels dissolve -->
        <update_interval>100</update_interval>

        <phase0>
            <name>calcite_grain</name>
            <material_number>0</material_number>  <!-- the solid grains in geometry.dat -->
            <substrate>2</substrate>              <!-- which substrate holds its inventory -->
            <full_density>27.1</full_density>     <!-- mol/L: 1000 / 36.9 cm3 per mol -->
            <initial_fill>27.1</initial_fill>     <!-- starts completely solid -->
        </phase0>
    </dissolution>
""" + K.io(vtk=200))

    case(14, "dissolution", xml, readme(
        14, "dissolution and pore reopening",
        "Acid enters from the left and eats into four calcite grains. As each\n"
        "grain thins, its outer voxels drop below `reopen_fraction` of full\n"
        "density and reopen to flow.",
        "- the acid front eating into the **upstream** faces of the grains\n"
        "  first, because that is where the acid arrives\n"
        "- Ca2+ leaving on the right\n"
        "- the porosity in the log **rising in steps**, the mirror of\n"
        "  example 13\n"
        "- **mass conserved**: every mole of calcite removed appears as a mole\n"
        "  of Ca2+ in the water",
        "`<dissolution><enabled>true`, one `<phase0>` block with a non zero\n"
        "`<initial_fill>`, and `defineDissolutionRate()` in\n"
        "`defineAbioticKinetics.hh`.",
        "`<initial_fill>27.1</initial_fill>` is what makes the ORIGINAL grains\n"
        "dissolvable rather than only material that precipitated first. A\n"
        "precipitate phase would use `<initial_fill>0</initial_fill>` and\n"
        "`<is_precipitate>true</is_precipitate>` instead.\n\n"
        "Nothing precipitates in this case. `defineAbioticRxnKinetics` is\n"
        "empty on purpose, so the only thing moving mineral is the dissolution\n"
        "hook and the mass balance is unambiguous."),
        "grains.dat", abiotic=KIN.ABIOTIC_DISSOL)

    # -----------------------------------------------------------------------
    # 15  both directions at once
    # -----------------------------------------------------------------------
    xml = (K.header(15, "PRECIPITATION AND DISSOLUTION TOGETHER", [
        "Both directions active on the same mineral, so a voxel can seal and",
        "later reopen.",
        "",
        "    precipitation   Fe2+ + HS-  ->  FeS(s)",
        "    dissolution     FeS(s) + H+ ->  Fe2+ + HS-",
        "",
        "The two blocks are independent, which is why they can both be on.",
        "The <phase0> block here has is_precipitate true and initial_fill 0:",
        "there is no mineral anywhere at the start, it grows by reaction, and",
        "the voxels it seals are stamped so they can later reopen.",
        "",
        "WHY reopen_fraction EARNS ITS KEEP HERE AND NOWHERE ELSE. With both",
        "directions active a voxel can sit near the sealing threshold. A",
        "single threshold would make it flip open and shut every update",
        "interval, and every flip re-solves the flow. 0.9 puts a gap between",
        "sealing at full density and reopening below 0.9 times full, so it",
        "settles instead of chattering. Watch the porosity: it should move in",
        "steps and then hold.",
        "",
        "THE TWO RATE LAWS ARE NOT INVERSES. Precipitation is second order in",
        "the dissolved ions; dissolution is first order in acid and in",
        "remaining mineral. They are separate laws that happen to move the",
        "same mineral, which is how the code is meant to be used: the",
        "equilibrium solver does not compute solubility products, so there is",
        "no saturation state to drive a single reversible rate.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>false</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>true</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=1, ade_max_iT=3000) + """
    <chemistry>
        <number_of_substrates>4</number_of_substrates>
""" + K.substrate(0, "Fe2", "0.", left=("Dirichlet", "5."), right=("Neumann", "0."))
       + K.substrate(1, "HS", "0.", left=("Neumann", "0."), right=("Dirichlet", "5."))
       + K.substrate(2, "H", "0.", left=("Dirichlet", "1e-2"), right=("Neumann", "0."),
                     comment="acid, which drives the reverse reaction")
       + K.substrate(3, "FeS", "0.", D="0.", left=("Neumann", "0."), right=("Neumann", "0."),
                     immobile=True, comment="the mineral") + """
    </chemistry>

    <microbiology>
        <number_of_microbes>0</number_of_microbes>
    </microbiology>
""" + K.OFF_EQ + """
    <precipitation>
        <enabled>true</enabled>
        <solid_substrate>3</solid_substrate>
        <max_precipRho>48.9</max_precipRho>
        <surface_only>1</surface_only>
        <perm_ratio>0</perm_ratio>
        <update_interval>100</update_interval>
    </precipitation>

    <dissolution>
        <enabled>true</enabled>
        <reopen_fraction>0.9</reopen_fraction>
        <surface_only>1</surface_only>
        <update_interval>100</update_interval>

        <phase0>
            <name>FeS_precipitate</name>
            <material_number>0</material_number>
            <substrate>3</substrate>
            <full_density>48.9</full_density>
            <initial_fill>0</initial_fill>          <!-- nothing there at the start -->
            <is_precipitate>true</is_precipitate>   <!-- voxels sealed by precipitation -->
        </phase0>
    </dissolution>
""" + K.io(vtk=200))

    case(15, "precip_and_dissolution", xml, readme(
        15, "precipitation and dissolution together",
        "Both directions on the same mineral. FeS forms where the iron and\n"
        "sulfide fronts meet and seals its voxels; acid arriving from the left\n"
        "eats it back out again and those voxels reopen.",
        "- porosity **falling, then rising** as the acid front catches up with\n"
        "  the seal\n"
        "- the porosity settling rather than oscillating. Chattering means\n"
        "  `reopen_fraction` is too close to 1\n"
        "- **mass conserved across both directions**: iron and sulfur are only\n"
        "  moved between the water and the mineral, never created",
        "`<precipitation>` and `<dissolution>` both enabled, on one mineral,\n"
        "with `<is_precipitate>true</is_precipitate>` linking them.",
        "This is the case that justifies `reopen_fraction`. Set it to 0.99 and\n"
        "rerun: you should see the porosity flip back and forth every update\n"
        "interval, with a flow solve wasted on each flip. That is the failure\n"
        "the hysteresis exists to prevent."),
        "grains.dat", abiotic=KIN.ABIOTIC_BOTH)



    # -----------------------------------------------------------------------
    # 16. a rate law read from a file, instead of compiled in
    # -----------------------------------------------------------------------
    xml = (K.header(16, "symbolic rate law from a .sym file", [
        "The same organism and the same two substrates as example 11, but the",
        "rate law is read from input/rates.sym at start-up instead of being",
        "compiled into defineKinetics.hh.",
        "",
        "Edit rates.sym and rerun. No rebuild. That is the whole point: a",
        "fitting run produces a formula, the formula goes in a text file, and",
        "the next simulation uses it, so comparing four candidate mechanisms is",
        "four files and four runs rather than four edits and four compiles.",
        "",
        "The file is echoed into the log before the run starts, so a result",
        "carries the rate law that produced it.",
    ]) + K.path() + """
    <simulation_mode>
        <biotic_mode>true</biotic_mode>
        <enable_kinetics>false</enable_kinetics>
        <enable_abiotic_kinetics>false</enable_abiotic_kinetics>
        <enable_fba_glpk>false</enable_fba_glpk>
        <enable_fba_cobrapy>false</enable_fba_cobrapy>
        <enable_surrogate>false</enable_surrogate>
        <enable_validation_diagnostics>true</enable_validation_diagnostics>
    </simulation_mode>
""" + K.domain("geometry.dat", peclet=0, ade_max_iT=500,
               materials="\n                <microbe0>3</microbe0>") + """
    <chemistry>
        <number_of_substrates>2</number_of_substrates>
""" + K.substrate(0, "acetate", "1.0", left=("Dirichlet", "1.0"), right=("Neumann", "0."))
         + K.substrate(1, "Fe3", "0.5", left=("Dirichlet", "0.5"), right=("Neumann", "0.")) + """
    </chemistry>

    <microbiology>
        <number_of_microbes>1</number_of_microbes>
        <maximum_biomass_density>50</maximum_biomass_density>
        <thrd_biofilm_fraction>0.01</thrd_biofilm_fraction>
        <CA_method>fraction</CA_method>
""" + K.microbe(0, "Geobacter", "FD", "symbolic", "1.0", decay="0.01",
                visc="0", Dbio="3e-13", Kc="0.05 0.01") + """
    </microbiology>
""" + K.OFF_EQ + K.OFF_PRECIP + K.OFF_DISSOL + """
    <!-- The names in rates.sym must match the substrate and microbe names above
         exactly. A name that does not resolve stops the run and says which one,
         rather than binding silently to the wrong lattice. -->
    <symbolic>
        <enabled>true</enabled>
        <expressions_file>input/rates.sym</expressions_file>
    </symbolic>
""" + K.io(vtk=100))

    case(16, "learned_rate_laws", xml, readme(
        16, "a rate law read from a file",
        "Example 11's system with the rate law moved out of C++ and into\n"
        "input/rates.sym, read at start-up.",
        "- growth of the same order as example 11: same organism, same two\n"
        "  substrates\n"
        "- **Fe(III) drawn down 4.8 times faster than acetate, in molar terms.**\n"
        "  That ratio is the reaction stoichiometry, and nobody tells the solver\n"
        "  it at run time -- it falls out of the three lines in rates.sym. If it\n"
        "  comes back at 4.8, the file parsed and every name bound correctly\n"
        "- the whole rate law echoed into the log before iteration 0",
        "`<symbolic>` with `<expressions_file>`, and\n"
        "`<reaction_type>symbolic</reaction_type>` on the organism.",
        "Change 0.0527 in rates.sym and rerun without rebuilding: that is the\n"
        "capability this case exists to show. Push a concentration outside a\n"
        "`range` line and the end-of-run report tells you what fraction of\n"
        "evaluations were clamped, and warns above 5%.\n"
        "\n"
        "`fitting/` holds the worked material behind this: the sweeps, the fitted\n"
        "files, and a checker that holds one of them against the law it came from."),
        "one_microbe.dat", extra_files={"input/rates.sym": RATES_SYM})

    print("\n16 examples written.")


if __name__ == "__main__":
    build()
