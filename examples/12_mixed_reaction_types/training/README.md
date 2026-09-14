# Preparing a run where organisms use different rate paths

**In:** a genome-scale model for the organism that needs one, and a kinetics
header for the organism that does not.
**Out:** a converted model the solver reads at start-up, and a compiled kinetics
block.
**Run:** `MODEL=yourmodel.xml OBJECTIVE=<biomass reaction id> ./offline.sh`
**Needed for:** this case only if you bring your own model. It ships everything
it needs.

---

## What this case is for

Every other case gives all its organisms the same kind of rate law. This one
does not. One organism has its metabolism solved as a linear program; another
has a saturating law compiled in. They share the same pore space, the same
substrates and the same time step.

```
    organism 1   <reaction_type>fba</reaction_type>        needs a model file
    organism 2   <reaction_type>kinetics</reaction_type>   needs a kinetics block
                              |
                              v
              one substrate budget per voxel, settled once
```

The preparation therefore has **two halves that have nothing to do with each
other**, which is the thing to hold on to when something goes wrong: a failure
is nearly always in one half or the other, not in the coupling.

---

## Half one: the model, for the flux balance organism

```bash
MODEL=yourmodel.xml OBJECTIVE=Biomass_Ecoli_core ./offline.sh
```

which runs `training/extractMM.py` and prints the exchange reaction table.
Identical to example 09; see that case's training README for what the tool does
and why to bind exchanges by name rather than index.

Then in `CompLaB.xml`, under that organism's block only:

```xml
<model_filename>input/converted_model</model_filename>
<exchange_reaction_names>EX_ac_e EX_o2_e</exchange_reaction_names>
```

## Half two: the kinetics, for the compiled organism

```bash
python3 training/makeKinetics.py --help
```

writes the C++ block that gets compiled in. Nothing is fitted here either: you
state the law and its constants, and the tool turns them into code so the rate
is evaluated without a file read or a parser in the inner loop.

Changing those constants means **rebuilding**. That is the trade this path
makes, and it is why examples 17 and 18 exist: a law in a file can be changed
without touching the compiler.

---

## What to check before you run

**The substrate names have to agree across both halves.** The exchange reactions
you bind for the flux balance organism and the names the kinetics block uses
must both resolve to entries in `<name_of_substrates>`. They are the only thing
connecting the two organisms to the same chemistry.

**Both organisms must be listed in `<name_of_microbes>`,** each with its own
`<reaction_type>`. An organism with no reaction type stops the run at start-up,
which is the behaviour you want: silently giving it a default would have it grow
on a law nobody chose.

**Read the start-up log.** It prints one line per organism naming the path it
was given. Two organisms both showing the same path means an edit went into the
wrong block, and that is far easier to see there than in the output fields.

---

## The one that catches people

**They share one substrate budget per voxel, and it is settled once.** If both
organisms want more donor than is present, neither gets what it asked for. That
is correct, and it means a rate you compute by hand for one organism in
isolation will not match what the run reports once the other is competing for
the same substrate. The closing report gives the drawn amounts per species, and
that is the number to compare against, not the rate law's own output.

---

## Related

- Example 09, the flux balance half on its own.
- Examples 05 to 07, the compiled kinetics half on its own.
- Example 08, two organisms both on compiled kinetics, which is the simpler case
  to debug first if the competition is behaving oddly.
