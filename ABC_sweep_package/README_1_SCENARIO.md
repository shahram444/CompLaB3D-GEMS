# The scenario: A + B makes C, in a pore space, with nothing in it to start

This is the reactive transport problem the sweep simulates, with every number in
it and the reason for each. The second README is the instructions; this one is
what the instructions are for.

---

## 1. In one paragraph

A pore space 320 micrometres across is full of water and nothing else. From one
end we push in a solution of **A**. From the other end, against the flow, we hold
a solution of **B**. Where the two meet they react, irreversibly and in one step,
to make **C**:

```
A + B  ->  C          R = k [A][B]        mol/L/s
```

Nothing happens where only one of them is present, because the rate is a product
and a product with a zero in it is zero. So the reaction lives on the surface
where the two plumes overlap, and **where that surface sits, and what shape it
has, is set by the flow and by the pore structure together.** That surface is
what the network is being asked to learn.

Two dimensionless numbers move it. Péclet decides how far B can work upstream
before the flow sweeps it back, which slides the meeting surface from the middle
of the domain towards the outlet. Damköhler decides how sharp the surface is:
slow reaction and the two plumes interpenetrate over most of the domain, fast
reaction and they annihilate on contact in a thin sheet.

---

## 2. The chemistry

One reaction, abiotic, volumetric, in every pore voxel. No microbes, no minerals,
no equilibrium speciation, and a pore space that never changes shape.

| species | index | role | initial | fed from |
|---|---|---|---|---|
| A | 0 | reactant | **0.0** | x = 0, the inlet |
| B | 1 | reactant | **0.0** | x = nx-1, the outlet |
| C | 2 | product | **0.0** | nowhere; it is made |

```
R      = k [A][B]           mol/L/s
dA/dt  = -R    dB/dt = -R    dC/dt = +R
```

Written in `kinetics/defineAbioticKinetics.hh`. `k` is read from the environment
(`PRT_KABIO`) so one binary serves the whole sweep with no rebuild between cases.

### Why second order and not a first-order decay

A first-order sink `R = -k[A]` is separable. Along a streamline the answer is an
exponential, the second species is decoration, and a network can fit it from the
Péclet number alone without ever looking at the pore space. There is nothing to
transfer and nothing a geometry-aware model could be better at.

`R = k[A][B]` cannot be solved one species at a time. The rate is zero wherever
either reactant is absent, however much of the other is present, so the reaction
is confined to the overlap, and the overlap is a geometric object. Change the
pore space and it moves.

### Why nothing starts in the domain

Every `<initial_concentration>` is `0.0`. That was the requirement, and it has one
consequence worth stating plainly: **both reactants have to enter through a
boundary, and CompLaB3D offers Dirichlet boundaries on the two x-normal planes
and nowhere else.** Two species, two available faces. So they come in from
opposite ends and travel toward each other. This is a counter-current reactor,
and it is the only configuration the constraint allows.

It is also the better problem. In a co-current feed, A and B entering together
with the same diffusivity satisfy identical equations with identical boundary
conditions and are therefore the same field; B carries no information and the
problem collapses to one species reacting with itself. Counter-current is the
version where both species matter.

---

## 3. The domain

| | |
|---|---|
| grid | 32 x 32 x 32 voxels |
| voxel | dx = 10 micrometres |
| domain | L = 320 micrometres |
| flow lattice | D3Q19, tau = 0.8 |
| transport lattice | D3Q7, one per species |
| diffusivity | D = 1e-9 m2/s, **the same for A, B and C** |
| timestep | 0.0086 to 0.0090 s, measured per case, not assumed |

`<nx>32</nx>` means the geometry file holds 32 x-slices: `complab_functions.hh`
does `nx += 2` after reading the tag and `readGeometry()` fills x = 1..nx-2 from
the file, duplicating the first and last slice into two ghost columns. The
transport domain is 32 voxels long.

**One diffusivity for all three species, on purpose.** It makes three combinations
of the fields exactly reaction-free (section 9), which is what
`tools/analytic_channel.py` checks against closed form. It also means no species
outruns another, so the meeting surface is set by the flow and the pore structure
alone rather than by a diffusivity contrast.

### Why the timestep is measured rather than computed

`complab.cpp:1276` sets `refNu = PoreMeanU * characteristic_length / Pe`, where
`PoreMeanU` is the mean speed of the flow field **in that geometry**, and then
`ade_dt = refNu * dx^2 / D`. So the timestep depends on the pore space and does
not exist until the flow has been solved. `make_campaign.py build --complab PATH`
runs the solver once per geometry and Péclet for a single step, reads the
`[ADE] dt=` line it prints, and sizes every case from that. Without it the run
lengths are an estimate and the tool says so everywhere a reader might look.

---

## 4. The boundary conditions

Every face is **held**. That is what the sweep runs, and it is the default:
`make_campaign.py build --boundaries dirichlet`.

| face | A | B | C |
|---|---|---|---|
| x = 0 | Dirichlet **1.0** | Dirichlet 0.0 | Dirichlet 0.0 |
| x = nx-1 | Dirichlet 0.0 | Dirichlet **1.0** | Dirichlet 0.0 |
| y = 0, y = ny-1 | wall, drawn into the geometry | | |
| z = 0, z = nz-1 | wall, drawn into the geometry | | |

Read it as two well-mixed reservoirs with the rock between them: one reservoir
holds A feed, the other holds B feed, and each species enters from its own end and
works toward the other. The zeros are not artificial. "No A in the B reservoir" is
true by construction, and the same for B in the A reservoir and for C in both,
since C is a product and neither reservoir is stirring any of it back in. Nothing
can be manufactured at a held plane, and `A - B`, `A + C` and `B + C` then have
closed-form steady profiles, which is what `tools/analytic_channel.py` checks.

### Why not the conventional open outflow

The ordinary flow-through arrangement holds each species at the face it is fed
from and leaves it open at the face it leaves by. That reads as the obvious choice
here, and `make_campaign.py build --boundaries mixed` still builds it:

| face | A | B | C |
|---|---|---|---|
| x = 0, inlet | Dirichlet **1.0** | Neumann | Neumann |
| x = nx-1, outlet | Neumann | Dirichlet **1.0** | Neumann |

**It was measured on this problem and it fails.** A full 60-run campaign was built
this way, ran to completion, and passed every operational check: 60 of 60 exit
clean, the Damköhler banner exact, the achieved Péclet within 0.3 per cent of
target, the positivity clamp at 0.000 per cent, no negative-concentration flags.
The concentrations are still wrong.

Nothing in this problem can concentrate. No species is fed above 1.0 mol/L and the
reaction only consumes, so any value above the feed was manufactured. Taking the
last frame of every run and masking to pore voxels:

| Pe | runs | median worst value | worst value | runs over 1.05 x feed |
|---|---|---|---|---|
| 0.1 | 20 | 1.001 | 1.018 | 0 |
| 0.5 | 20 | 1.014 | 1.165 | 4 |
| 2.0 | 20 | 1.462 | **67.3** | 20 |

The mechanism is visible in the worst case. In run_0010 (Pe 2.0, Da 0.1) B reaches
**67.3 on plane x = 0** and decays inward 67.3, 48.1, 29.9, 11.2. Plane x = 0 is
B's Neumann face **and** the face the flow enters through. A zero-gradient
condition at an inflow is ill-posed: `FlatAdiabaticBoundaryFunctional3D` implements
Neumann by calling `defineDensity()` on the plane with the density of the layer
inside it, the plane's value is then advected inward, and the plane re-reads it
from the layer it just fed. At Pe 2.0 the gain around that loop exceeds one and the
face runs away.

An open face at a genuine **outflow** is not safe either. In run_0004 (Pe 2.0,
Da 0.1) species A, a near-passive tracer at that Damköhler, has a pore mean of
**1.051** and **52 per cent of its pore voxels above 1.05 x feed**, with its peak
in the interior at x = 26. `defineDensity()` discards whatever the plane was
holding and replaces it with the interior value, which for a field still rising
from zero tops the plane up from nothing every step and streams the surplus back
in. The same failure was measured at 32 per cent of the water-phase mass on
example 14 of the CompLaB distribution, where calcium starts at zero and can only
come from calcite. Every field here rises from zero.

Use `--boundaries mixed` only knowing all of that. The source carries the same
numbers in a comment above `BOUNDARY_SETS`.

### The check that catches it

`tools/verify_dataset.py` reports the largest value anywhere in the file against
the feed. Over 1.05 x feed it warns. Over 1.25 x feed it **fails the dataset**,
because at that point the field is not a noisy version of the right answer, it is a
different answer. The one clean example from the failed campaign shows what the
gate is protecting: run_0001 (Pe 0.1, Da 10) has A falling 1.003 to 0.168 across
x, B rising 0.051 to 0.997 the other way, C peaking at 0.632 around x = 26, and
zero voxels over 1.05 x feed. That is the designed physics, and it is what every
run should look like.

### Why the y and z walls are drawn into the geometry file

Those four faces get no boundary condition at all. The advection-diffusion
lattices stream off the block wherever they are open pore and read back an
envelope nothing updates. `tools/make_geometries.py` writes bounce-back on all
four in every file, which is what the solver's own note recommends over
`<outer_faces>sealed</outer_faces>`: the wall is then part of the geometry the
network is shown, visible in the file and in the preview, rather than something
done to the geometry afterwards.

The pore space is therefore a duct of 30 x 30 open voxels packed with grains, and
**porosity is quoted over that interior**, because the interior is the packing and
the walls are the container.

---

## 5. The two dimensionless numbers

```
Pe   = U L / D                 how fast the flow carries, against how fast diffusion spreads
Da_r = k_r L^2 / D             how fast the reaction runs, against how fast diffusion feeds it
```

both built on the same L, which is what makes them comparable.

**Damköhler needs a first-order rate constant**, because `k L^2 / D` is
dimensionless only for one; a second-order `k` has units that do not cancel. The
pseudo-first-order collapse is `k_r = R(C_ref)/C_ref` at a reference composition,
and the reference here is the **feed**, because nothing starts in the domain and
the feed is the only composition the reaction ever sees at full strength:

```
R(A0, A0) = k A0^2        k_r = k A0        Da = k A0 L^2 / D
```

so a campaign that wants a given Da sets

```
k = Da * D / (A0 * L^2)
```

which is exactly what `make_campaign.py` computes and exactly what the startup
banner reports back, once `patch/apply_damkohler_feed_reference.py` is applied.
Stock v1.3.2 builds its reference from `<initial_concentration>`, which is zero in
every case here, so without the patch it reports nothing at all for the whole
sweep.

### What the sweep covers

```
Pe  in  {0.02, 0.2, 2.0}         3 values, a decade apart
Da  in  {0.1, 1.0, 10}           3 values, a decade apart
geometries                      30   (10 porosities x 3 packings)
                             ------
runs                           270
snapshots each                  21   (one at t = 0, which the collector drops)
training pairs                5400
```

Both ladders are evenly spaced in the log of the quantity, because the log is the
axis the surrogate is conditioned on. Pe stops at 0.02 rather than 0 because
`complab.cpp:1156` switches the flow solver off entirely at Pe = 0: there would be
no velocity field to feed the surrogate and no mean pore speed for the timestep to
be derived from. At 0.02 advection contributes about two per cent of the transport,
which is the diffusive limit in every way that matters, with the machinery
unchanged.

### What Péclet does, physically

B has to work upstream. Its steady profile is
`B0 (e^{Pe x/L} - 1)/(e^{Pe} - 1)`, which falls by a factor of e over a distance
`L/Pe`. Solving A and B as conservative tracers and asking where they cross gives
the meeting surface:

| Pe | meeting surface at porosity 0.32 | at porosity 0.59 | B's reach from its feed face |
|---:|---:|---:|---|
| 0.02 | 0.508 | 0.504 | the whole domain, many times over |
| 0.2 | 0.577 | 0.542 | about 51 voxels, still longer than the domain |
| 2.0 | 0.889 | 0.805 | about 5 voxels |

Two things to read off that. Down a column is the Péclet effect the surrogate is
conditioned on. **Across each row is the porosity effect at fixed Péclet**, and it
grows from under a voxel at Pe = 0.02 to about 2.7 voxels at Pe = 2.0, so the
geometry dependence is strongest at the top of the ladder. That is where the
held-out porosities are really being tested.

The reason porosity moves the front at all is that CompLaB's Péclet is built on the
mean speed over the whole box, a superficial speed, so the speed a solute actually
feels is larger by about 1/porosity.

**The sweep stops at Pe = 2, and that is a resolution limit rather than a
preference.** B's e-folding reach upstream is `nx * porosity / Pe` voxels. At
Pe = 2 that is about 5 voxels in the tightest pack: thin but resolved. At Pe = 20
it is 0.51 voxels, the cell Péclet reaches 1.95 where the D3Q7 lattice begins to
oscillate and the positivity clamp starts firing, and the whole reaction collapses
onto the feed plane where the pore structure has nothing left to say about it.
Fixing that by refinement would need nx near 190 rather than 32. A parameter value
that destroys the thing being learned is not a wider sweep, it is a wasted third of
one.

### What Damköhler does

| Da | k, L/(mol s) | regime | what the field looks like |
|---:|---:|---|---|
| 0.1 | 9.766e-4 | transport limited | the two plumes interpenetrate over most of the domain and C is spread out and weak |
| 1.0 | 9.766e-3 | balanced | reaction and transport times comparable, and the least predictable of the three |
| 10 | 9.766e-2 | reaction limited | A and B annihilate where they meet, C is a sheet, and each reactant is confined to its own side |

**Three values rather than two.** The previous campaign used 0.1 and 10 only: two
decades in a single jump, one interval, and no way for a surrogate conditioned on
log Da to show whether it had learned a trend or memorised two labels. 1.0 is the
transition point where neither the reaction nor the transport is the bottleneck,
which is also the hardest of the three to predict.

---

## 6. The concentration unit

The feed is **1.0 mol/L** for both reactants. That is high for a real groundwater
and it is chosen deliberately, for one reason that is worth being explicit about:

`/samples/conc` in the dataset is **float16**, which is what
`3D/tools/dataset_reader.py` reads and what every PRT-DeepONet dataset holds.
float16's smallest normal value is 6.1e-5. A feed of 1e-3 mol/L would put the
leading edge of the front, the part at a thousandth of the feed and below, into
the subnormal range where the relative error is several per cent and rising. That
edge is exactly where the interesting physics is.

A feed of 1.0 puts five decades of the field inside the normal range, and it makes
the stored number its own dimensionless concentration: `C/A0 = C`. The chemistry
is a model chemistry and Pe and Da are what define the problem, so nothing
physical is given up by choosing the unit this way.
`tools/verify_dataset.py` checks the resulting headroom rather than assuming it.

---

## 7. The pore spaces

**One family, ten porosities, three packings each.** Thirty packs of overlapping
spheres, grain radius 3 voxels:

| target porosity | packings | achieved | tortuosity | gids | split |
|---:|---:|---|---|---|---|
| 0.32 | 3 | 0.3225, 0.3227, 0.3208 | 1.406, 1.376, 1.366 | 0, 1, 2 | train |
| 0.35 | 3 | 0.3523, 0.3497, 0.3516 | 1.416, 1.433, 1.325 | 3, 4, 5 | **validation** |
| 0.38 | 3 | 0.3821, 0.3785, 0.3790 | 1.325, 1.298, 1.277 | 6, 7, 8 | train |
| 0.41 | 3 | 0.4120, 0.4105, 0.4110 | 1.336, 1.259, 1.236 | 9, 10, 11 | **test** |
| 0.44 | 3 | 0.4416, 0.4395, 0.4415 | 1.300, 1.264, 1.322 | 12, 13, 14 | train |
| 0.47 | 3 | 0.4718, 0.4685, 0.4685 | 1.273, 1.245, 1.215 | 15, 16, 17 | train |
| 0.50 | 3 | 0.5001, 0.4992, 0.5039 | 1.206, 1.244, 1.234 | 18, 19, 20 | **validation** |
| 0.53 | 3 | 0.5319, 0.5331, 0.5313 | 1.188, 1.202, 1.181 | 21, 22, 23 | train |
| 0.56 | 3 | 0.5609, 0.5591, 0.5604 | 1.155, 1.161, 1.157 | 24, 25, 26 | **test** |
| 0.59 | 3 | 0.5900, 0.5917, 0.5909 | 1.166, 1.140, 1.159 | 27, 28, 29 | train |

plus one **straight duct** (gid 900), run separately, which is the case with
closed-form answers and is not training data.

One family rather than a mixture of families, because a mixture at this size gives
three or four examples of each structure, and three examples of a structure teach a
network that the structure is a label rather than a variable. Ten porosities of one
family is a dimension the network can actually move along.

**Three packings per level rather than two.** Two cannot separate a porosity effect
from one packing's quirks: if the two happen to agree there is nothing to compare
them against, and if they disagree there is no way to tell which is the outlier.
Three gives a spread at fixed porosity, which is the only thing that says how much
of a held-out error is porosity and how much is packing.

**The split is held out by porosity, not drawn at random.** Two packings at the same
porosity are close enough that splitting between them would leak: the test score
would measure how well the network recognised a porosity it had already been
trained on. Every held-out level here sits strictly inside the training range
[0.32, 0.59]:

```
validation   0.35   between training 0.32 and 0.38
             0.50   between training 0.47 and 0.53
test         0.41   between training 0.38 and 0.44
             0.56   between training 0.53 and 0.59
```

so this asks for interpolation, which is the honest question for an operator
network. Extrapolating below 0.32 or above 0.59 is a different experiment and
deserves its own geometries.

**Two held-out levels per split rather than one.** With a single held-out porosity
there is no way to tell a real generalisation gap from that one level happening to
be easy. Two disagree or they do not, and either answer is informative.

Each geometry is delivered as CompLaB's material map (0 grain interior, 1
bounce-back shell, 2 pore) plus the two distance fields PRT-DeepONet's trunk uses:
the geodesic distance from the inlet through the pore space, and the Euclidean
distance to the nearest solid. Achieved tortuosity runs from 1.14 at the open end
to 1.43 at the tight end, decreasing with porosity across the whole family.

---

## 8. The 270 runs, and which is which

**270 runs. 162 train, 54 validation, 54 test.** Each run writes 21 snapshots, the
collector drops the one at iteration 0 because it is identically zero, so each run
contributes **20 training pairs**:

| | porosity levels | geometries | runs | training pairs | share |
|---|---|---:|---:|---:|---:|
| **train** | 0.32, 0.38, 0.44, 0.47, 0.53, 0.59 | 18 | 162 | 3240 | 60 per cent |
| **validation** | 0.35, 0.50 | 6 | 54 | 1080 | 20 per cent |
| **test** | 0.41, 0.56 | 6 | 54 | 1080 | 20 per cent |
| total | 10 | **30** | **270** | **5400** | |

plus the straight duct, gid 900, which is run separately in its own campaign and is
**not** part of the dataset. It is the closed-form check, not training data.

### The full matrix

Every geometry gets the same nine parameter pairs, in the same order, so a run id
tells you everything about the case:

```
run id = 9 * gid + (3 * Peclet index) + Damkohler index
```

| gid | porosity | split | runs | Pe = 0.02 | Pe = 0.2 | Pe = 2.0 |
|---:|---:|---|---|---|---|---|
| | | | | Da 0.1 / 1 / 10 | Da 0.1 / 1 / 10 | Da 0.1 / 1 / 10 |
| 0 | 0.3225 | train | 0000-0008 | 0000 / 0001 / 0002 | 0003 / 0004 / 0005 | 0006 / 0007 / 0008 |
| 1 | 0.3227 | train | 0009-0017 | 0009 / 0010 / 0011 | 0012 / 0013 / 0014 | 0015 / 0016 / 0017 |
| 2 | 0.3208 | train | 0018-0026 | 0018 / 0019 / 0020 | 0021 / 0022 / 0023 | 0024 / 0025 / 0026 |
| 3 | 0.3523 | **validation** | 0027-0035 | 0027 / 0028 / 0029 | 0030 / 0031 / 0032 | 0033 / 0034 / 0035 |
| 4 | 0.3497 | **validation** | 0036-0044 | 0036 / 0037 / 0038 | 0039 / 0040 / 0041 | 0042 / 0043 / 0044 |
| 5 | 0.3516 | **validation** | 0045-0053 | 0045 / 0046 / 0047 | 0048 / 0049 / 0050 | 0051 / 0052 / 0053 |
| 6 | 0.3821 | train | 0054-0062 | 0054 / 0055 / 0056 | 0057 / 0058 / 0059 | 0060 / 0061 / 0062 |
| 7 | 0.3785 | train | 0063-0071 | 0063 / 0064 / 0065 | 0066 / 0067 / 0068 | 0069 / 0070 / 0071 |
| 8 | 0.3790 | train | 0072-0080 | 0072 / 0073 / 0074 | 0075 / 0076 / 0077 | 0078 / 0079 / 0080 |
| 9 | 0.4120 | **test** | 0081-0089 | 0081 / 0082 / 0083 | 0084 / 0085 / 0086 | 0087 / 0088 / 0089 |
| 10 | 0.4105 | **test** | 0090-0098 | 0090 / 0091 / 0092 | 0093 / 0094 / 0095 | 0096 / 0097 / 0098 |
| 11 | 0.4110 | **test** | 0099-0107 | 0099 / 0100 / 0101 | 0102 / 0103 / 0104 | 0105 / 0106 / 0107 |
| 12 | 0.4416 | train | 0108-0116 | 0108 / 0109 / 0110 | 0111 / 0112 / 0113 | 0114 / 0115 / 0116 |
| 13 | 0.4395 | train | 0117-0125 | 0117 / 0118 / 0119 | 0120 / 0121 / 0122 | 0123 / 0124 / 0125 |
| 14 | 0.4415 | train | 0126-0134 | 0126 / 0127 / 0128 | 0129 / 0130 / 0131 | 0132 / 0133 / 0134 |
| 15 | 0.4718 | train | 0135-0143 | 0135 / 0136 / 0137 | 0138 / 0139 / 0140 | 0141 / 0142 / 0143 |
| 16 | 0.4685 | train | 0144-0152 | 0144 / 0145 / 0146 | 0147 / 0148 / 0149 | 0150 / 0151 / 0152 |
| 17 | 0.4685 | train | 0153-0161 | 0153 / 0154 / 0155 | 0156 / 0157 / 0158 | 0159 / 0160 / 0161 |
| 18 | 0.5001 | **validation** | 0162-0170 | 0162 / 0163 / 0164 | 0165 / 0166 / 0167 | 0168 / 0169 / 0170 |
| 19 | 0.4992 | **validation** | 0171-0179 | 0171 / 0172 / 0173 | 0174 / 0175 / 0176 | 0177 / 0178 / 0179 |
| 20 | 0.5039 | **validation** | 0180-0188 | 0180 / 0181 / 0182 | 0183 / 0184 / 0185 | 0186 / 0187 / 0188 |
| 21 | 0.5319 | train | 0189-0197 | 0189 / 0190 / 0191 | 0192 / 0193 / 0194 | 0195 / 0196 / 0197 |
| 22 | 0.5331 | train | 0198-0206 | 0198 / 0199 / 0200 | 0201 / 0202 / 0203 | 0204 / 0205 / 0206 |
| 23 | 0.5313 | train | 0207-0215 | 0207 / 0208 / 0209 | 0210 / 0211 / 0212 | 0213 / 0214 / 0215 |
| 24 | 0.5609 | **test** | 0216-0224 | 0216 / 0217 / 0218 | 0219 / 0220 / 0221 | 0222 / 0223 / 0224 |
| 25 | 0.5591 | **test** | 0225-0233 | 0225 / 0226 / 0227 | 0228 / 0229 / 0230 | 0231 / 0232 / 0233 |
| 26 | 0.5604 | **test** | 0234-0242 | 0234 / 0235 / 0236 | 0237 / 0238 / 0239 | 0240 / 0241 / 0242 |
| 27 | 0.5900 | train | 0243-0251 | 0243 / 0244 / 0245 | 0246 / 0247 / 0248 | 0249 / 0250 / 0251 |
| 28 | 0.5917 | train | 0252-0260 | 0252 / 0253 / 0254 | 0255 / 0256 / 0257 | 0258 / 0259 / 0260 |
| 29 | 0.5909 | train | 0261-0269 | 0261 / 0262 / 0263 | 0264 / 0265 / 0266 | 0267 / 0268 / 0269 |

The porosities are what the generator actually achieved, not the targets; they are
in `geometries/index.csv` and in every `params.json`.

### Who decides the split, and where it is written

`tools/make_geometries.py` assigns it from the porosity and writes it into
`geometries/index.csv`. `make_campaign.py` copies it into each case's
`params.json`. `collect_to_h5.py` writes it into the dataset as
`/samples/split`, and `tools/train_prt3d.py` reads it from there, so the split
travels with the data and nothing downstream has to re-derive it.

**Nothing is split by sample or by run.** Whole geometries are held out. Two
packings at the same porosity are close enough to each other that splitting
between them would leak the pore structure into the test set, and the score would
then measure how well the network recognised a porosity it had already trained
on. `tools/verify_dataset.py` checks that the three splits share no geometry.

A note if you use the PRT-DeepONet-GUI's own trainer, `3D/model/train.py`: it
draws its test split at random by geometry with `--test-frac` rather than reading
`/samples/split`. That is still a geometry split and still honest; it is just not
this one. `tools/train_prt3d.py` uses the porosity split above.

---

## 9. What can be checked against a closed form

Three combinations of the three fields obey a **reaction-free** equation exactly,
for any rate constant at all. The reaction is one-for-one-for-one and all three
species share one diffusivity, so

| combination | why the reaction cancels |
|---|---|
| A - B | by subtraction: both lose R |
| A + C | by addition: A loses R, C gains it |
| B + C | likewise |

Each is therefore a conservative tracer. **Carrying WHICH boundary values depends
on the variant**: with every face held, each combination's boundary values are
known in advance, and at steady state in a duct each equals the plain
advection-diffusion solution

```
S(x) = S(0) + (S(L) - S(0)) * (e^{Pe x/L} - 1) / (e^{Pe} - 1)
```

which at Pe = 0 is a straight line. Nothing about the rate law enters, so a
disagreement is a **transport** fault, and these would be broken by an inequality
anywhere in the stoichiometry, by the three species not sharing a diffusivity, or
by the positivity clamp cutting one species and not its partner.

With an open outflow the boundary value is instead whatever the flow leaves
there, so there is nothing to compare against. `tools/analytic_channel.py` refuses
a case built any other way rather than reporting a meaningless number. Since the
sweep is all Dirichlet too, the duct case and the sweep share a boundary set and
the duct result speaks directly to the campaign.

It is the only part of the package that checks the physics rather than the
plumbing.

---

## 10. Run length and snapshots

```
t_end = 1.25 * (L^2 / D) / (Pe + 1)
```

One crossing of the domain, then a quarter again. The two terms are the diffusive
time `L^2/D = 102.4 s` and the advective time `L^2/(Pe D)`; combining them as a
harmonic sum is what makes the expression right at both ends and continuous
between.

It is deliberately conservative, because CompLaB's Péclet is superficial and the
interstitial front arrives sooner, by roughly the porosity. And it is
geometry-independent, which is the point: **every case at one Péclet covers the
same physical interval**, so `t_norm = 0.5` means one thing across the whole
dataset.

| Pe | one crossing | run to | snapshot every | steps at the estimate | steps at the measured dt |
|---:|---:|---:|---:|---:|---:|
| 0.02 | 100.4 s | 125.6 s | 6.3 s | 10921 | about 14400 |
| 0.2 | 85.4 s | 106.7 s | 5.3 s | 9281 | about 12300 |
| 2.0 | 34.2 s | 42.8 s | 2.1 s | 3721 | about 4900 |

The last column is the honest one: the estimator's timestep of 0.0115 s is about
24 per cent larger than the 0.0087 s the solver actually reported on the previous
campaign, so the calibrated step counts come out that much higher. The build step
measures it rather than guessing, and the numbers it prints replace both.

21 snapshots, evenly spaced, the first at iteration 0. `ade_max_iT` is set to
`interval * 20 + 1` so the last write lands on the last step and every run holds
the same number of frames. `<ade_converge_iT>` is read by the solver and never
used, so no run stops early and the time axes cannot drift apart.

The collector drops the frame at iteration 0: nothing starts in the domain, so it
is identically zero in every case, and keeping it would spend a twenty-first of
the dataset on a sample whose answer is "all zero" for every geometry and every
parameter. That leaves 20 usable snapshots per run and 5400 training pairs.

**Cost.** 270 runs, between 4900 and 14400 transport steps each on 32768 voxels
with three transport lattices. That is **23 to 31 core-hours in total**, the two
figures bracketing the answer at the estimated and the measured timestep. With the
Slurm array throttled to 12 concurrent tasks it finishes in **2 to 3 hours of wall
clock**, and no single run takes more than about 15 minutes, so the 2 hour limit
per task in the documented command has a wide margin.

About **21 GB of output**: the concentration volumes, a rate volume beside each
one, the flow field and the closing checkpoints, at roughly 77 MB per run. The
collected `.h5` is about **1.0 GB**, most of it the float16 concentration array at
(270, 20, 3, 32, 32, 32).

---

## 11. Every input value in CompLaB.xml, annotated

This is `campaign02/runs/run_0000/CompLaB.xml` with the long reasoning removed and
one short comment per value. The shipped `xml/CompLaB.xml.template` is the same
file with the reasoning left in; `make_campaign.py` fills the `@TOKEN@`
placeholders. **Do not edit a rendered copy** or the case and the `params.json`
beside it stop agreeing.

```xml
<parameters>

  <path>
    <src_path>src</src_path>          <!-- where the COBRApy bridge would be; unused here -->
    <input_path>input</input_path>    <!-- geometry.dat lives here -->
    <output_path>output</output_path> <!-- .vti, run.log and summary.csv land here -->
  </path>

  <simulation_mode>
    <biotic_mode>false</biotic_mode>              <!-- no microbes; skips <microbiology> entirely -->
    <enable_kinetics>false</enable_kinetics>      <!-- defineKinetics.hh never called -->
    <enable_abiotic_kinetics>true</enable_abiotic_kinetics>  <!-- defineAbioticKinetics.hh: A + B to C -->
    <abiotic_rate_scale>1.0</abiotic_rate_scale>  <!-- PINNED. The rate comes from PRT_KABIO;
                                                       both would multiply and double it -->
    <enable_validation_diagnostics>false</enable_validation_diagnostics>  <!-- per-voxel tracing, very slow -->
  </simulation_mode>

  <LB_numerics>
    <domain>
      <nx>32</nx>                     <!-- x-SLICES IN THE FILE. The solver adds 2 ghost columns -->
      <ny>32</ny>
      <nz>32</nz>
      <dx>10</dx>                     <!-- voxel size, in <unit> -->
      <unit>um</unit>                 <!-- so dx = 10 micrometres and the domain is 320 -->
      <characteristic_length>320</characteristic_length>
                                      <!-- L, same unit as dx. Divided by dx at parse time, so it
                                           becomes 32 lattice units. The L in BOTH Pe and Da -->
      <filename>geometry.dat</filename>   <!-- read from <input_path>; 32*32*32 integers, no header -->
      <material_numbers>
        <pore>2</pore>                <!-- open water; the only place chemistry happens -->
        <solid>0</solid>              <!-- grain interior; NO DYNAMICS, not simulated at all -->
        <bounce_back>1</bounce_back>  <!-- grain surface and the confining wall -->
                                      <!-- NO <microbe0>: there are no microbes -->
      </material_numbers>
      <outer_faces>open</outer_faces> <!-- leave the y/z faces as written: the generator already
                                           drew bounce-back there, so there is nothing to seal -->
    </domain>

    <delta_P>1e-05</delta_P>          <!-- A SEED ONLY. The solver measures permeability and
                                           back-solves the real drop to hit <Peclet>. But a value
                                           below threshold silently sets Pe = 0 and the case
                                           becomes pure diffusion, so it must not be zero -->
    <Peclet>0.02</Peclet>             <!-- THE TARGET. 0.02, 0.2 or 2.0 across the sweep -->
    <tau>0.8</tau>                    <!-- flow relaxation time; nu = (tau - 1/2)/3 -->
    <track_performance>false</track_performance>   <!-- true would time the kernels and write no VTI -->

    <iteration>
      <ns_max_iT1>60000</ns_max_iT1>  <!-- cap on the first flow solve -->
      <ns_max_iT2>60000</ns_max_iT2>  <!-- cap on any re-solve; there are none, the geometry is frozen -->
      <ns_converge_iT1>1e-8</ns_converge_iT1>   <!-- flow convergence tolerance, first solve -->
      <ns_converge_iT2>1e-6</ns_converge_iT2>
      <ns_update_interval>100000000</ns_update_interval>
                                      <!-- larger than any run, so the flow is solved once -->
      <ade_max_iT>10121</ade_max_iT>  <!-- transport steps. = save_VTK_interval * 20 + 1, so the
                                           last write lands on the last step and every run holds
                                           exactly 21 frames. Set per case from Pe and the
                                           MEASURED timestep -->
      <ade_converge_iT>0</ade_converge_iT>
                                      <!-- read by the solver and never used, so no run stops
                                           early and the time axes cannot drift apart -->
      <ade_update_interval>1</ade_update_interval>   <!-- rebuild diffusivities every step -->
    </iteration>
  </LB_numerics>

  <chemistry>
    <number_of_substrates>3</number_of_substrates>   <!-- A, B, C -->

    <!-- SUBSTRATE 0 IS SPECIAL: its <in_pore> sets the global timestep
         ade_dt = refNu * dx^2 / D[0], and every other species is normalised by it -->
    <substrate0>
      <name_of_substrates>A</name_of_substrates>     <!-- ALSO the .vti base name: A_0000000.vti -->
      <initial_concentration>0.0</initial_concentration>   <!-- nothing starts in the domain -->
      <substrate_diffusion_coefficients>
        <in_pore>1e-09</in_pore>      <!-- m2/s. Same for all three, which is what makes
                                           A-B, A+C and B+C exactly reaction-free -->
        <in_biofilm>1e-09</in_biofilm><!-- no biofilm here; still parsed, so still set -->
      </substrate_diffusion_coefficients>
      <left_boundary_type>Dirichlet</left_boundary_type>    <!-- x = 0, the A reservoir: A is fed here -->
      <left_boundary_condition>1</left_boundary_condition>  <!-- A0 = 1.0 mol/L, the feed -->
      <right_boundary_type>Dirichlet</right_boundary_type>  <!-- x = nx-1, the B reservoir: no A in it,
                                                                 true by construction -->
      <right_boundary_condition>0.0</right_boundary_condition>
    </substrate0>

    <substrate1>
      <name_of_substrates>B</name_of_substrates>
      <initial_concentration>0.0</initial_concentration>
      <substrate_diffusion_coefficients>
        <in_pore>1e-09</in_pore>
        <in_biofilm>1e-09</in_biofilm>
      </substrate_diffusion_coefficients>
      <left_boundary_type>Dirichlet</left_boundary_type>       <!-- x = 0, the A reservoir: no B in it -->
      <left_boundary_condition>0.0</left_boundary_condition>
      <right_boundary_type>Dirichlet</right_boundary_type>     <!-- B's FEED face, the far end -->
      <right_boundary_condition>1</right_boundary_condition>   <!-- B0 = 1.0, fed at x = nx-1, so B
                                                                    works back against the flow -->
    </substrate1>

    <substrate2>
      <name_of_substrates>C</name_of_substrates>
      <initial_concentration>0.0</initial_concentration>
      <substrate_diffusion_coefficients>
        <in_pore>1e-09</in_pore>
        <in_biofilm>1e-09</in_biofilm>
      </substrate_diffusion_coefficients>
      <left_boundary_type>Dirichlet</left_boundary_type>
      <left_boundary_condition>0.0</left_boundary_condition>
      <right_boundary_type>Dirichlet</right_boundary_type>
      <right_boundary_condition>0.0</right_boundary_condition>
      <!-- C is never fed and neither reservoir stirs any back in, so both faces are
           held at zero. C leaves whichever way the flow and the gradient take it -->
    </substrate2>
  </chemistry>

  <!-- NO <microbiology> block: <biotic_mode>false</biotic_mode> skips it -->
  <!-- NO <precipitation> and NO <dissolution>: the pore space is frozen on purpose -->

  <equilibrium>
    <enabled>false</enabled>          <!-- no aqueous speciation; the four tags below are the
                                           minimum the parser accepts while disabled -->
    <components>none</components>
    <stoichiometry><species0>0.0</species0></stoichiometry>
    <logK><species0>0.0</species0></logK>
  </equilibrium>

  <diagnostics>
    <enabled>true</enabled>
    <summary_csv>summary.csv</summary_csv>   <!-- one scalar row per snapshot, beside the volumes -->
    <interval>506</interval>                 <!-- same cadence as the VTI output -->
    <tolerance>1e-6</tolerance>
    <!-- NO <conserve> entries: every face is Dirichlet, so mass crosses the boundary by design
         and no species total is constant. A conservation check here would fail every interval -->
  </diagnostics>

  <IO>
    <read_NS_file>false</read_NS_file>       <!-- start the flow from scratch, not a checkpoint -->
    <ns_rerun_iT0>0</ns_rerun_iT0>
    <read_ADE_file>false</read_ADE_file>
    <ns_filename>nsLattice</ns_filename>     <!-- the flow field: nsLattice_0000000.vti -->
    <mask_filename>maskLattice</mask_filename>
    <subs_filename>subsLattice</subs_filename>
    <!-- NOTE: this names the .chk CHECKPOINTS only. The concentration .vti files are named
         after the SPECIES, so they are A_*.vti, B_*.vti, C_*.vti -->
    <bio_filename>bioLattice</bio_filename>
    <save_VTK_interval>506</save_VTK_interval>   <!-- a snapshot every 506 steps: 21 in all -->
    <save_CHK_interval>0</save_CHK_interval>     <!-- no checkpoints; they are large and unused -->
  </IO>

</parameters>
```

### The values that change from case to case

Everything above is fixed across the sweep except these, which `make_campaign.py`
writes per case:

| tag or variable | what it is | where it comes from |
|---|---|---|
| `<Peclet>` | 0.02, 0.2 or 2.0 | the sweep |
| `<ade_max_iT>` | run length in steps | `1.25 (L^2/D)/(Pe+1)` divided by the measured timestep |
| `<save_VTK_interval>` | snapshot spacing | `ade_max_iT` over 20 |
| `<diagnostics><interval>` | CSV row spacing | the same number |
| `input/geometry.dat` | the pore space | one of the ten packs |
| `PRT_KABIO` in `env.sh` | the rate constant | `Da D / (A0 L^2)`, so Da 0.1, 1.0 or 10 |
| `PRT_DT` in `env.sh` | the timestep | measured by `build --complab` |

and one thing that is fixed per CAMPAIGN rather than per case: the boundary set,
`--boundaries dirichlet` for both the sweep and the analytic duct. It is recorded
in `campaign.json`, in every `params.json`, and as an attribute on the collected
`.h5`.

---

## 12. Every input value in the rate law, annotated

`kinetics/defineAbioticKinetics.hh`, which is the entire chemistry. Three numbers
come from the environment so one binary serves the whole sweep with no rebuild;
`make_campaign.py` writes them into each case's `env.sh` and `run_one.sh` sources
it before launching the solver.

| variable | value in this sweep | what it is |
|---|---|---|
| `PRT_KABIO` | 9.766e-4 at Da 0.1, 9.766e-3 at Da 1.0, 9.766e-2 at Da 10 | the second-order rate constant k, in L/(mol s). `k = Da * D / (A0 * L^2)` |
| `PRT_DT` | 0.0086 to 0.0090 s, measured per case | the transport timestep. Only used by the cap below, and it has to match the `[ADE] dt` the run prints |
| `PRT_MAXFRAC` | 0.25 | no single step may consume more than this fraction of the scarcer reactant |

The rate law itself, with the reasoning stripped out:

```cpp
void defineAbioticRxnKinetics(
    std::vector<double> C,          // C[0] = A, C[1] = B, C[2] = C, in mol/L
    std::vector<double>& subsR,     // OUTPUT: one rate per species, mol/L/s
    plb::plint mask                 // the voxel's material code
) {
    (void) mask;                    // deliberately unused: the reaction is volumetric,
                                    //   not gated to the grain surface
    const AbioticParams::Params& K = AbioticParams::get();   // reads PRT_* once, then caches

    for (size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;

    if (C.size() < 3 || subsR.size() < 3) return;   // size guard: a short list would
                                                    //   otherwise be read past the end

    const double A = std::max(C[0], 0.0);           // floor at zero; a clamped voxel can
    const double B = std::max(C[1], 0.0);           //   arrive at a tiny negative

    if (A <= K.MIN_CONC || B <= K.MIN_CONC) return; // the plumes have not met here yet.
                                                    //   True for most of the domain for
                                                    //   most of the run, so it is worth
                                                    //   returning before the multiply

    double R = K.k_abio * A * B;                    // THE RATE LAW. mol/L/s

    if (K.dt_kinetics > 0.0) {                      // the per-step cap
        const double scarcest = (A < B) ? A : B;    //   one cap for both, because they are
        const double capR = scarcest * K.MAX_RATE_FRACTION / K.dt_kinetics;
        if (R > capR) R = capR;                     //   consumed one for one. Capping them
    }                                               //   separately would break stoichiometry
                                                    //   the moment it stopped being 1:1

    subsR[0] = -R;      // A consumed
    subsR[1] = -R;      // B consumed
    subsR[2] = +R;      // C produced

    AbioticKineticsStats::accumulate(R);            // the per-interval summary in the log

    for (size_t i = 0; i < subsR.size(); ++i)       // a NaN here would be applied silently
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}
```

**The cap is not the clamp.** The v1.3.2 positivity clamp sits at the apply site
in the solver and refuses an increment larger than what the voxel holds; it is the
last line of defence and it is counted in the `[CLAMP]` line. The cap above keeps
the step from overshooting so far that the clamp has to fire at all, because a
clamped step is no longer the reaction the case was configured for. With these
numbers one step consumes about 0.1 per cent of the feed, so neither should ever
trigger, and `check_campaign.py` verifies that before you submit.

**The inert dissolution hook.** The file also defines `defineDissolutionRate()`,
which sets `mineralR = 0` and touches nothing else. v1.3 calls it unconditionally
from `dissolutionVOP.hh`, so the symbol has to exist even in a case with no
mineral, and this sweep has none because the geometry is frozen on purpose.

**`kinetics/defineKinetics.hh` is the biotic file and every rate in it is zero.**
`<biotic_mode>false</biotic_mode>` means it is never called, and it is still
mandatory: `complab3d_processors_part1.hh` includes it unconditionally at line 47
and `complab.cpp` calls `KineticsStats::getStats()` at every output interval
whether or not an organism exists. A tree assembled without it does not compile;
one assembled without its `KineticsStats` namespace compiles and fails to link.

**The one line to look for in every run log.** On its first reaction step the
header prints

```
[KIN] abiotic A+B->C  k_abio=0.0009766 L/(mol s)  dt_kinetics=0.0115 s  maxfrac=0.25
```

If that line is missing, the binary was built against the CompLaB defaults, the
rate constant is compiled in, `PRT_KABIO` is being ignored, and every case in the
campaign ran the same chemistry while recording a different Damköhler.

---

## 13. What this scenario does not do

Worth saying, so nobody discovers it later and thinks it was hidden.

- **The pore space never changes.** No precipitation, no dissolution.
  PRT-DeepONet learns an operator from a fixed geometry to a concentration field,
  so the geometry has to be an input rather than a variable.
- **No biology.** `<biotic_mode>false</biotic_mode>`. `defineKinetics.hh` is
  present and inert, because the solver includes it unconditionally and calls
  `KineticsStats` at every output interval whether or not any organism exists.
- **Péclet stops at 2**, for the reason in section 5. A dataset that needs
  advection-dominated transport needs a different configuration, not a wider
  sweep of this one.
- **One reaction, one stoichiometry, one diffusivity.** The sweep varies Pe, Da
  and the pore space, and holds everything else fixed. That is three axes, which
  is what 270 runs can populate honestly.
- **Two geometries per porosity** is thin. It is enough to see whether the network
  transfers between pore spaces; it is not enough to quantify how well.
