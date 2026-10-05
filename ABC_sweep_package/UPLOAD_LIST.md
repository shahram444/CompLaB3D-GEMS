# What to upload to the cluster, and where

Everything goes under `/scratch/sa01687/ABC_sweep/`, into the folder of the
same name it has on the PC. Overwrite what is there.

## Upload these 8. The run needs them.

| On the PC, under `ABC_sweep_package\` | Goes to | What changed |
|---|---|---|
| `tools\make_geometries.py` | `/scratch/sa01687/ABC_sweep/tools/` | 10 porosity levels x 3 packings = 30 geometries, and the new seed blocking |
| `tools\make_campaign.py` | `/scratch/sa01687/ABC_sweep/tools/` | Pe 0.02 0.2 2.0, Da 0.1 1.0 10, dirichlet default, throttle 12 |
| `tools\check_campaign.py` | `/scratch/sa01687/ABC_sweep/tools/` | the gate, with the boundary-set wording corrected |
| `tools\collect_to_h5.py` | `/scratch/sa01687/ABC_sweep/tools/` | streams each run to disk instead of holding all 270 in RAM |
| `tools\verify_dataset.py` | `/scratch/sa01687/ABC_sweep/tools/` | consecutive-pair sensitivity, held-zero face checks, spread sampling |
| `tools\selftest.py` | `/scratch/sa01687/ABC_sweep/tools/` | builds a 24-case subset, and pins all four x planes |
| `tools\analytic_channel.py` | `/scratch/sa01687/ABC_sweep/tools/` | the duct rebuild hint now names campaign02 values |
| `xml\CompLaB.xml.template` | `/scratch/sa01687/ABC_sweep/xml/` | boundary comments now describe the all-Dirichlet default |

## Upload these 2 as well, to keep the copies identical

| On the PC | Goes to | Note |
|---|---|---|
| `kinetics\defineAbioticKinetics.hh` | `/scratch/sa01687/ABC_sweep/kinetics/` | comment only. Upload to keep the two copies identical; no rebuild needed |
| `patch\apply_damkohler_feed_reference.py` | `/scratch/sa01687/ABC_sweep/patch/` | comment only. Same |

**No rebuild of the solver.** The only change to a compiled file is a comment,
so the binary at `~/abc_build/complab` stays as it is.

## Do not upload these. They are already correct on the cluster.

- `kinetics/defineKinetics.hh`, unchanged
- `tools/train_prt3d.py`, unchanged
- `tools/setup_build.sh`, unchanged

## Do not upload these. They stay on the PC.

- `tools/make_design_docx.py`, needs Node and the docx package. It rebuilds the design document; it is not part of a run
- `README_1_SCENARIO.md`, documentation
- `README_2_STEP_BY_STEP.md`, documentation
- `campaign02_design.docx`, documentation

---

## Check the upload before running anything

Manual upload is where a file arrives truncated, or in the wrong folder, and
nothing downstream tells you: the campaign builds, the array runs, and the
mistake shows up as a strange dataset a day later. These are the checksums of
the files as they left the PC. On the cluster:

```bash
cd /scratch/sa01687/ABC_sweep
md5sum -c UPLOAD_MANIFEST.txt
```

Upload `UPLOAD_MANIFEST.txt` too, into `/scratch/sa01687/ABC_sweep/`. Every
line must say `OK`. A `FAILED` line names the file to upload again.

```
bcb37f0d592cbf1ddd71b416eeabd0bc  tools/make_geometries.py
6f43754bb1e48dc433e82a0a8f3059cb  tools/make_campaign.py
a1406b1d77c2e3402e68dfcff0355292  tools/check_campaign.py
04bf7177b4d32b15cd18365d08de241b  tools/collect_to_h5.py
d0fa59d40a41187b3828021e77b37dc2  tools/verify_dataset.py
cbda7ff491e52fd3a8afb5326072b886  tools/selftest.py
9991065b85620c0c1578a4a0ef8cfcbe  tools/analytic_channel.py
5ba7e2fbe5b58e89f40cac1ce015e83a  xml/CompLaB.xml.template
439756cee4f34c185108086dbf0949b5  kinetics/defineAbioticKinetics.hh
e070414072be7a65b024643f5ce89552  patch/apply_damkohler_feed_reference.py
0b796b455d63ff3381a34cbdfe95a657  kinetics/defineKinetics.hh
4791f8433baca5a840f95408cd15664b  tools/train_prt3d.py
02719906fe1c1118b43ee47a2b96f519  tools/setup_build.sh
```
