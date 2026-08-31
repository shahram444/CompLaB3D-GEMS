# Bundled genome-scale models — provenance and licence

These three models are **not part of CompLB3D**. They are third-party scientific
data, bundled so that the examples run with no network access, which matters
because cluster compute nodes generally cannot reach the internet.

| file | model | metabolites × reactions | source |
|---|---|---|---|
| `e_coli_core.xml.gz` | *E. coli* core metabolism | 72 × 95 | BiGG |
| `iJO1366.xml.gz` | *E. coli* K-12 MG1655 | 1805 × 2583 | BiGG |
| `STM_v1_0.xml.gz` | *Salmonella enterica* Typhimurium LT2 | 2436 × 3357 | BiGG |

Obtained from the copies distributed with [cobrapy](https://github.com/opencobra/cobrapy),
which are themselves from [BiGG Models](http://bigg.ucsd.edu/).

---

## Read this before publishing or redistributing

> **The licence on these files has not been verified.**
>
> Someone needs to read <http://bigg.ucsd.edu/license> and confirm that
> redistribution inside an AGPL-3.0 repository is permitted, and under what
> attribution conditions, **before** this repository is cited in a paper or
> archived on Zenodo.
>
> This is the same category of question as vendoring Palabos, and it was
> answered there by *not* bundling it.

If the answer turns out to be no, nothing breaks:

```bash
rm models/*.xml.gz
```

The manifest still lists every model with its URL, so `<allow_download>` fetches
them on demand, and `tools/` and `extractMM.py` are unaffected. **The bundling is
a convenience, not a dependency.**

## Citing them

Each model comes from its own publication and carries its own citation
obligation, separate from BiGG's and from CompLB3D's. If you run a simulation on
one of these, cite the model's source paper. BiGG's page for each model names it.

- BiGG Models: King ZA *et al.* (2016) *Nucleic Acids Research* 44:D515–D522.
- iJO1366: Orth JD *et al.* (2011) *Molecular Systems Biology* 7:535.
- *E. coli* core: Orth JD, Fleming RMT, Palsson BØ (2010) *EcoSal Plus*.
- STM_v1_0: Thiele I *et al.* (2011) *BMC Systems Biology* 5:8.

## Why they are gzipped

Unpacked they are 21 MB; gzipped they are under 1 MB. The solver decompresses
into the cache directory on first use and reads the plain XML afterwards, so the
cost is paid once.

## Adding your own

Append a row to `manifest.txt`:

```
name  nmet  nrxn  objective  hash  bytes  url
```

`hash` may be `0` to skip the checksum. `nmet` and `nrxn` are the important
columns: they are what lets the solver notice that a model has been revised
under you, which silently shifts every positional `<exchange_reaction_indices>`
entry. Set them to `0` to skip that check too, but read the warning in
`complab3d_fetch.hh` first.
