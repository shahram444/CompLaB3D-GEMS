# Putting this repository on GitHub

Every step, in order, from the folder on your disk to a citable release. Nothing
here is specific to your machine except the two places marked **your**.

Throughout, the repository is called `CompLaB3D-GEMS` and the account
`shahram444`. If you choose a different name or account, change it in three
places at the end — step 9 says where.

---

## Before you start

You need `git` (`git --version`) and either the GitHub web interface or the
GitHub CLI (`gh --version`). Both routes are given; use whichever you have.

You also need the repository folder itself — the one containing `README.md`,
`src/`, `examples/`. Everything below is run from inside it.

```bash
cd /path/to/CompLaB3D-GEMS      # your folder
ls                              # you should see README.md, src/, examples/
```

---

## 1. Check the tree before you publish anything

```bash
./tests/check_repo.sh
```

Eight checks, one second, no dependencies. It should end with
`All structural checks passed.` If it does not, fix what it names first — it is
much easier than fixing it after the repository is public.

---

## 2. Tell git who you are

Only needed once per machine. Skip if `git config user.name` already answers.

```bash
git config --global user.name  "Shahram Asgari"
git config --global user.email "shahram.asgari@uga.edu"
```

The email you use here appears in every commit. Use the one attached to your
GitHub account, or GitHub will not link the commits to you.

---

## 3. Start the repository

```bash
git init
git branch -M main
```

`git init` creates the hidden `.git` folder. `git branch -M main` names the
branch `main`, which is what GitHub expects.

---

## 4. Check what git is about to include

`.gitignore` already excludes build output, `run/`, VTI files and Python caches.
Confirm that before you commit, because a stray 2 GB of results is very hard to
remove afterwards.

```bash
git add -A
git status --short | head -50          # the first fifty files
git status --short | wc -l             # how many in total
```

You should see roughly 180 files and nothing under `run/`, `build/` or `*.vti`.
If something unwanted appears, add a line for it to `.gitignore`, then

```bash
git reset
git add -A
```

and look again.

---

## 5. The first commit

```bash
git commit -m "CompLaB3D-GEMS v1.0

Three-dimensional pore-scale reactive transport with an evolving pore space
and six interchangeable rate paths: compiled kinetics, flux balance analysis
through GLPK and through COBRApy, a fitted surrogate network, a symbolic rate
law, and a graph network.

Sixteen examples, the offline pipelines they need, and the regression suite."
```

---

## 6. Create the repository on GitHub

**With the web interface.** Go to <https://github.com/new>. Set the name to
`CompLaB3D-GEMS`, the description to *Three-dimensional pore-scale reactive
transport with an evolving pore space and six rate paths*, and choose **Public**.
Do **not** tick "Add a README", "Add .gitignore" or "Choose a license" — you have
all three already, and ticking them creates a conflicting first commit.

Then connect your folder to it:

```bash
git remote add origin https://github.com/shahram444/CompLaB3D-GEMS.git
git push -u origin main
```

**With the GitHub CLI.** One command does both:

```bash
gh auth login                   # once per machine
gh repo create CompLaB3D-GEMS --public --source=. --remote=origin --push \
   --description "Three-dimensional pore-scale reactive transport with an evolving pore space and six rate paths"
```

### If the push asks for a password

GitHub stopped accepting account passwords over HTTPS. Either use `gh auth login`
(which handles it), or create a personal access token at
<https://github.com/settings/tokens> with the `repo` scope and paste that as the
password. To avoid retyping it:

```bash
git config --global credential.helper store     # Linux
git config --global credential.helper osxkeychain   # macOS
```

---

## 7. Set the repository up so people can find it

On the repository page, click the gear beside **About** and add topics:
`lattice-boltzmann`, `reactive-transport`, `pore-scale`, `palabos`,
`flux-balance-analysis`, `biogeochemistry`, `computational-fluid-dynamics`,
`surrogate-model`, `symbolic-regression`, `graph-neural-network`.

Still on that panel, tick **Releases** and **Packages** off if you do not want
them in the sidebar, and paste the description if you did not set it in step 6.

Then, under **Settings → General → Features**, turn **Issues** on and
**Wikis** off — the manual is in `docs/`, and a wiki that repeats it will drift.

---

## 8. Make a release, and get a DOI

A release is what people cite. It also freezes a version, so that someone who
reproduces your results a year from now gets the same code.

**Connect Zenodo first, before tagging.** Zenodo only sees releases made after
you switch it on:

1. Sign in at <https://zenodo.org> with your GitHub account.
2. Go to <https://zenodo.org/account/settings/github/> and flip the switch
   beside `shahram444/CompLaB3D-GEMS`.

Then tag and release:

```bash
git tag -a v1.0 -m "CompLaB3D-GEMS v1.0"
git push origin v1.0
```

On GitHub, go to **Releases → Draft a new release**, choose the `v1.0` tag,
title it `CompLaB3D-GEMS v1.0`, and paste the `## v1.0` section of
[`../CHANGELOG.md`](../CHANGELOG.md) as the body. Publish.

Within a few minutes Zenodo will mint a DOI. Copy it, and put it in two places:

- `CITATION.cff` — the `doi:` field;
- `README.md` — a badge under the title:

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

Commit and push that change. It does not need a new release.

---

## 9. If you chose a different name or account

Three files hard-code `shahram444/CompLaB3D-GEMS`. Change all three together:

| File | What to change |
|---|---|
| `README.md` | the `git clone` line under **Start here** |
| `CITATION.cff` | `repository-code:` and `url:` |
| `docs/PUBLISHING.md` | this file, so it stays true |

```bash
grep -rn "shahram444/CompLaB3D-GEMS" .    # find every occurrence
```

---

## From now on: the edit loop

This is the whole of day-to-day use. Four commands.

```bash
git pull                        # start from what is on GitHub
# ... edit ...
./tests/check_repo.sh           # cheap, catches structural mistakes
git add -A
git commit -m "what changed, and why"
git push
```

`git pull` first matters as soon as you work from two machines, or once anyone
else has push access.

### Working on something risky

Branch, so `main` always builds:

```bash
git checkout -b transient-solver
# ... edit, commit as often as you like ...
git push -u origin transient-solver
```

Then open a pull request on GitHub, read your own diff, and merge. Reading the
diff before merging catches more than anything else in this document.

### If you commit something you should not have

Not yet pushed:

```bash
git reset --soft HEAD~1         # undo the commit, keep the edits
```

Already pushed, and it was large or private: rewriting published history is
awkward, and a file that has been public should be treated as public. Delete it,
commit the deletion, and rotate anything secret. `git filter-repo` can scrub the
history if the file was genuinely never fetched, but do not rely on that.

---

## What not to commit

`.gitignore` covers the usual cases. The ones worth knowing by name:

| Do not commit | Why |
|---|---|
| `run/` | Assembled working directories. `setup_case.sh` recreates them exactly. |
| `build/`, `*.o`, `complab` | Build output, machine-specific. |
| `*.vti`, `*.dat` results | Output fields. These grow to gigabytes and belong in a data archive with its own DOI, not in git. |
| Large genome-scale models | The three in `models/` are gzipped and small. A full BiGG model is not; link to it instead. |
| Cluster paths, usernames, allocation numbers | They stop working the moment anyone else clones this. |

For output too large for git and too important to lose, put it in Zenodo or in
your institution's archive, and reference the DOI from the README. GitHub
refuses single files over 100 MB, and warns above 50 MB.
