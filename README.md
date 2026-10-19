# erdos-problems

Work against the [Erdős Problems](https://www.erdosproblems.com) database. This public repository contains the five problem folders; the tracker and Lean corpora stay local.

```
p40/                 the density that forces unbounded representation counts
p200/                longest AP of primes in {1,...,N}
p1086/               triangles of one area
p1162/               number of subgroups of S_n
p1186/               monochromatic k-APs in 2-colourings
tracker-and-solver/  local working folder, git-ignored and not published here
formal/              local third-party Lean corpora, git-ignored and not published here
```

`p<N>/` is sustained work on a single problem N. `p200/` and `p1162/` are the worked examples; `p1086/`, `p1186/` and `p40/` also have written READMEs.

The five problems were separate public repositories (`incrediblecrab/erdos-40`, `incrediblecrab/erdos-200`, `incrediblecrab/erdos-1086`, `incrediblecrab/erdos-1162` and `incrediblecrab/erdos-1186`) until October 19, 2026, when they were merged here with their histories.

The local `tracker-and-solver/` holds the database: all 1221 problems as erdosproblems.com publishes them, a plain-English summary of each, a private layer for recording what we think about them, and for the 586 solved ones the answer, the papers that settled it, and a reconstructed route.

If you have the local tracker, start at `tracker-and-solver/README.md`:

```bash
cd tracker-and-solver
./erdos.py stats
```

`tracker-and-solver/data/erdos.db` is the shared source of truth locally for problem text, status, tags and bibliography; it is not in this repository. Local work should read it, not re-scrape the site.

Everything here runs on the shared interpreter at `~/.venvs/main/bin/python`.

## Attribution and third-party material

The problem statements, status, tags and bibliography in the local `tracker-and-solver/data/erdos.db` come from [erdosproblems.com](https://www.erdosproblems.com), curated by Thomas Bloom. The community database mirroring that site, [teorth/erdosproblems](https://github.com/teorth/erdosproblems), is released under the Apache License 2.0; the local tracker reuses that content under those terms, but the database is not redistributed here. Problem numbering follows erdosproblems.com.

No third-party papers are redistributed here. Each problem's `refs/` directory keeps fetched papers git-ignored and carries a `README.md` naming the sources, plus a `fetch.sh` where the sources can be re-downloaded. Run those scripts to populate them locally. `tracker-and-solver/data/papers/` belongs to the local tracker and is not part of this repository.

Nothing in this repository claims to have solved an open Erdős problem. `p40/README.md`, `p200/README.md`, `p1086/README.md`, `p1162/README.md` and `p1186/README.md` state the status of their problem and what the directory does and does not establish.
