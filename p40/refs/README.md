# `refs/` — the sources, fetched and not committed

**Objective.** Local copies of the sources quoted in [`../NOTES.md`](../NOTES.md), for `../src/final_check.py` to check each quotation and fact against. Only this README and `fetch.sh` are committed, since the pages and papers keep their own terms.

**Inputs.** `sh fetch.sh` downloads them, formal-conjectures at the pin in `../lean/lakefile.toml`, and prints their hashes. The documentation was checked against the versions fetched on September 23, 2026; a later fetch of the site can differ.

| files | support |
|---|---|
| `site_*`, `thread_*`, `bib_*` | the pages of problems 28, 39, 40, 158 and 1191, the threads of §1 and §6, and the entries of the problem's two sources |
| `fc_*.lean` | the excerpts of §1.2 and §3 |
| `pliego_*`, `arxiv_*`, `cilleruelo_*`, `obryant_*`, `zenodo_*` | the papers of §3, §4 and §6 and their records |
| `oeis_A005282.txt`, `b005282.txt` | the known answer in §6 |
| `license_*.txt` | the licences in `../README.md` |

Text is compared after NFC normalisation without whitespace, since PDF text layers break lines and write ő as o plus an accent. The Zenodo record spells Erdős with U+02DD. The note's reference [2] attributes arXiv:math/0407117, O'Bryant's bibliography, to Ruzsa.
