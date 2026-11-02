# `refs/` — the sources, fetched and not committed

**Objective.** Local copies of the sources quoted in [`../NOTES.md`](../NOTES.md), against which `../src/final_check.py` checks each quotation and fact. Only this README and `fetch.sh` are committed.

**Inputs.** `sh fetch.sh` downloads them and prints their hashes. Since October 6, 2026 the site shows no statuses, so they come from `teorth_problems.yaml`.

| files | support |
|---|---|
| `site_*`, `thread_*`, `bib_*`, `blog_9.html` | the problem pages, the threads of §1 and §6, the problem's two sources, the status change |
| `teorth_problems.yaml` | the statuses of §1.1 and §3 |
| `fc_*.lean` | the excerpts of §1.2 and §3, at the pin |
| `fchead_*` | formal-conjectures' main branch and its commit lists, §1.2 and §3 |
| `pliego_*`, `arxiv_*`, `cilleruelo_*`, `obryant_*`, `zenodo_*` | the papers of §3, §4 and §6; `obryant_*.tex` has the reviews the PDF omits |
| `oeis_A005282.txt`, `b005282.txt` | the known answer in §6 |
| `license_*.txt` | the licences in `../README.md` |

Text is compared NFC-normalised, without whitespace: PDF text layers break lines and write ő as o plus an accent. The Zenodo record spells Erdős with U+02DD. The note's reference [2] attributes arXiv:math/0407117, O'Bryant's bibliography, to Ruzsa.
