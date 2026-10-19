# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) and GitHub Copilot CLI when working in this repository.

The public repository contains `p40/`, `p200/`, `p1086/`, `p1162/` and `p1186/`, each for one problem, plus root documentation.

Read each problem folder's own `CLAUDE.md` before working there; it holds that folder's commands and rules. Work on one problem must not change another.

`tracker-and-solver/` and `formal/` are local, git-ignored working folders and must never be committed: this repository is public, the tracker holds private notes, and `formal/` clones third-party Lean corpora.

`refs/` directories hold third-party papers fetched locally, never committed; only their source lists and fetch scripts are tracked.

The shared interpreter is `~/.venvs/main/bin/python`.
