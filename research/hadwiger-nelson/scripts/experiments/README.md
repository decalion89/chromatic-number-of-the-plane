# Exploratory scripts

These 730 scripts are the one-off experiments behind the research log
([`../../docs/research-log.md`](../../docs/research-log.md)). They are kept as a
record of what was tried, not as maintained tools. The maintained tools are in
[`../`](../).

Things to know before running one:

- **Paths.** Each script finds the project root from its own location, so the
  `hn` package, `data/` and the tools in `scripts/` resolve in any checkout.
- **Missing inputs.** Many scripts read intermediate files from a temporary
  working directory (paths under `/tmp/...`). Those files were never committed.
  Such scripts document a computation but will not run as they stand.
- **Superseded results.** A script may encode a hypothesis that the log later
  corrected or withdrew. The log, not the script, says which results stand.
