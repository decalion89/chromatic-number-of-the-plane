# Exploratory scripts

These 733 scripts are the one-off experiments behind the research log
([`../../docs/research-log.md`](../../docs/research-log.md)). They are kept as a
record of what was tried, not as maintained tools. The maintained tools are in
[`../`](../).

Things to know before running one:

- **Paths.** Each script finds the project root from its own location, so the
  `hn` package, `data/` and the tools in `scripts/` resolve in any checkout.
- **Missing inputs.** Many scripts read intermediate files from a temporary
  working directory, written `/tmp/hn`. Those files were never committed, so
  such scripts document a computation but will not run as they stand. (The
  scripts first named the temporary directory of the session that ran them;
  only that path string has since been changed.)
- **Superseded results.** A script may encode a hypothesis that the log later
  corrected or withdrew. The log, not the script, says which results stand.
