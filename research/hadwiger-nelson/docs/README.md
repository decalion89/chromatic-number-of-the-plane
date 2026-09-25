# Docs

[`note/`](note/) holds the three-page note on the two theorems,
[`planes-4-chromatic.pdf`](note/planes-4-chromatic.pdf), and its HTML source. The PDF is rendered
with headless Chromium:
`chromium --headless --no-pdf-header-footer --print-to-pdf=planes-4-chromatic.pdf planes-4-chromatic.html`.

[`research-log.md`](research-log.md) is the project's complete lab notebook, about 8 000 lines. It
records every experiment, result, dead end, correction and withdrawn claim, in the order they
happened.

[`figures/`](figures/) holds the SVG figures used in the READMEs. `scripts/make_figures.py` regenerates
them from the data.

Read the log to see how a result was reached or why an approach was abandoned. For what is established,
start with the project [README](../README.md) and the [notes](../notes/).
