# The note

[`planes-4-chromatic.pdf`](planes-4-chromatic.pdf) is a three-page note proving χ(ℚ(√2, √3)²) = 4
and giving a short proof of K. G. Fischer's theorem χ(ℚ(√3, √11)²) = 4 (1994).
[`planes-4-chromatic.html`](planes-4-chromatic.html) is its source. The PDF is rendered with headless
Chromium, run in this folder:

```sh
chromium --headless --no-pdf-header-footer --print-to-pdf=planes-4-chromatic.pdf planes-4-chromatic.html
```

Chromium does not write an author into the PDF; version 4 has it added afterwards with `pypdf`
(`PdfWriter.add_metadata`).

The fuller account, with the checks behind each step, is
[`notes/local_colourings.md`](../../notes/local_colourings.md) §8 and §10; the tests are
`tests/test_q311.py` and `tests/test_q23.py`. The note has not been refereed.
