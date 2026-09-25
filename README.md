# The Hadwiger–Nelson problem over number fields

**Sergi Galán** · research repository, 2026

[![tests](https://github.com/decalion89/darwin-50/actions/workflows/tests.yml/badge.svg)](https://github.com/decalion89/darwin-50/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)

The Hadwiger–Nelson problem asks for the chromatic number χ(ℝ²) of the plane:
the least number of colours such that no two points at distance exactly 1 share
a colour. It is known that 5 ≤ χ(ℝ²) ≤ 7. This repository studies the problem
through exact arithmetic:
- every point has algebraic coordinates, so "distance 1" is decided exactly;
- every claim that a graph cannot be coloured comes with a proof certificate
  checked by an independent program.

## Main result

> **Theorem.** χ(ℚ(√3, √11)²) = χ(ℚ(√2, √3)²) = 4.

The plane with coordinates in ℚ(√3, √11) can be coloured with four colours.
That field is the smallest one whose plane contains a Moser spindle. The same
holds for the plane over ℚ(√2, √3).

- **Background.** Madore (2015) proved 4 ≤ χ(ℚ(√3, √11)²) ≤ 5. Moorhouse
  (2010) and Exoo–Ismailescu (2018) left the value open. Voronov (Polymath16,
  2021) conjectured both values.
- **Proof idea.** Change coordinates to α = x + y/√3, β = 2y/√3. The squared
  distance becomes α² − αβ + β². This form is anisotropic modulo a prime above
  2 with residue field 𝔽₂, so Madore's reduction argument applies. Reducing
  (α, β) modulo that prime gives a proper 4-colouring of the whole plane. Speyer
  used reduction modulo 2 in these coordinates in 2018 to 4-colour the Moser
  ring.
- **Status.** Proved, and checked by computer. Not yet refereed.

**Read:** [the three-page note (PDF)](research/hadwiger-nelson/docs/note/planes-4-chromatic.pdf) ·
[full details](research/hadwiger-nelson/notes/local_colourings.md) (§8 and §10) ·
[comparison with the literature](research/hadwiger-nelson/notes/literature.md)

## Other results

| result | evidence |
|---|---|
| χ(ℝ²) ≥ 4: the Moser spindle has no 3-colouring | DRAT proof, checked by `drat-trim` |
| χ(ℝ²) ≥ 5: de Grey's 1581-vertex graph, rebuilt from its 39-point seed, has no 4-colouring | DRAT proof of 13.1 M lemmas, checked by `drat-trim` |
| χ(ℚ(√−3, √−11)) = 4 and χ(ℚ(√−3, √−11, √−247)) = 5, for the whole complex fields | proofs in the notes, with unit tests |
| Two-distance graphs with no 5-colouring: 187 points and 72 points. The first needs a single unit-distance gadget to prove χ(ℝ²) ≥ 6 | four SAT solvers, and DRAT proofs checked by `drat-trim` |

**χ(ℝ²) ≥ 6 has not been proved.** The search for the missing gadget is
described in the [project page](research/hadwiger-nelson/README.md).

## Repository layout

| path | contents |
|---|---|
| [`research/hadwiger-nelson/`](research/hadwiger-nelson/README.md) | project page: results, evidence and references |
| [`research/hadwiger-nelson/hn/`](research/hadwiger-nelson/hn/) | Python library: exact number fields, unit-distance graphs, SAT colouring, certificates |
| [`research/hadwiger-nelson/notes/`](research/hadwiger-nelson/notes/) | technical notes and the literature comparison |
| [`research/hadwiger-nelson/docs/`](research/hadwiger-nelson/docs/) | the note, figures and the chronological research log |
| [`research/hadwiger-nelson/data/`](research/hadwiger-nelson/data/) | graphs and witnesses in exact coordinates |
| [`research/hadwiger-nelson/certificates/`](research/hadwiger-nelson/certificates/) | colourings and verified proofs |
| [`research/hadwiger-nelson/tests/`](research/hadwiger-nelson/tests/) | test suite |
| [`research/hadwiger-nelson/scripts/`](research/hadwiger-nelson/scripts/) | verification, search and figure tools |

Every folder has a README describing its contents.

## Reproducing

```sh
cd research/hadwiger-nelson
python3 -m pip install -r requirements.txt
python3 -m pytest -q tests/test_q311.py tests/test_q23.py   # the two theorems, in seconds
```

`scripts/verify_pair.py` rebuilds any graph from its data file, recomputes every
edge exactly and runs the solvers. GitHub Actions runs the fast part of the test
suite on every push and pull request.

## Citing

Use GitHub's "Cite this repository" button, which reads
[`CITATION.cff`](CITATION.cff). Please also cite the original papers listed in
the project page.

## How this work was done

This is AI-assisted research. The code, the computations and most of the text
were produced with Claude (Anthropic), under the direction of Sergi Galán. Every
computational claim is machine-checked as described above. No result has been
peer reviewed yet.

## Reporting an error

Corrections are welcome. If a claim fails to reproduce, please open an issue
giving the file, the command you ran, and what you saw.

## License

Code, data and text are released under the [MIT License](LICENSE).

## Resumen en español

Este repositorio estudia el problema de Hadwiger–Nelson: el número cromático del
plano, que se sabe que está entre 5 y 7. Se trabaja con aritmética exacta en
cuerpos de números. Toda afirmación de que un grafo no se puede colorear va
acompañada de un certificado verificado por un programa independiente.

Resultado principal: los planos con coordenadas en ℚ(√3, √11) y en ℚ(√2, √3)
tienen número cromático exactamente 4. Lo primero era una pregunta abierta
desde 2010, y Voronov conjeturó ambos casos en 2021. La prueba cambia de
coordenadas para que el argumento de reducción de Madore funcione módulo 2, y
extiende a todo el plano un coloreado de Speyer (2018). Está explicada en una
[nota de tres páginas](research/hadwiger-nelson/docs/note/planes-4-chromatic.pdf).
Es un trabajo hecho con ayuda de IA y todavía no ha sido revisado por pares.
