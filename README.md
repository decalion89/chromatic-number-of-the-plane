# Research in combinatorial geometry: the Hadwiger–Nelson problem

[![tests](https://github.com/decalion89/darwin-50/actions/workflows/tests.yml/badge.svg)](https://github.com/decalion89/darwin-50/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)

This repository holds an ongoing computational and arithmetic attack on the
**Hadwiger–Nelson problem**. The problem asks for the chromatic number χ(ℝ²) of
the plane: the least number of colours such that no two points at distance
exactly 1 share a colour. Since 2018 it has been known that χ(ℝ²) ∈ {5, 6, 7}.
The goal here is the lower bound **χ(ℝ²) ≥ 6**.

**→ [`research/hadwiger-nelson/`](research/hadwiger-nelson/README.md)**: results,
code, data, certificates and the full research log.

## Highlights

- **χ(ℚ(√3, √11)²) = χ(ℚ(√2, √3)²) = 4.** Both planes are 4-chromatic.
  - The first question was left open by Moorhouse (2010), Madore (2015) and
    Exoo–Ismailescu (2018).
  - Voronov (2021) conjectured both values.
  - The proof extends a 2-adic colouring used by Speyer in Polymath16 (2018).
  - It has not been refereed.
- **Machine-checked reproductions** of χ(ℝ²) ≥ 4 and χ(ℝ²) ≥ 5. De Grey's
  1581-vertex graph was rebuilt from its seed and certified with a DRAT proof
  verified by `drat-trim`. The proof fixes the colours of one triangle, a
  standard symmetry step.
- **Two-distance witnesses** for the six-colour route, each verified by four SAT
  solvers and `drat-trim`. One of them needs a single unit-distance gadget,
  thanks to a Galois symmetry of its field.
- **χ(ℝ²) ≥ 6 has not been proved.** The searches for the missing gadget are
  running.

## Repository layout

```
research/
  hadwiger-nelson/
    README.md          overview of results, with evidence and references
    hn/                Python library: exact number fields, unit-distance graphs,
                       SAT colouring, certificates, local (adelic) colourings
    tests/             pytest suite
    scripts/           maintained tools, indexed in scripts/README.md;
                       scripts/experiments/ keeps the one-off runs
    data/              graphs and witnesses in exact coordinates (data/README.md)
    certificates/      colourings and DRAT verification logs (certificates/README.md)
    notes/             technical notes and literature comparison
    docs/              the chronological research log
```

## Method and provenance

All coordinates are exact algebraic numbers. Every non-colourability claim is
rebuilt from the coordinates and decided by independent SAT solvers. It is then
certified by a DRAT proof that `drat-trim` checks.

This is AI-assisted research, carried out with Claude (Anthropic) through
Claude Code under the direction of Sergi Galán. No result here has been
peer reviewed.

## Citing

GitHub's "Cite this repository" button uses [`CITATION.cff`](CITATION.cff). Please cite the
original papers listed in the research README as well.

## Reporting an error

Corrections are welcome. If a claim here fails to reproduce, open an issue naming:
- the file;
- the command you ran;
- what you saw.

Every graph can be rebuilt with `research/hadwiger-nelson/scripts/verify_pair.py`.

## License

Code, data and text are released under the [MIT License](LICENSE).

## Resumen en español

Investigación sobre el problema de Hadwiger–Nelson: el número cromático del
plano, que se sabe que es 5, 6 o 7. El objetivo es demostrar que es al menos 6.

Resultado principal hasta ahora: los planos con coordenadas en ℚ(√3, √11) y en
ℚ(√2, √3) son 4-cromáticos. Lo primero era una pregunta abierta desde 2010, y
Voronov conjeturó ambos casos en 2021. La prueba extiende una idea 2-ádica de
Speyer (Polymath16, 2018).

Todo cálculo se verifica con varios SAT solvers y con pruebas DRAT. El trabajo
se ha hecho con ayuda de IA y aún no ha sido revisado por pares.
