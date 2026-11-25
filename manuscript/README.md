# Singular squared-distance Laplacians

This note resolves both clauses of Das–Mohanty's concluding conjecture in *Electronic Journal of Linear Algebra* 40 (2024), p. 488, for the Laplacian-like matrix in equation (3.7). Its arXiv counterpart is on p. 13 of version 2, with matrix equation (3.6).

For positive integer part sizes and a singular squared-distance matrix, the positive eigenvalues are exactly one quarter with multiplicity N−t and one with multiplicity h−1; the inertia, ordered positive/zero/negative, is (N−s−1, 1, s). The note also proves the general nonzero-cofactor inertia formula and corrects the source's Theorem 4.4(ii), with the exact counterexample K_{1,1,3,3}.

## Additional mathematical coverage

The conjectured positive spectrum holds throughout the full negative-ratio domain R<0, rather than just on the singular locus R=−1; zero multiplicities are omitted. This follows by matching the explicit eigenspace dimensions to the general positive index.

The characteristic polynomial determines every remaining negative eigenvalue in the singular case: repeated nonsingleton sizes supply repeated diagonal eigenvalues, while the other negative roots interlace strictly between their distinct diagonal values, with one root below them all. Further corollaries give the singular distance nullspace and adjugate, inverse formulas including the zero-cofactor boundary, singleton bounds and their equality cases, and the previously established squared-distance inertia and energy formulas.

The same congruence argument answers both parts of [Lucid's MathOverflow question 513954, revision 2](https://mathoverflow.net/revisions/513954/2): the structured-downdate inertia formula holds even for singular inputs, and a continuous positive-semidefinite family on a finitely punctured open simplex has a constant pivot sign and, in the negative branch, a fixed negative diagonal index. The note proves the inertia identity for arbitrary symmetric inputs and the family assertion on any nonempty connected parameter space. Lucid's question already contains the decisive hyperplane congruences; this attribution concerns the mathematical mechanism, not a priority claim.

## Files and reproduction

- `note.tex`: self-contained LaTeX source with its bibliography.
- `note.pdf`: readable mathematical note.
- `metadata.json`: title, creators, coverage and source identifiers.
- `verify.py`: exact rational matrix checks, requiring Python 3.9 or later and no third-party modules.
- `verification.json`: deterministic output of the finite checks, supplied alongside the source archive.

The source archive contains the TeX, this README, metadata, the checker, its deterministic output and a hash manifest. The PDF and deterministic checker output also accompany it separately.

Build the PDF using Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, with a standard TeX Live installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded; BibTeX and external figures are unnecessary. Required LaTeX packages are geometry, amsmath, amssymb, amsthm, booktabs, array, longtable, fontenc, lmodern, hyperref and xurl.

Reproduce the finite checks:

```sh
python3 verify.py > verification.json
sha256sum verification.json
```

Expected SHA256 for `verification.json`:

```text
37a34d9db113c95d1ee8c7d97719996f9905ea60664279db076cb7c8556f4ec7
```

The checker enumerates all 2,834 nondecreasing lists of two through eight part sizes in the range one through eight, with at most 24 vertices. The 16 lists with vanishing cofactor sum are outside the definition of L. For the other 2,818 lists it verifies the original omitted-list entries, general inertia and exact multiplicities. It checks full matrices for all 11 singular examples and the characteristic polynomial of K_{1,1,3,3}. These are finite checks; the complete arguments for arbitrary positive integer sizes are in the note.

Additional exact checks cover all 2,834 distance quotients, all 2,823 nonsingular distance inverses including the 16 zero-cofactor examples, all 2,818 defined Laplacian characteristic polynomials, and full characteristic polynomials in the 11 singular cases. Exact Sturm sequences count the residual negative roots in their asserted intervals. The diagonal structured-downdate check enumerates 3,615 pairs in dimensions one to three, using diagonal entries −2 through 2 and vector entries −1 through 1; it checks all 2,816 nonzero-pivot pairs, including 1,272 singular inputs, for both inertia and the direct-sum kernel formula. Connected-family assertions have written proofs rather than finite sampling arguments.

## Mathematical source coverage

[The pinned Lean file](https://github.com/the-omega-institute/trureturing/blob/a2d277ea0980f1230eff2c7e27838372effccb30/D5/S3/Combinatorics/SpectralGraph/DasMohantySingularSquaredDistance.lean) contains `result` in namespace `D5.S3.Combinatorics.SpectralGraph.DasMohantySingularSquaredDistance`, proving both clauses of the singular conjecture on the original nonzero-cofactor domain. `singular_source_downdate` proves its rank-one identity, and `singular_source_reduction` gives inertia transfer and the two positive-eigenvector families. The general nonzero-cofactor formula, correction of Theorem 4.4(ii), and extensions of cofactors and exact multiplicities have written proofs in the note. The script supplies the finite evidence described above.

[StructuredRankOneInertia.lean](https://github.com/the-omega-institute/trureturing/blob/a2d277ea0980f1230eff2c7e27838372effccb30/D5/S3/QuadraticForms/StructuredRankOneInertia.lean), in namespace `D5.S3.QuadraticForms.StructuredRankOneInertia`, contains `inertia_downdate`, `kernel_split` and `nullity_downdate` for arbitrary real symmetric inputs with nonzero pivot.

[StructuredRankOneFamily.lean](https://github.com/the-omega-institute/trureturing/blob/a2d277ea0980f1230eff2c7e27838372effccb30/D5/S3/QuadraticForms/StructuredRankOneFamily.lean), in namespace `D5.S3.QuadraticForms.StructuredRankOneFamily`, contains `mo513954_revision2`, covering the diagonal-input inertia formula and the original finite-puncture simplex family statement; `punctured_simplex_pathConnected` and `continuous_family` supply its two topological steps. The extension to all connected parameter spaces is a written proof in the note. These two linked source files are distinct from the graph-specific declaration.

## Creators and license

1. Haobo Ma — ChronoAI Pte Ltd, Singapore.
2. Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

CC BY 4.0 (`cc-by-4.0`).
