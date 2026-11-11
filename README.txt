SCHMEISSER'S CONJECTURE FOR POLYNOMIALS OF DEGREE AT LEAST 6251

Qiyue Tang
Zhili College, Tsinghua University, Beijing 100084, China
tangqy24@mails.tsinghua.edu.cn
November 11, 2026

MAIN RESULT

If p is a complex polynomial of degree n >= 6251 and all its zeros lie
in the closed unit disk, every point in the convex hull of its zeros is
at distance at most one from a zero of p'. Zeros and critical points are
counted with multiplicity wherever sums or products are used.

The radius one is sharp, as shown by p(z) = z^n - 1.
All degree-dependent arguments cover the entire ray n >= 6251.

CONTENTS

Schmeisser_6251.pdf
    The complete article, including proofs, exact parameter tables,
    the finite-verification predicates, and references.

Schmeisser_6251.tex
    A self-contained LaTeX source file, including every section,
    bibliography entry, formula, and printed table.

manuscript_modular/
    The equivalent section-based source, for convenient editing.

verification/
    All fixed data, verification programs, the component index,
    SHA-256 input manifest, and fresh results for the 30 principal components.

requirements.txt
    The Python dependency versions used for the supplied computations.

BUILDING THE ARTICLE

With a normal TeX Live or MiKTeX installation, run from this directory:

    latexmk -pdf -interaction=nonstopmode -halt-on-error Schmeisser_6251.tex

Alternatively run pdflatex twice, or until cross-references stabilize.
The article uses the standard amsart class and the packages named in
its preamble. No BibTeX or external download is part of compilation.

REPRODUCING THE FINITE VERIFICATIONS

The supplied computations used Python 3.12.14, python-flint 0.9.0, and
SymPy 1.14.0. Install the pinned dependencies in a suitable environment:

    python3 -m pip install -r requirements.txt

List every executable component:

    python3 verification/run_all.py --list

Run the 30 principal components:

    python3 verification/run_all.py --jobs 3

Run all 45 components, including 15 additional implementations and
independent targeted cross-checks:

    python3 verification/run_all.py --extended --jobs 3

Run selected components, or preserve the supplied results separately:

    python3 verification/run_all.py --select origin_t,origin_x --jobs 2
    python3 verification/run_all.py --extended --jobs 3 --output fresh_results

The runner verifies every fixed input hash before execution. Each
component then runs in its own temporary copy of the fixed files.
The temporary directories contain the fixed inputs listed in the manifest.
New logs and receipts are written
to verification/results/ unless --output specifies another directory.
SUMMARY.json records the outcome and runtime environment.

The directory hierarchy under verification/source/ is retained because
some independent programs resolve their fixed data by relative paths.
It is included in this portable package.

THE ROLE OF THE ELECTRONIC APPENDIX

The article gives the analytic arguments, parameter domains, complete
tails, closed-boundary comparisons, coverage rules, and uniform degree
monotonicities. The electronic appendix supplies the larger finite
rational lists and executable checks. Together they provide the analytic proof and its finite verifications.

Exact polynomial cancellations use integer or rational arithmetic.
Polynomial inequalities are established on whole intervals or boxes,
using rational Bernstein coefficients or Sturm chains. Transcendental
inequalities use explicit complete series remainders or directed
Arb/Acb balls; an unresolved sign is rejected. No numerical optimizer
or sampled plot is an acceptance criterion.

INDEX.json identifies the checker for each component. MANIFEST.json
identifies the fixed inputs by SHA-256. Individual receipts provide
the finite counts and numerical enclosures produced by each checker.
Displayed decimal diagnostics do not replace the underlying exact or
outward-rounded comparisons.
