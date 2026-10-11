# Online Resource 2: exact computational certificates

Global certificates for endogenous mortgage entry under joint collateral-law constraints
Journal of Global Optimization
Author and corresponding author: Letong Zhou
Affiliation: EDHEC Business School, Paris, France
Email: letong.zhou@edhec.com

The reader guide is Online Resource 1 (ESM_1.pdf), supplied separately.
This archive contains the exact synthetic inputs embedded in Python sources,
lossless rational results and the routines that reproduce them.

Install dependencies and run ten suites under normal and optimized Python:

    python3 -m pip install -r requirements.txt
    python3 reproduce.py --math-only

The mathematical routines require no external data, earlier report or TeX
installation. The conservative version floor is Python 3.11. Python 3.12.14,
SymPy 1.14.0 and mpmath 1.3.0 were used for the preparation checks; the
compatibility audit also ran clean environments with Python 3.11.16 and
3.13.5 and the same pinned dependencies. The expected canonical hashes
include every report field. The kernel report retains legacy traceback
frames, the exact failure class and message, exit status and absence of a
success report; it removes interpreter-only range underlines and the added
source line of the -c command. The raw trace remains in the execution log.
SOURCE_MANIFEST.json records the shipped checker hashes. Reproduction creates
its own execution ledger and console records in the extracted directory.

Additional audits:

    python3 validation/coverage_audit.py
    python3 validation/decimal_serializer_regression.py
    python3 -O validation/decimal_serializer_regression.py

Rational decisions use Fraction or SymPy Rational; displayed decimal endpoints
are rounded outward. Parameter cells and type cells enclose their entire
continuous domains. The analytic proofs in the manuscript establish the
unrestricted law/menu/lottery conclusions. The figure uses floating coordinates
only for rendering and is supplied with the editable manuscript.

See IDENTIFICATION.json for resource metadata and FILE_MANIFEST.json for
the included files. The manuscript source and editorial correspondence are
separate submission files.
