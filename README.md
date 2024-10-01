# siganture_gb

## DESCRIPTION

SageMath package for signature-based Gröbner basis computations in the free algebra.

## LICENCSE

Distributed under the terms of the GNU General Public License, either version 2 or (at your option) any later version (https://www.gnu.org/licenses/)

## REQUIREMENTS

- SageMath 9.3 or later is recommended

## DEPENDENCIES

- The Python library `pyahocorasick` (https://pyahocorasick.readthedocs.io/en/latest/).
  Will be installed automatically when `signature_gb` is installed as described below.
  

## INSTALLATION

### With SageMath built from source or binaries from sagemath.org

**Note**: This way of installing the package also automatically installs the `pyahocorasick` library.
Thus, except executing the command below, no additional work is required.

To download and install the latest version on a system where SageMath
was built from source or installed from official packages, run

    sage -pip install [--user] git+https://github.com/ClemensHofstadler/signature_gb.git
    
The optional `--user` flag causes the package to be installed in your `.sage`
directory instead of the SageMath installation tree.

Alternatively, run (square brackets indicate optional flags)

    sage -pip install [--user] [--editable] .

from the root of a local git checkout. The `--editable` flag causes the
"installed" version to point to your local checkout, making it easier,
if you edit the code, to run the modified version. See the pip documentation
for more installation options.

Microsoft Windows users should run the above commands in a "SageMath shell", see

- https://wiki.sagemath.org/SageWindows

Apple macOS users may need additional steps before they are able to add external
packages to their Sage installations. See

- https://github.com/3-manifolds/fix_mac_sage/releases
- https://ask.sagemath.org/question/51130

for more information.

### Using the package without installation

To use `signature_gb` directly from a git checkout (without installation), run

    AHOCORASICK_BYTES=yes sage -pip install ./pyahocorasick-master
from the checkout, followed by

    sage -python setup.py build_ext --inplace

After that, add the `src/` directory to your Python `sys.path`.

The package contains compiled (Cython) modules which are automatically built as
part of the installation procedure. Installation will fail if they cannot be
built.


## USAGE

More information on how to use the package will be provided here soon together with the first actual release.

Note: The benchmark examples included in the papers *Short proofs of ideal membership* and *Signature Gröbner bases in the free algebras over rings* can be run using this beta version. For further information, we refer to the respective README files that come together with the benchmark examples.

In the following, we describe some basic functionality.

Make sure the path to the folder signature_gb is visible to SageMath, for example, by running
```
  sys.path.append(PATH)
```
where PATH is the path to the directory in which the folder signature_gb is.

Then the package can be loaded by calling
```
  from signature_gb import *
```
The basic data strcuture provided by the package is that of a *LabelledModule*.
It can be constructed from a list of SageMath noncommutative polynomials as follows.

*Remark*: So far, only polynomials with rational coefficients are supported.
Furthermore, variable names are restricted to single (lower- and uppercase) characters.

```
  F.<a,b,c,d> = FreeAlgebra(QQ)
  gens = [a*b*a-a, b*a*b-b, a*b-c*d, b*a-d*c, c*d*c-c, d*c*d-d]
  M = LabelledModule(gens,[a,b,c,d])
```
The second argument provides the monomial order w.r.t. which the computations will be exectued.
A list $[x_1,x_2,\dots, x_n]$ yields a degree-lexicographic ordering $x_1 < x_2 < \dots x_n$.
By default, the signature ordering is degree-over-position-over-term (dpot).
It can be changed to degree-over-term-over-position (dtop) as follows.
```
  LabelledModule(gens,[a,b,c,d],signature_order='dtop')
```

The main functionality provided by a LabelledModule is computing signature and labelled Gröbner bases.
A signature Gröbner basis of a LabelledModule can be enumerated as follows.
```
  M = LabelledModule(gens,[a,b,c,d])
  G, H = M.signature_GB(100)
```
This runs 100 iterations of a signature-based algorithm and outputs a (partial) signature Gröbner basis *G* and a (partial) basis *H* of the leading term module of the syzgy module.
To compute a signature basis up to some fixed signature, a *sig_bound* can be provided in form of a positive integer $N$.
Then the algorithm computes a signature basis up to degree $N$ (if the number of iterations is chosen large enough).

```
  M = LabelledModule(gens,[a,b,c,d])
  G, H = M.signature_GB(100,sig_bound=3)
```
To reconstruct a (partial) labelled Gröbner basis and a (partial) basis of the syzygy module, run the following commands in the given order.

```
  G,H = M.signature_GB(100)
  G2 = M.reconstruct_labelled_basis()
  H2 = M.reconstruct_syzygies()
```

Once a labelled Gröbner basis is reconstructed, a LabelledModule also provides the possibility to test ideal membership of noncommutative polynomials.
If an ideal membership can be verified, it outputs a cofactor representation.
```
  F.<a,b,c,d,e> = FreeAlgebra(QQ)
  gens = [1-a*b,1-b*a,a*e*a-a, e*a*e-e, a*e-c*d, e*a-d*c, c*d*c-c, d*c*d-d]
  M = LabelledModule(gens,[a,b,c,d,e])
  G,H = M.signature_GB(100)
  M.reconstruct_labelled_basis()
  M.membership_test(b - e)
```

