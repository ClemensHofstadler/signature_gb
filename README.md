# SignatureGB_SageMath

## DESCRIPTION

A beta version of the SageMath package SignatureGB for computing signature Gröbner bases in the free algebra.
An actual release with documentation, user guide and examples will be coming soon.

## LICENCSE

Distributed under the terms of the GNU General Public License, either version 2 or (at your option) any later version (https://www.gnu.org/licenses/)

## REQUIREMENTS

- SageMath 9.1 or later

- The Python library pyahocorasick (https://pyahocorasick.readthedocs.io/en/latest/)

## INSTALLATION

Download the source code from this page.

Then install the Python library pyahocorasick as decribed [here](https://pyahocorasick.readthedocs.io/en/latest/).

After this, the SignatureGB package can be installed. As this package is written mostly in Cython,
it needs to be compiled before it can be used. This can be done by running the command
```
sage setup.py build_ext --inplace
```
inside the folder SignatureGB.

## USAGE

More information on how to use the package will be provided here soon together with the first actual release.

Note: The benchmark examples included in the papers *Short proofs of ideal membership* and *Signature Gröbner bases in the free algebras over rings* can be run using this beta version. For further information, we refer to the respective README files that come together with the benchmark examples.

In the following, we describe some basic functionality.

Make sure the path to the folder SignatureGB is visible to SageMath, for example, by running
```
  sys.path.append(PATH)
```
where PATH is the path to the directory in which the folder SignatureGB is.

Then the package can be loaded by calling
```
  from SignatureGB import *
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
Similarly, a nested list $[[x_1,\dots,x_{n_1}],[y_1,\dots,y_{n_2}],\dots,[z_1,\dots,z_{n_k}]]$ yields a block ordering
$x_1 < \dots x_{n_1} \ll y_1 < \dots y_{n_2} \ll \dots \ll z_1 < \dots < z_{n_k}$, where each block is ordered degree-lexicographically.
For example
```
  LabelledModule(gens,[[a,b],[c,d]])
```
yields a LabelledModule with monomial order $a < b \ll c < d$.
By default, the signature ordering is degree-over-position-over-term (dpot).
It can be changed to degree-over-term-over-position (dtop) as follows.
```
  LabelledModule(gens,[[a,b],[c,d]],signature_order='dtop')
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

One a labelled Gröbner basis is reconstructed, a LabelledModule also provides the possibility to test ideal membership of noncommutative polynomials.
If an ideal membership can be verified, it outputs a cofactor representation.
```
  F.<a,b,c,d,e> = FreeAlgebra(QQ)
  gens = [1-a*b,1-b*a,a*e*a-a, e*a*e-e, a*e-c*d, e*a-d*c, c*d*c-c, d*c*d-d]
  M = LabelledModule(gens,[a,b,c,d,e])
  G,H = M.signature_GB(100)
  M.reconstruct_labelled_basis()
  M.membership_test(b - e)
```

