# signature_gb

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

To use `signature_gb` directly from a git checkout (only installing `pyahocorasick` but not `signature_gb`), run

    AHOCORASICK_BYTES=yes sage -pip install ./pyahocorasick-master
    
from the checkout, followed by

    sage -python setup.py build_ext --inplace

After that, add the `src/` directory to your Python `sys.path`.

The package contains compiled (Cython) modules which are automatically built as
part of the installation procedure. Installation will fail if they cannot be
built.


## USAGE

An actual documentation on how to use the package will be provided here soon.

In the following, we describe some basic functionality.

After installing the package as described above, its functionality can be 
imported into a Sage session as follows
```
  from signature_gb import *
```
The basic data strcuture provided by the package is that of a `LabelledModule`.
It can be constructed from a list of SageMath noncommutative polynomials as follows.

*Remark*: So far, only polynomials with rational coefficients are supported.

```
  F.<a,b,c,d> = FreeAlgebra(QQ)
  gens = [a*b*a-a, b*a*b-b, a*b-c*d, b*a-d*c, c*d*c-c, d*c*d-d]
  M = LabelledModule(gens)
```
A `LabelledModule` comes together with a monomial ordering and a module ordering for the signature-based
computations. By default, the monomial ordering is inherited from the parent of the generators, i.e.,
from the free algebra. The default module ordering is degree-over-position-over-term (DoPoT).

To change the monomial ordering, an optional argument `monomial_order` can be provided.
This argument can be set to a list of variables `monomial_order = [x1,x2,...,xn]`, 
which yields a degree-lexicographic ordering $x_1 < x_2 < \dots x_n$.
```
  LabelledModule(gens, monomial_order=[d,c,b,a])
```
It can also be set to a list of list of variables `monomial_order = [[x1,...,xi],...,[y1,...,yj]]`,
which yields a block-ordering $x_1 < \dots x_i \ll \dots \ll y_1 < \dots < y_j$ (each block is compared with a degree-lexicographic ordering).
```
  LabelledModule(gens, monomial_order=[[a,b],[c,d]])
```
To change the module ordering, the optional argument `signature_order` can be used.
By default, this is set to `signature_order = 'DoPoT'`. It can be changed to a 
degree-over-term-over-position (DoToP) as follows
```
  LabelledModule(gens, monomial_order=[[a,b],[c,d]], signature_order='DoToP')
```

The main functionality provided by a `LabelledModule` is computing signature Gröbner bases.
A signature Gröbner basis of a `LabelledModule` can be enumerated via the method `signature_basis`.
The following optional arguments can be provided to this method:
- `maxiter` (default: `10`): The maximal number of iterations of a signature-based Gröbner basis algorithm
- `maxdeg` (default: `-1`): The maximal degree of ambiguities that are considered during the computation.
  The default value `-1` causes all ambiguities to be considered.
- `sig_bound` (default: `-1`): The maximal degree of signatures that are considered during the computation.
  The default value `-1` causes all signatures to be considered.
- `verbose` (default: `0`): Determines how much information about the computational progress is displayed.
  The higher the value, the more information is printed.

So, for example, the command
```
  G, H = M.signature_basis(maxiter=100,sig_bound=3)
```
runs at most 100 iterations of a signature-based Gröbner basis algorithm and outputs a signature Gröbner basis
`G` up to degree 3 (i.e., the result is a signature Gröbner basis up to the smallest signature with degree 3)
as well as a (partial) basis `H` of the leading term module of the syzgy module.

The package also provides a modular algorithm for computing signature Gröbner bases.
It can be called via the method `modular_signature_basis`, which takes the following optional arguments:
- all optional arguments that can also be given to `signature_gb` with the same meaning.
- `threads` (default: `None`): Number of threads to be used for the parallel computations. If the default value
  `None` is used, as many threads will be used as the method `os.cpu_count()` returns.
- `num_primes` (default: `None`): Number of primes to be used before the first reconstruction attempt.
  If the default value `None` is used, `num_primes` is set to the maximimum of 8 and `threads`.
  If `num_primes` is not enough to reconstruct a basis, this value is increased by `threads` until a successful reconstruction is possible.
- `verification` (default: `'rigorous'`): Sets the verification procedure. By default, a rigouros verification over the rationals is done.
  If set to `'probabilistic'`, a verification in positive characteristic is done, which is usually faster but only yields the correct result with
  high probability. If set to any other value, no verification is performed at all.
```
  G,H = M.modular_signature_basis(maxiter=100, threads=4, num_primes=4)
```

## More examples

To compute, for instance, a signature Gröbner basis up to degree 9 for the benchmark `cyclic5` using the modular algorithm, we can proceed as follows:

```
from signature_gb import *

F = FreeAlgebra(QQ,['v','w','x','y','z'])
v,w,x,y,z = F.gens()
gens = [ 
    v+w+x+y+z,
    v*w+w*x+x*y+v*z+y*z,
    v*w*x+w*x*y+v*w*z+v*y*z+x*y*z,
    v*w*x*y+v*w*x*z+v*w*y*z+v*x*y*z+w*x*y*z,
    v*w*x*y*z-1,
    ]
cyclic5 = LabelledModule(gens)
G, H = cyclic5.modular_signature_basis(maxiter=100,sig_bound=9,threads=2)
```

Using more threads and changing the verification to `probabilistic` speeds up the computation:
```
G, H = cyclic5.modular_signature_basis(maxiter=100,sig_bound=9,threads=4,verification="probabilistic")
```

Here is another example using the `P6` benchmark. We use 5 primes before the first reconstruction attempt.
```
from signature_gb import *

F = FreeAlgebra(QQ,['a','b','c'])
a,b,c = F.gens()
gens = [
    c*c*c + 2*c*c*b + 3*c*c*a + 5*b*c*c + 7*a*c*a, 
    b*c*c + 11*b*a*b + 13*a*a*a
    ]
P6 = LabelledModule(gens)
G, H = P6.modular_signature_basis(maxiter=50,sig_bound=10, num_primes=5, verbose=1)
```
The option `verbose=1` prints some information about the computational progress (in particular, it shows that 5 primes are not sufficient to correctly reconstruct the signature basis)


