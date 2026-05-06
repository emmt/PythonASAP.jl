# PythonASAP [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://emmt.github.io/PythonASAP.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://emmt.github.io/PythonASAP.jl/dev/) [![Build Status](https://github.com/emmt/PythonASAP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/PythonASAP.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/emmt/PythonASAP.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/PythonASAP.jl)

`PythonASAP` is a thin layer over [`ASAP`](https://github.com/emmt/ASAP.jl) to facilitate the
use of `ASAP` from Python.

# Installation

On Julia side, `PythonASAP` requires:

- Julia (of course), [`juliaup`](https://github.com/JuliaLang/juliaup) is the recommended
  method to do this;

- a few packages, among others [`ASAP`](https://github.com/emmt/ASAP.jl) and
  [`PythonCall`](https://github.com/JuliaPy/PythonCall.jl).

In the following, we only assume that Julia has been installed and that the Git repositories
of [`ASAP`](https://github.com/emmt/ASAP.jl) and
[`PythonASAP`](https://github.com/emmt/PythonASAP.jl) have been cloned, all these on the
same machine.

On Python side (i.e, to call Julia from Python), the `juliacall` package is needed. This
package may be installed in a *virtual environment* and, as an illustration, we show below
the steps to install `juliacall` in a virtual environment on a Unix machine:

1. Make sure Python package `venv` is installed. On Ubuntu-based systems, this is done by:

   ``` shell
   apt install python3.12-venv
   ```

2. Create a new virtual Python environment, for example:

   ``` shell
   python3 -m venv juliacall
   ```

   will create a virtual environment in directory `juliacall` of the current directory.

3. Activate the virtual Python environment. Assuming Bash, this is done by:

   ``` shell
   . ./juliacall/bin/activate
   ```

   Your prompt changes to include the name of the environment. To exit form the environment,
   type `deactivate` (this is an alias in your shell).

4. While the virtual Python environment is active, install needed `pip` packages and start
   the Python interpreter:

   ``` shell
   pip3 install juliacall # to be done only once for this environment
   ipython3
   ```

5. To check that Julia can be called from Python, start a Python interpreter and type:

   ``` python
   from juliacall import Main as jl
   jl.seval("println(\"Hello world!\")")
   ```

6. To isolate things in Julia, it is good to work within a specific environment. The
   following lines will create such an environment (called `asap+python` in our example with
   settings saved in a project and manifest files in a sub-directory named `asap+python`),
   add dependencies for Julia (assuming that Git repositories for ASAP and PythonASAP are in
   local directories named `ASAP.jl` and `PythonASAP`), precompile, and load the packages.
   All commands are to be executed from the Python interpreter:

   ``` python
   from juliacall import Main as jl
   jl.seval("using Pkg") # we will use Julia's package manager
   jl.Pkg.Registry.add("General") # add general registry of Julia's packages
   jl.Pkg.Registry.add(jl.Pkg.RegistrySpec(url = "https://github.com/emmt/EmmtRegistry")) # add Eric's registry
   jl.Pkg.add("Revise") # add the Revise package
   jl.Pkg.activate("asap+python") # create and activate Julia environment
   jl.Pkg.develop(path="./ASAP.jl") # add the ASAP package for a local Git repository
   jl.Pkg.develop(path="./PythonASAP.jl") # idem for PythonASAPrepository
   jl.Pkg.resolve() # resolve dependencies and precompile
   jl.seval("using Revise, ASAP, PythonASAP") # load packages
   ```

# Usage

The following sub-sections demonstrate the usage of ASAP in Python.


## Starting the session

After the above installation steps, it is sufficient to do the following:

1. In the shell, activate Python virtual environment and start the interpreter:

   ``` shell
   . ./juliacall/bin/activate
   ipython3
   ```

2. In the Python interpreter, import `juliacall` package, activate Julia's environment, and load
   Julia's packages:

   ``` python
   from juliacall import Main as jl
   jl.seval("using Pkg") # load Julia's package manager
   jl.Pkg.activate("asap+python") # activate Julia environment
   jl.seval("using Revise, ASAP, PythonASAP") # load packages, Revise must come first
   ```

All these should be done in the `__init__.py` script of a Python package (or module) that
does not yet exist...

## Creating an ASAP sparse structure

To create an ASAP sparse structure the inputs are:

- `fmt` the format for the triangular sparse matrix: `RowWiseLower`, `RowWiseUpper`,
  `ColumnWiseLower`, or `ColumnWiseUpper` are the 4 possibilities.

- `msk`, a Boolean mask to indicate the size of a data frame and the location of valid nodes
  (where the mask is true).

- `perm`, the strategy to build the permutation: `None` (or `"none"`) to use no
  permutation (lexicographic order of nodes), `"FRiM"` for the original algorithm (require
  mask with all dimensions equal to `2^p + 1` and no invalid nodes), or `"multiscale"` for
  the generalization of FRiM. `perm` may also be specified as a vector of indices to
  directly provide the permutation.

- a parameter to determine the number of structural non-zeros:

  - if a permutation is used, `m`, the maximum number of structural non-zeros per row for a
    row-wise format, per column for a column-wise format;

  - if no permutation is used, `dmax` the maximal distance of the preceding nodes that are
    linked to a given node;

- the multi-scale strategy also requires to indicate the first valid node to start the
  permutation.


``` python
# Setup.
import numpy as np
from juliacall import Main as jl, convert as jlconvert
jl.sevals("using Revise, ASAP, PythonASAP")
jl_asap = jl.PythonASAP # alias

# Create the mask.
dims = (300, 200)
mask = np.ones(dims, np.bool_)

# Build the sparse structure.
S = jl_asap.build("RowWiseLower", mask, "multiscale", 5, jl.CartesianIndex((150,100)))
```


## Learning the coefficients of an ASAP sparse factor

Once the structure of the ASAP sparse factor has been built (the object `S` in the preceding
example), the structural non-zeros may be learned from a given symmetric positive-definite
matrix or from a sample of data frames. In the latter case, the ASAP model will be learned
so as to approximate the empirical covariance matrix of the sample.

For the fast (and recommended) learning method, there are two possibilities:

``` python
# Learn coefficients.
mdl1 = jl_asap.learn("inv(R'*R)", S, dataset)
mdl2 = jl_asap.learn("inv(R*R')", S, dataset)
```

The model `mdl1` (resp. `m2`) corresponds to `M ≈ inverse(transpose(R)⋅R)` (resp. `M ≈
inverse(R⋅transpose(R))`) with `M` the target matrix and `R` a sparse triangular matrix
whose structure is that of `S`. Argument `dataset` is a sequence of data frames.


# Known issues

## Storage order of dimensions

Multi-dimensional `Numpy` arrays have versatile storage order: C-like (row-major),
FORTRAN-like (column-major), or even mixed orders. The storage order is specified by the
array strides. Furthermore array elements are not necessarily contiguous. In `ASAP`, arrays
must have fast linear indexing with 1-based indices and contiguous elements. This is
required by the way permutations are represented (as an ordered list of linear indices).
This may change in the future but until then, to avoid copies, multi-dimensional `Numpy`
arrays must have contiguous elements stored in either column-major or row-major order. In
the latter case, the array dimensions are permuted. For these reasons, `PyArray` objects,
the Julia counterparts of Python arrays provided by `PythonCall`, are converted into
instances of `FastNumPyArray`. This solution is acceptable if all `NumPy` arrays used with
ASAP have the same storage order (this constraint is not asserted in the current code). In
the future, better solutions may be devised perhaps involving temporary copies.

## Multi-threading

`ASAP` may benefit from multi-threading to accelerate some computations like learning the
coefficients of the approximation.
