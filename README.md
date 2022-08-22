# GMatElastoPlastic

[![CI](https://github.com/tdegeus/GMatElastoPlastic/workflows/CI/badge.svg)](https://github.com/tdegeus/GMatElastoPlastic/actions)
[![Doxygen -> gh-pages](https://github.com/tdegeus/GMatElastoPlastic/workflows/gh-pages/badge.svg)](https://tdegeus.github.io/GMatElastoPlastic)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/gmatelastoplastic.svg)](https://anaconda.org/conda-forge/gmatelastoplastic)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/python-gmatelastoplastic.svg)](https://anaconda.org/conda-forge/python-gmatelastoplastic)

Elasto-plastic material model.
An overview of the theory can be found in `docs/theory/readme.tex`
conveniently compiled to this [PDF](docs/theory/readme.pdf).

# Disclaimer

This library is free to use under the
[MIT license](https://github.com/tdegeus/GMatElastoPlastic/blob/master/LICENSE).
Any additions are very much appreciated, in terms of suggested functionality, code,
documentation, testimonials, word-of-mouth advertisement, etc.
Bug reports or feature requests can be filed on
[GitHub](https://github.com/tdegeus/GMatElastoPlastic).
As always, the code comes with no guarantee.
None of the developers can be held responsible for possible mistakes.

Download:
[.zip file](https://github.com/tdegeus/GMatElastoPlastic/zipball/master) |
[.tar.gz file](https://github.com/tdegeus/GMatElastoPlastic/tarball/master).

(c - [MIT](https://github.com/tdegeus/GMatElastoPlastic/blob/master/LICENSE))
T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me |
[github.com/tdegeus/GMatElastoPlastic](https://github.com/tdegeus/GMatElastoPlastic)

# Python implementation

## Partial example

```python
import GMatElastoPlastic.Cartesian3d as GMat

shape = [...]
K = np.empty(shape)
G = np.empty(shape)
sigy0 = np.empty(shape)
H = np.empty(shape)

GMat.LinearHarderningXd model(K, G, sigy0, H)
...

Eps = np.empty(shape + [3, 3])
...

model.Eps = Eps
print(model.Sig)

...
model.increment()
```

## Installation

### Using conda

```bash
conda install -c conda-forge python-gmatelastoplastic
```

Note that *xsimd* and hardware optimisations are **not enabled**.
To enable them you have to compile on your system, as is discussed next.

### From source

>   You need *xtensor*, *xtensor-python* and optionally *xsimd* as prerequisites.
>   The easiest is to use *conda* to get the prerequisites:
>
>   ```bash
>   conda install -c conda-forge xtensor xsimd xtensor-python
>   ```
>
>   If you then compile and install with the same environment you should be good to go.
>   Otherwise, a bit of manual labour might be needed to treat the dependencies.

```bash
git checkout https://github.com/tdegeus/GMatElastoPlastic.git
cd GMatElastoPlastic

# Only if you want to use hardware optimisation:
export SKBUILD_CONFIGURE_OPTIONS="-DUSE_SIMD=1"

python -m pip install . -v
```

# C++ implementation

## Partial example

```cpp
#include <GMatElastoPlastic.Cartesian3d.h>

namespace GMat = GMatElastoPlastic::Cartesian3d;

int main()
{
    static const size_t rank = ...;

    xt::xtensor<double, rank> K = ...;
    xt::xtensor<double, rank> G = ...;
    xt::xtensor<double, rank> sigy0 = ...;
    xt::xtensor<double, rank> H = ...;

    GMat::LinearHarderningXd model(K, G, sigy0, H);
    ...

    xt::xtensor<double, rank + 2> Eps;
    ...

    // all necessary computation are done at this point
    model.set_Eps(Eps);
    ...

    // get reference to stress
    auto Sig = model.Sig();

    // update history varibles for the next indcrement
    model.increment();

    return 0;
}
```

## Installation

### Using conda

```bash
conda install -c conda-forge gmatelastoplastic
```

### From source

```bash
git checkout https://github.com/tdegeus/GMatElastoPlastic.git
cd GMatElastoPlastic

cmake -Bbuild
cd build
cmake --install .
```

## Compiling

## Using CMake

### Example

Your `CMakeLists.txt` can be as follows

```cmake
cmake_minimum_required(VERSION 3.1)
project(example)
find_package(GMatElastoPlastic REQUIRED)
add_executable(example example.cpp)
target_link_libraries(example PRIVATE GMatElastoPlastic)
```

### Targets

The following targets are available:

*   `GMatElastoPlastic`
    Includes the library and its dependencies.

*   `GMatElastoPlastic::assert`
    Enables IO-assertions by defining `GMATELASTOPLASTIC_ENABLE_ASSERT`.

*   `GMatElastoPlastic::debug`
    Enables assertions of all dependencies.

*   `GMatElastoPlastic::compiler_warings`
    Enables compiler warnings (generic).

### Optimisation

It is advised to think about compiler optimisation and enabling *xsimd*.
Using *CMake* this can be done using the `xtensor::optimize` and `xtensor::use_xsimd` targets.
The above example then becomes:

```cmake
cmake_minimum_required(VERSION 3.1)
project(example)
find_package(GMatElastoPlastic REQUIRED)
find_package(xtensor REQUIRED)
find_package(xsimd REQUIRED)
add_executable(example example.cpp)
target_link_libraries(example PRIVATE
    GMatElastoPlastic
    xtensor::optimize
    xtensor::use_xsimd)
```

See the [documentation of xtensor](https://xtensor.readthedocs.io/en/latest/).

## By hand

Presuming that the compiler is `c++`, compile using:

```
c++ -I/path/to/GMatElastoPlastic/include ...
```

Note that you have to take care of the *xtensor* dependency, the C++ version, optimisation,
enabling *xsimd*, ...

## Using pkg-config

Presuming that the compiler is `c++`, compile using:

```
c++ `pkg-config --cflags GMatElastoPlastic` ...
```

Note that you have to take care of the *xtensor* dependency, the C++ version, optimization,
enabling *xsimd*, ...

# Upgrading instructions

## Upgrading to >v0.3.*

The individual material point and the array of material points was fully integrated.
In addition, the number of copies was reduced.

### C++

There is only a single class `LinearHarderning`. It's functions where renamed:

*   `.setStrain(...)` -> `.set_Eps(...)`
*   `.Stress()` -> `.Sig()` (now returns a reference).
*   `.stress(...)`: deprecated.
*   `.Tangent()` -> `.C()` (now returns a reference).
*   `.tangent(...)`: deprecated.

### Python

There is only a single class `LinearHarderning`. It's functions are converted to properties:

*   `.setStrain(...)` -> `.Eps = ...`
*   `.Stress()` -> `.Sig` (now returns a reference).
*   `.stress(...)`: deprecated.
*   `.Tangent()` -> `.C` (now returns a reference).
*   `.tangent(...)`: deprecated.

## Upgrading to >v0.2.*

`xtensor_fixed` was completely deprecated in v0.2.0, as were the type aliases
`Tensor2` and `Tensor4`.
Please update your code as follows:

*   `Tensor2` -> `xt::xtensor<double, 2>`.
*   `Tensor4` -> `xt::xtensor<double, 4>`.

**Tip:** Used `auto` as return type as much as possible.
This simplifies implementation, and renders is less subjective to library
return type changes.

Compared to v0.1.0, v0.2.0 has some generalisations and efficiency updates.
This requires the following changes:

*   `Matrix` has been generalised to `Array<rank>`. Practically this requires changing:
    -   `Matrix` to `Array<2>` in C++.
    -   `Matrix` to `Array2d` in Python.
        Note that `Array1d`, `Array3d`, are also available.

*   `Array<rank>.check` ->
    ```cpp
    if (xt::any(xt::equal(array.type(), Type::Unset))) {
        throw std::runtime_error("Please set all points");
    }
    ```
    Note however that it is no longer required to set all points,
    unset points are filled-up with zeros.

*   Strain is now stored as a member.
    Functions like `stress` now return the state based on the last specified strain,
    specified using `setStrain(Esp)`. This leads to the following changes:
    - `stress`: no argument.
    - `tangent`: no argument, single return value (no longer returns stress).

# Change-log

## v0.3.0

Complete API overhaul.

## v0.2.1

*   Using scikit-build, setuptools_scm, xtensor-python (#21)
*   CMake clean-up (#21)

## v0.2.0

Compared to v0.1.0, v0.2.0 has some generalisations and efficiency updates.
This requires the following changes:

*   `Matrix` has been generalised to `Array<rank>`. Practically this requires changing:
    -   `Matrix` to `Array<2>` in C++.
    -   `Matrix` to `Array2d` in Python.
        Note that `Array1d`, `Array3d`, are also available.

*   `Array` now sets zeros for all `Type::Unset` points.
    The function `check` is deprecated accordingly.

*   Strain is now stored as a member.
    Functions like `stress` now return the state based on the last specified strain,
    specified using `setStrain(Esp)`. This leads to the following changes:
    - `stress`: no argument.
    - `tangent`: no argument, single return value (no longer returns stress).

*   Tensor operations are now provided centrally in the GMat eco-system,
    by GMatTensor
