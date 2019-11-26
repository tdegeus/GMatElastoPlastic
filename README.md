# GMatElastoPlastic

[![Travis](https://travis-ci.com/tdegeus/GMatElastoPlastic.svg?branch=master)](https://travis-ci.com/tdegeus/GMatElastoPlastic)

Elasto-plastic material model. An overview of the theory can be found in `docs/` in particular in this [PDF](docs/readme.pdf).

# Contents

<!-- MarkdownTOC levels="1,2" -->

- [Implementation](#implementation)
- [Installation](#installation)
    - [C++ headers](#c-headers)
    - [Python module](#python-module)
- [Compiling](#compiling)
    - [By hand](#by-hand)
    - [Using pkg-config](#using-pkg-config)
    - [Using `CMakeLists.txt`](#using-cmakeliststxt)

<!-- /MarkdownTOC -->

# Implementation

The headers are meant to be self-explanatory, please check them out:

* [Cartesian3d.h](include/GMatElastoPlastic/Cartesian3d.h)

Only a tiny example is presented here, that is meant to understand the code's structure:

```cpp
#include <GMatElastoPlastic/Cartesian3d.h>

int main()
{
    // a single material point
    // - create class
    GMatElastoPlastic::Cartesian3d::Elastic elastic(K, G);
    GMatElastoPlastic::Cartesian3d::LinearHardening plastic(K, G, sigy0, H);
    // - compute stress [allocate result]
    Sig = elastic.Stress(Eps);
    ...
    // - compute stress [no allocation]
    elastic.stress(Eps, Sig); 
    ...

    // a "matrix" of material points
    // - create class
    GMatElastoPlastic::Cartesian3d::Elastic matrix(nelem, nip);
    // - set material
    matrix.setElastic(I, K, G);
    matrix.setLinearHardening(I, K, G, sigy0, H);
    // - compute stress [allocate result]
    Sig = matrix.Stress(Eps);
    ...
    // - compute stress [no allocation]
    matrix.stress(Eps, Sig); 
    ...
}
```

# Installation

## C++ headers

### Using conda

```bash
conda install -c conda-forge gmatelastoplastic
```

### From source

```bash
# Download GMatElastoPlastic
git checkout https://github.com/tdegeus/GMatElastoPlastic.git
cd GMatElastoPlastic

# Install headers, CMake and pkg-config support
cmake .
make install
```

## Python module

### From source

> To get the prerequisites you *can* use conda
> 
> ```bash
> conda install -c conda-forge pyxtensor
> conda install -c conda-forge xsimd
> ```

```bash
# Download GMatElastoPlastic
git checkout https://github.com/tdegeus/GMatElastoPlastic.git
cd GMatElastoPlastic

# Compile and install the Python module
python setup.py build
python setup.py install
```

# Compiling

## By hand

Presuming that the compiler is `c++`, compile using:

```
c++ -I/path/to/GMatElastoPlastic/include ...
```

## Using pkg-config

Presuming that the compiler is `c++`, compile using:

```
c++ `pkg-config --cflags GMatElastoPlastic` ...
```

## Using `CMakeLists.txt`

Using *GMatElastoPlastic* the `CMakeLists.txt` can be as follows

```cmake
cmake_minimum_required(VERSION 3.1)

project(example)

find_package(xtensor REQUIRED)
find_package(GMatElastoPlastic REQUIRED)

add_executable(example example.cpp)

target_link_libraries(example
    PRIVATE
    xtensor
    GMatElastoPlastic)
```

Compilation can then proceed using 

```bash
cmake .
make
```
