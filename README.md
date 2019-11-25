# GMatElastoPlastic

Elasto-plastic material model. An overview of the theory can be found in `docs/` in particular in this [PDF](docs/readme.pdf).

# Contents

<!-- MarkdownTOC -->

- [Implementation](#implementation)
- [Installation](#installation)
    - [Linux / macOS](#linux--macos)
        - [Install systemwide \(depends on your privileges\)](#install-systemwide-depends-on-your-privileges)
        - [Install in custom location \(user\)](#install-in-custom-location-user)

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

## Linux / macOS

### Install systemwide (depends on your privileges)

1.  Proceed to a (temporary) build directory. For example:

    ```bash
    cd /path/to/GMatElastoPlastic
    mkdir build
    cd build
    ```

2.  'Install' `GMatElastoPlastic`. For the path in **1.**:

    ```bash
    cmake .. 
    make install
    ```

> One usually does not need any compiler arguments after following this protocol.

### Install in custom location (user)

1.  Proceed to a (temporary) build directory. For example:

    ```bash
    cd /path/to/GMatElastoPlastic
    mkdir build
    cd build
    ```

2.  'Install' `GMatElastoPlastic`, to install it in a custom location. For the path in **1.**:

    ```bash
    mkdir /custom/install/path
    cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/custom/install/path
    make install
    ```

3.  Add the appropriate paths to for example your ``~/.bashrc`` (or ``~/.zshrc``). For the path in **2.**: 

    ```bash
    export PKG_CONFIG_PATH=/custom/install/path/share/pkgconfig:$PKG_CONFIG_PATH
    export CPLUS_INCLUDE_PATH=$HOME/custom/install/path/include:$CPLUS_INCLUDE_PATH
    ```

> One usually has to inform the CMake or the compiler about `${CPLUS_INCLUDE_PATH}`.

