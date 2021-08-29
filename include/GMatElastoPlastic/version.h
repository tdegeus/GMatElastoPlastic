/**
Version information.

\file
\copyright Copyright 2018. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTOPLASTIC_VERSION_H
#define GMATELASTOPLASTIC_VERSION_H

#include "config.h"
#include <GMatElastic/version.h>
#include <GMatTensor/version.h>

/**
Current version.

Either:

-   Configure using CMake at install time. Internally uses::

        python -c "from setuptools_scm import get_version; print(get_version())"

-   Define externally using::

        -DGMATELASTOPLASTIC_VERSION="`python -c "from setuptools_scm import get_version; print(get_version())"`"

    From the root of this project. This is what ``setup.py`` does.

Note that both ``CMakeLists.txt`` and ``setup.py`` will construct the version using ``setuptools_scm``.
Tip: use the environment variable ``SETUPTOOLS_SCM_PRETEND_VERSION``
to overwrite the automatic version.
*/
#ifndef GMATELASTOPLASTIC_VERSION
#define GMATELASTOPLASTIC_VERSION "@PROJECT_VERSION@"
#endif

namespace GMatElastoPlastic {

/**
Return version string, e.g.

    "0.8.0"

\return std::string
*/
inline std::string version();

/**
Return versions of this library and of all of its dependencies.
The output is a list of strings:

    "gmatelastic=0.2.0",
    "gmatelastoplastic=0.7.0",
    "gmattensor=0.8.0",
    "xtensor=0.20.1"

\return List of strings.
*/
inline std::vector<std::string> version_dependencies();

} // namespace GMatElastoPlastic

#include "version.hpp"

#endif
