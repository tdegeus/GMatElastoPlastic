/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTOPLASTIC_CONFIG_H
#define GMATELASTOPLASTIC_CONFIG_H

/**
All assertions are implementation as:

    GMATELASTOPLASTIC_ASSERT(...)

They can be enabled by:

    #define GMATELASTOPLASTIC_ENABLE_ASSERT

(before including GMatElastic).
The advantage is that:

-   File and line-number are displayed if the assertion fails.
-   Assertions can be enabled/disabled independently from those of other libraries.

\throw std::runtime_error
*/
#ifdef GMATELASTOPLASTIC_ENABLE_ASSERT
#define GMATELASTOPLASTIC_ASSERT(expr) GMATTENSOR_ASSERT_IMPL(expr, __FILE__, __LINE__)
#else
#define GMATELASTOPLASTIC_ASSERT(expr)
#endif

/**
Linear elastoplastic material model.
*/
namespace GMatElastoPlastic {

/**
Define container type.
The default `xt::xtensor` can be changed using:

-   `#define GMATELASTOPLASTIC_USE_XTENSOR_PYTHON` -> `xt::pytensor`
*/
namespace array_type {

#ifdef GMATELASTOPLASTIC_USE_XTENSOR_PYTHON

/**
Fixed (static) rank array.
*/
template <typename T, size_t N>
using tensor = xt::pytensor<T, N>;

#else

/**
Fixed (static) rank array.
*/
template <typename T, size_t N>
using tensor = xt::xtensor<T, N>;

#endif

} // namespace array_type

} // namespace GMatElastoPlastic

#endif
