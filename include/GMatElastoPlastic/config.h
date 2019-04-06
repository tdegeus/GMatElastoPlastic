/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot

================================================================================================= */

#ifndef GMATELASTOPLASTIC_CONFIG_H
#define GMATELASTOPLASTIC_CONFIG_H

// -------------------------------------------------------------------------------------------------

#include <tuple>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xmath.hpp>

// -------------------------------------------------------------------------------------------------

#ifndef NDEBUG
#define GMATELASTOPLASTIC_ENABLE_ASSERT
#endif

#ifdef GMATELASTOPLASTIC_ENABLE_ASSERT
#define GMATELASTOPLASTIC_ASSERT(expr) GMATELASTOPLASTIC_ASSERT_IMPL(expr, __FILE__, __LINE__)
#define GMATELASTOPLASTIC_ASSERT_IMPL(expr, file, line)                                                                   \
    if (!(expr))                                                                                                          \
    {                                                                                                                     \
        throw std::runtime_error(std::string(file) + ':' + std::to_string(line) + ": assertion failed (" #expr ") \n\t"); \
    }
#else
#define GMATELASTOPLASTIC_ASSERT(expr)
#endif

// -------------------------------------------------------------------------------------------------

#define GMATELASTOPLASTIC_WORLD_VERSION 0
#define GMATELASTOPLASTIC_MAJOR_VERSION 0
#define GMATELASTOPLASTIC_MINOR_VERSION 1

#define GMATELASTOPLASTIC_VERSION_AT_LEAST(x,y,z) \
  (GMATELASTOPLASTIC_WORLD_VERSION>x || (GMATELASTOPLASTIC_WORLD_VERSION>=x && \
  (GMATELASTOPLASTIC_MAJOR_VERSION>y || (GMATELASTOPLASTIC_MAJOR_VERSION>=y && \
                                         GMATELASTOPLASTIC_MINOR_VERSION>=z))))

#define GMATELASTOPLASTIC_VERSION(x,y,z) \
  (GMATELASTOPLASTIC_WORLD_VERSION==x && \
   GMATELASTOPLASTIC_MAJOR_VERSION==y && \
   GMATELASTOPLASTIC_MINOR_VERSION==z)

// -------------------------------------------------------------------------------------------------

#endif
