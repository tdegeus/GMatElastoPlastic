/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

*/

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_ELASTIC_HPP
#define GMATELASTOPLASTIC_CARTESIAN3D_ELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlastic {
namespace Cartesian3d {

inline Elastic::Elastic(double K, double G) :
    GMatElastic::Cartesian3d::Elastic(K, G)
{
}

} // namespace Cartesian3d
} // namespace GMatElastoPlastic

#endif
