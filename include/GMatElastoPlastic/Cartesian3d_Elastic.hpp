/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

================================================================================================= */

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_ELASTIC_HPP
#define GMATELASTOPLASTIC_CARTESIAN3D_ELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Elastic::Elastic(double K, double G) : m_K(K), m_G(G)
{
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::K() const
{
  return m_K;
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::G() const
{
  return m_G;
}

// -------------------------------------------------------------------------------------------------

template <class T>
inline void Elastic::stress(const Tensor2& Eps, T&& Sig) const
{
  auto I    = Cartesian3d::I2();
  auto epsm = trace(Eps) / 3.0;
  auto Epsd = Eps - epsm * I;
  xt::noalias(Sig) = 3.0 * m_K * epsm * I + 2.0 * m_G * Epsd;
}

// -------------------------------------------------------------------------------------------------

inline Tensor2 Elastic::Stress(const Tensor2& Eps) const
{
  Tensor2 Sig;
  this->stress(Eps, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

template <class T, class S>
inline void Elastic::tangent(const Tensor2& Eps, T&& Sig, S&& C) const
{
  auto I    = Cartesian3d::I2();
  auto II   = Cartesian3d::II();
  auto I4d  = Cartesian3d::I4d();
  auto epsm = trace(Eps) / 3.0;
  auto Epsd = Eps - epsm * I;
  xt::noalias(Sig) = 3.0 * m_K * epsm * I + 2.0 * m_G * Epsd;
  xt::noalias(C) = m_K * II + 2.0 * m_G * I4d;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<Tensor2,Tensor4> Elastic::Tangent(const Tensor2& Eps) const
{
  Tensor2 Sig;
  Tensor4 C;
  this->tangent(Eps, Sig, C);
  return std::make_tuple(Sig, C);
}

// -------------------------------------------------------------------------------------------------

inline void Elastic::increment() const
{
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::epsp() const
{
  return 0.0;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
