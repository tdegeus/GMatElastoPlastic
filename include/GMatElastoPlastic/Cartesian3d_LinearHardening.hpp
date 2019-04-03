/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

================================================================================================= */

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_LINEARHARDENING_HPP
#define GMATELASTOPLASTIC_CARTESIAN3D_LINEARHARDENING_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline LinearHardening::LinearHardening(double K, double G, double sigy0, double H) :
  m_K(K), m_G(G), m_sigy0(sigy0), m_H(H)
{
}

// -------------------------------------------------------------------------------------------------

inline double LinearHardening::K() const
{
  return m_K;
}

// -------------------------------------------------------------------------------------------------

inline double LinearHardening::G() const
{
  return m_G;
}

// -------------------------------------------------------------------------------------------------

inline double LinearHardening::sigy0() const
{
  return m_sigy0;
}

// -------------------------------------------------------------------------------------------------

inline double LinearHardening::H() const
{
  return m_H;
}

// -------------------------------------------------------------------------------------------------

template <class T>
inline void LinearHardening::stress(const T2& Eps, T&& Sig)
{
  // strain increment
  T2 dEps = Eps - m_Eps_t;

  // trial elastic strain: presume the strain increment to be fully elastic
  xt::noalias(m_Epse) = m_Epse_t + dEps;

  // decompose trial elastic strain
  T2     I     = Cartesian3d::I();
  double epsem = trace(m_Epse) / 3.0;
  T2     Epsed = m_Epse - epsem * I;

  // trial stress
  double sigm  = 3.0 * m_K * epsem;
  T2     Sigd  = 2.0 * m_G * Epsed;
  double sigeq = std::sqrt(1.5 * ddot22(Sigd,Sigd));

  // evaluate the yield surface
  double phi = sigeq - (m_sigy0 + m_H * m_epsp_t);

  // return map
  if ( phi > 0 ) {
    // - plastic flow
    double dgamma = phi / (3.0 * m_G + m_H);
    // - update stress, elastic strain, and equivalent plastic strain
    xt::noalias(Sigd) = (1.0 - 3.0 * m_G * dgamma / sigeq) * Sigd;
    xt::noalias(Epsed) = Sigd / (2.0 * m_G);
    xt::noalias(m_Epse) = epsem * I + Epsed;
    m_epsp = m_epsp_t + dgamma;
  }

  // compute stress
  xt::noalias(Sig) = sigm * I + Sigd;

  // store history
  xt::noalias(m_Eps) = Eps;
}

// -------------------------------------------------------------------------------------------------

inline T2 LinearHardening::Stress(const T2& Eps)
{
  T2 Sig;
  this->stress(Eps, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

template <class T, class S>
inline void LinearHardening::tangent(const T2& Eps, T&& Sig, S&& C)
{
  // unit tensors
  T2 I   = Cartesian3d::I();
  T4 II  = Cartesian3d::II();
  T4 I4d = Cartesian3d::I4d();

  // elastic tangent
  xt::noalias(C) = m_K * II + 2.0 * m_G * I4d;

  // strain increment
  T2 dEps = Eps - m_Eps_t;

  // trial elastic strain: presume the strain increment to be fully elastic
  xt::noalias(m_Epse) = m_Epse_t + dEps;

  // decompose trial elastic strain
  double epsem = trace(m_Epse) / 3.0;
  T2     Epsed = m_Epse - epsem * I;

  // trial stress
  double sigm  = 3.0 * m_K * epsem;
  T2     Sigd  = 2.0 * m_G * Epsed;
  double sigeq = std::sqrt(1.5 * ddot22(Sigd,Sigd));

  // evaluate the yield surface
  double phi = sigeq - (m_sigy0 + m_H * m_epsp_t);

  // return map
  if ( phi > 0 ) {
    // - plastic flow (direction)
    double dgamma = phi / (3.0 * m_G + m_H);
    T2 N = 1.5 * Sigd / sigeq;
    // - update stress, elastic strain, and equivalent plastic strain
    xt::noalias(Sigd) = (1.0 - 3.0 * m_G * dgamma / sigeq) * Sigd;
    xt::noalias(Epsed) = Sigd / (2.0 * m_G);
    xt::noalias(m_Epse) = epsem * I + Epsed;
    m_epsp = m_epsp_t + dgamma;
    // - update tangent
    T4 NN;
    dyadic22(N, N, NN);
    T4 K = C; // use temporary to overcome xtensor bug: remove when fixed in xtensor
    K -= 6.0 * std::pow(m_G,2.0) * dgamma / sigeq * I4d;
    K += 4.0 * std::pow(m_G,2.0) * (dgamma / sigeq - 1.0 / (3.0 * m_G + m_H)) * NN;
    xt::noalias(C) = K;
  }

  // compute stress
  xt::noalias(Sig) = sigm * I + Sigd;

  // store history
  xt::noalias(m_Eps) = Eps;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<T2,T4> LinearHardening::Tangent(const T2& Eps)
{
  T2 Sig;
  T4 C;
  this->tangent(Eps, Sig, C);
  return std::make_tuple(Sig, C);
}

// -------------------------------------------------------------------------------------------------

inline void LinearHardening::increment()
{
  m_epsp_t = m_epsp;
  xt::noalias(m_Epse_t) = m_Epse;
  xt::noalias(m_Eps_t) = m_Eps;
}

// -------------------------------------------------------------------------------------------------

inline double LinearHardening::epsp() const
{
  return m_epsp;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
