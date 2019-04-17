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
inline void LinearHardening::stress(const Tensor2& Eps, T&& Sig)
{
  // strain increment
  Tensor2 dEps = Eps - m_Eps_t;

  // trial elastic strain: presume the strain increment to be fully elastic
  xt::noalias(m_Epse) = m_Epse_t + dEps;

  // decompose trial elastic strain
  Tensor2 I     = Cartesian3d::I2();
  double  epsem = trace(m_Epse) / 3.0;
  Tensor2 Epsed = m_Epse - epsem * I;

  // trial stress
  double  sigm  = 3.0 * m_K * epsem;
  Tensor2 Sigd  = 2.0 * m_G * Epsed;
  double  sigeq = std::sqrt(1.5 * A2_ddot_B2(Sigd,Sigd));

  // evaluate the yield surface
  double phi = sigeq - (m_sigy0 + m_H * m_epsp_t);

  // return map
  if (phi > 0) {
    // - plastic flow
    double dgamma = phi / (3.0 * m_G + m_H);
    // - update trial stress (only the deviatoric part)
    Sigd *= (1.0 - 3.0 * m_G * dgamma / sigeq);
    // - update elastic strain (only the deviatoric part)
    xt::noalias(Epsed) = Sigd / (2.0 * m_G);
    xt::noalias(m_Epse) = epsem * I + Epsed;
    // - update equivalent plastic strain
    m_epsp = m_epsp_t + dgamma;
  }

  // compute stress
  xt::noalias(Sig) = sigm * I + Sigd;

  // store history
  xt::noalias(m_Eps) = Eps;
}

// -------------------------------------------------------------------------------------------------

inline Tensor2 LinearHardening::Stress(const Tensor2& Eps)
{
  Tensor2 Sig;
  this->stress(Eps, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

template <class T, class S>
inline void LinearHardening::tangent(const Tensor2& Eps, T&& Sig, S&& C)
{
  // unit tensors
  Tensor2 I   = Cartesian3d::I2();
  Tensor4 II  = Cartesian3d::II();
  Tensor4 I4d = Cartesian3d::I4d();

  // elastic tangent
  xt::noalias(C) = m_K * II + 2.0 * m_G * I4d;

  // strain increment
  Tensor2 dEps = Eps - m_Eps_t;

  // trial elastic strain: presume the strain increment to be fully elastic
  xt::noalias(m_Epse) = m_Epse_t + dEps;

  // decompose trial elastic strain
  double  epsem = trace(m_Epse) / 3.0;
  Tensor2 Epsed = m_Epse - epsem * I;

  // trial stress
  double  sigm  = 3.0 * m_K * epsem;
  Tensor2 Sigd  = 2.0 * m_G * Epsed;
  double  sigeq = std::sqrt(1.5 * A2_ddot_B2(Sigd,Sigd));

  // evaluate the yield surface
  double phi = sigeq - (m_sigy0 + m_H * m_epsp_t);

  // return map
  if (phi > 0) {
    // - plastic flow
    double dgamma = phi / (3.0 * m_G + m_H);
    // - direction of plastic flow
    Tensor2 N = 1.5 * Sigd / sigeq;
    // - update trial stress (only the deviatoric part)
    Sigd *= (1.0 - 3.0 * m_G * dgamma / sigeq);
    // - update elastic strain (only the deviatoric part)
    xt::noalias(Epsed) = Sigd / (2.0 * m_G);
    xt::noalias(m_Epse) = epsem * I + Epsed;
    // - update equivalent plastic strain
    m_epsp = m_epsp_t + dgamma;
    // - update tangent
    Tensor4 NN;
    A2_dyadic_B2(N, N, NN);
    Tensor4 K = C; // use temporary to overcome xtensor bug: remove when fixed in xtensor
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

inline std::tuple<Tensor2,Tensor4> LinearHardening::Tangent(const Tensor2& Eps)
{
  Tensor2 Sig;
  Tensor4 C;
  this->tangent(Eps, Sig, C);
  return std::make_tuple(Sig, C);
}

// -------------------------------------------------------------------------------------------------

inline void LinearHardening::increment()
{
  m_epsp_t = m_epsp;
  xt::noalias(m_Epse_t) = m_Epse;
  xt::noalias(m_Eps_t ) = m_Eps;
}

// -------------------------------------------------------------------------------------------------

inline double LinearHardening::epsp() const
{
  return m_epsp;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
