/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

*/

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_LINEARHARDENING_HPP
#define GMATELASTOPLASTIC_CARTESIAN3D_LINEARHARDENING_HPP

#include "Cartesian3d.h"
#include <xtensor/xio.hpp>

namespace GMatElastoPlastic {
namespace Cartesian3d {

inline LinearHardening::LinearHardening(double K, double G, double sigy0, double H)
    : m_K(K), m_G(G), m_sigy0(sigy0), m_H(H)
{
    m_epsp = 0.0;
    m_epsp_t = 0.0;
    m_Eps = xt::zeros<double>({3, 3});
    m_Eps_t = xt::zeros<double>({3, 3});
    m_Epse = xt::zeros<double>({3, 3});
    m_Epse_t = xt::zeros<double>({3, 3});
    m_Sig = xt::empty<double>({3, 3});
    m_C = xt::empty<double>({3, 3, 3, 3});
}

inline double LinearHardening::K() const
{
    return m_K;
}

inline double LinearHardening::G() const
{
    return m_G;
}

inline double LinearHardening::sigy0() const
{
    return m_sigy0;
}

inline double LinearHardening::H() const
{
    return m_H;
}

inline double LinearHardening::epsp() const
{
    return m_epsp;
}

inline void LinearHardening::increment()
{
    m_epsp_t = m_epsp;
    xt::noalias(m_Epse_t) = m_Epse;
    xt::noalias(m_Eps_t) = m_Eps;
}

template <class T>
inline void LinearHardening::setStrainPtr(const T* arg, bool tangent)
{
    namespace GT = GMatTensor::Cartesian3d::pointer;
    std::copy(arg, arg + 9, m_Eps.data());

    auto I = Cartesian3d::I2();
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();

    // trial elastic strain: presume the strain increment to be fully elastic
    xt::noalias(m_Epse) = xt::eval(m_Epse_t + (m_Eps - m_Eps_t));

    // decompose trial elastic strain
    auto Epsed = decltype(m_Epse)::from_shape(m_Epse.shape());
    double epsem = GT::Hydrostatic_deviatoric(m_Epse.data(), Epsed.data());

    // trial stress
    double sigm = 3.0 * m_K * epsem;
    auto Sigd = xt::eval(2.0 * m_G * Epsed);
    double sigeq = std::sqrt(1.5 * GT::A2s_ddot_B2s(Sigd.data(), Sigd.data()));

    // evaluate the yield surface
    double phi = sigeq - (m_sigy0 + m_H * m_epsp_t);

    // elastic tangent
    if (tangent) {
        xt::noalias(m_C) = m_K * II + 2.0 * m_G * I4d;
    }

    // return map
    if (phi > 0) {
        auto N = xt::eval(1.5 * Sigd / sigeq);
        double dgamma = phi / (3.0 * m_G + m_H);
        Sigd *= (1.0 - 3.0 * m_G * dgamma / sigeq);
        xt::noalias(Epsed) = Sigd / (2.0 * m_G);
        xt::noalias(m_Epse) = epsem * I + Epsed;
        m_epsp = m_epsp_t + dgamma;
        if (tangent) {
            auto NN = decltype(m_C)::from_shape(m_C.shape());
            GT::A2_dyadic_B2(N.data(), N.data(), NN.data());
            m_C -= 6.0 * std::pow(m_G, 2.0) * dgamma / sigeq * I4d;
            m_C += 4.0 * std::pow(m_G, 2.0) * (dgamma / sigeq - 1.0 / (3.0 * m_G + m_H)) * NN;
        }
    }

    xt::noalias(m_Sig) = sigm * I + Sigd;
}

template <class T>
inline void LinearHardening::strainPtr(T* ret) const
{
    std::copy(m_Eps.cbegin(), m_Eps.cend(), ret);
}

template <class T>
inline void LinearHardening::stressPtr(T* ret) const
{
    std::copy(m_Sig.cbegin(), m_Sig.cend(), ret);
}

template <class T>
inline void LinearHardening::tangentPtr(T* ret) const
{
    std::copy(m_C.cbegin(), m_C.cend(), ret);
}

template <class T>
inline void LinearHardening::setStrain(const T& arg, bool tangent)
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(arg, {3, 3}));
    return this->setStrainPtr(arg.data(), tangent);
}

template <class T>
inline void LinearHardening::strain(T& ret) const
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(ret, {3, 3}));
    return this->strainPtr(ret.data());
}

template <class T>
inline void LinearHardening::stress(T& ret) const
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(ret, {3, 3}));
    return this->stressPtr(ret.data());
}

template <class T>
inline void LinearHardening::tangent(T& ret) const
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(ret, {3, 3, 3, 3}));
    return this->tangentPtr(ret.data());
}

inline xt::xtensor<double, 2> LinearHardening::Strain() const
{
    return m_Eps;
}

inline xt::xtensor<double, 2> LinearHardening::Stress() const
{
    return m_Sig;
}

inline xt::xtensor<double, 4> LinearHardening::Tangent() const
{
    return m_C;
}

} // namespace Cartesian3d
} // namespace GMatElastoPlastic

#endif
