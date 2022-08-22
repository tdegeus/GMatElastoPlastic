/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_H
#define GMATELASTOPLASTIC_CARTESIAN3D_H

#include <GMatElastic/Cartesian3d.h>
#include <GMatTensor/Cartesian3d.h>

#include "config.h"
#include "version.h"

namespace GMatElastoPlastic {

/**
Implementation in a 3-d Cartesian coordinate frame.
*/
namespace Cartesian3d {

using GMatElastic::Cartesian3d::epseq;
using GMatElastic::Cartesian3d::Epseq;
using GMatElastic::Cartesian3d::sigeq;
using GMatElastic::Cartesian3d::Sigeq;

/**
Array of material points with a linear elasto-plastic constitutive response with linear hardening.

\tparam N Rank of the array.
*/
template <size_t N>
class LinearHardening : public GMatElastic::Cartesian3d::Elastic<N> {
private:
    array_type::tensor<double, N> m_sigy0; ///< Initial yield stress per item.
    array_type::tensor<double, N> m_H; ///< Hardening modulus per item.
    array_type::tensor<double, N> m_epsp; ///< Equivalent plastic strain tensor per item.
    array_type::tensor<double, N> m_epsp_t; ///< `epsp` as the last increment.
    array_type::tensor<double, N + 2> m_Eps_t; ///< `Eps` as the last increment.
    array_type::tensor<double, N + 2> m_Epse; ///< Elastic strain tensor per item.
    array_type::tensor<double, N + 2> m_Epse_t; ///< `Epse` as the last increment.

    using GMatElastic::Cartesian3d::Elastic<N>::m_K;
    using GMatElastic::Cartesian3d::Elastic<N>::m_G;
    using GMatElastic::Cartesian3d::Elastic<N>::m_Eps;
    using GMatElastic::Cartesian3d::Elastic<N>::m_Sig;
    using GMatElastic::Cartesian3d::Elastic<N>::m_C;

    using GMatTensor::Cartesian3d::Array<N>::m_ndim;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor4;
    using GMatTensor::Cartesian3d::Array<N>::m_size;
    using GMatTensor::Cartesian3d::Array<N>::m_shape;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor4;

    using GMatElastic::Cartesian3d::Elastic<N>::energy;

public:
    using GMatTensor::Cartesian3d::Array<N>::rank;

    LinearHardening() = default;

    /**
    Construct system.
    \param K Bulk modulus per item.
    \param G Shear modulus per item.
    \param sigy0 Initial yield stress per item.
    \param H Hardening modulus per item.
    */
    template <class T>
    LinearHardening(const T& K, const T& G, const T& sigy0, const T& H)
    {
        GMATELASTIC_ASSERT(xt::has_shape(K, sigy0.shape()));
        GMATELASTIC_ASSERT(xt::has_shape(K, H.shape()));

        this->init_Elastic(K, G);

        m_sigy0 = sigy0;
        m_H = H;

        m_Eps = xt::zeros<double>(m_shape_tensor2);
        m_Eps_t = xt::zeros<double>(m_shape_tensor2);
        m_Epse = xt::zeros<double>(m_shape_tensor2);
        m_Epse_t = xt::zeros<double>(m_shape_tensor2);
        m_epsp = xt::zeros<double>(m_shape);
        m_epsp_t = xt::zeros<double>(m_shape);

        this->refresh(true); // initialise tangent (elastic)
    }

    /**
    Initial yield stress per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& sigy0() const
    {
        return m_sigy0;
    }

    /**
    Hardening modulus per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& H() const
    {
        return m_H;
    }

    void refresh(bool compute_tangent = true) override
    {
        namespace GT = GMatTensor::Cartesian3d::pointer;

#pragma omp parallel
        {
            auto I = GMatTensor::Cartesian3d::I2();
            auto II = GMatTensor::Cartesian3d::II();
            auto I4d = GMatTensor::Cartesian3d::I4d();

            double K;
            double G;
            double sigy0;
            double H;
            double epsp_t;

            auto Eps = xt::adapt(m_Eps.data(), {m_ndim, m_ndim});
            auto Eps_t = xt::adapt(m_Eps_t.data(), {m_ndim, m_ndim});
            auto Epse = xt::adapt(m_Epse.data(), {m_ndim, m_ndim});
            auto Epse_t = xt::adapt(m_Epse_t.data(), {m_ndim, m_ndim});
            auto Sig = xt::adapt(m_Sig.data(), {m_ndim, m_ndim});
            auto C = xt::adapt(m_C.data(), {m_ndim, m_ndim, m_ndim, m_ndim});

            auto Epsed = xt::empty_like(I);

#pragma omp for
            for (size_t i = 0; i < m_size; ++i) {

                K = m_K.flat(i);
                G = m_G.flat(i);
                sigy0 = m_sigy0.flat(i);
                H = m_H.flat(i);
                epsp_t = m_epsp_t.flat(i);

                Eps.reset_buffer(&m_Eps.flat(i * m_stride_tensor2), m_stride_tensor2);
                Eps_t.reset_buffer(&m_Eps_t.flat(i * m_stride_tensor2), m_stride_tensor2);
                Epse.reset_buffer(&m_Epse.flat(i * m_stride_tensor2), m_stride_tensor2);
                Epse_t.reset_buffer(&m_Epse_t.flat(i * m_stride_tensor2), m_stride_tensor2);
                Sig.reset_buffer(&m_Sig.flat(i * m_stride_tensor2), m_stride_tensor2);
                C.reset_buffer(&m_C.flat(i * m_stride_tensor4), m_stride_tensor4);

                // trial elastic strain: presume the strain increment to be fully elastic
                xt::noalias(Epse) = Epse_t + (Eps - Eps_t);

                // decompose trial elastic strain
                double epsem = GT::Hydrostatic_deviatoric(Epse.data(), Epsed.data());

                // trial stress
                double sigm = 3.0 * K * epsem;
                auto Sigd = xt::eval(2.0 * G * Epsed);
                double sigeq = std::sqrt(1.5 * GT::A2s_ddot_B2s(Sigd.data(), Sigd.data()));

                // evaluate the yield surface
                double phi = sigeq - (sigy0 + H * epsp_t);

                // elastic tangent
                if (compute_tangent) {
                    xt::noalias(C) = K * II + 2.0 * G * I4d;
                }

                // return map
                if (phi > 0) {
                    auto N2 = xt::eval(1.5 * Sigd / sigeq);
                    double dgamma = phi / (3.0 * G + H);
                    Sigd *= (1.0 - 3.0 * G * dgamma / sigeq);
                    xt::noalias(Epsed) = Sigd / (2.0 * G);
                    xt::noalias(Epse) = epsem * I + Epsed;
                    m_epsp.flat(i) = epsp_t + dgamma;
                    if (compute_tangent) {
                        auto NN = xt::empty_like(II);
                        GT::A2_dyadic_B2(N2.data(), N2.data(), NN.data());
                        C -= 6.0 * std::pow(G, 2.0) * dgamma / sigeq * I4d;
                        C += 4.0 * std::pow(G, 2.0) * (dgamma / sigeq - 1.0 / (3.0 * G + H)) * NN;
                    }
                }

                xt::noalias(Sig) = sigm * I + Sigd;
            }
        }
    }

    /**
    Plastic strain per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& epsp() const
    {
        return m_epsp;
    }

    /**
    Update history variables.
    */
    void increment()
    {
        std::copy(m_epsp.cbegin(), m_epsp.cend(), m_epsp_t.begin());
        std::copy(m_Epse.cbegin(), m_Epse.cend(), m_Epse_t.begin());
        std::copy(m_Eps.cbegin(), m_Eps.cend(), m_Eps_t.begin());
    }
};

} // namespace Cartesian3d
} // namespace GMatElastoPlastic

#endif
