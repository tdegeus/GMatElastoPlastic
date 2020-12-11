/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

*/

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_ARRAY_HPP
#define GMATELASTOPLASTIC_CARTESIAN3D_ARRAY_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlastic {
namespace Cartesian3d {

template <size_t N>
inline Array<N>::Array(const std::array<size_t, N>& shape)
{
    this->init(shape);
    m_type = xt::ones<size_t>(m_shape) * Type::Unset;
    m_index = xt::empty<size_t>(m_shape);
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::K() const
{
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            ret.data()[i] = 0.0;
            break;
        case Type::Elastic:
            ret.data()[i] = m_Elastic[m_index.data()[i]].K();
            break;
        case Type::LinearHardening:
            ret.data()[i] = m_LinearHardening[m_index.data()[i]].K();
            break;
        }
    }

    return ret;
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::G() const
{
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            ret.data()[i] = 0.0;
            break;
        case Type::Elastic:
            ret.data()[i] = m_Elastic[m_index.data()[i]].G();
            break;
        case Type::LinearHardening:
            ret.data()[i] = m_LinearHardening[m_index.data()[i]].G();
            break;
        }
    }

    return ret;
}

template <size_t N>
inline void Array<N>::epsp(xt::xtensor<double, N>& ret) const
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(ret, m_shape));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            ret.data()[i] = 0.0;
            break;
        case Type::Elastic:
            ret.data()[i] = 0.0;
            break;
        case Type::LinearHardening:
            ret.data()[i] = m_LinearHardening[m_index.data()[i]].epsp();
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::increment()
{
    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            break;
        case Type::Elastic:
            break;
        case Type::LinearHardening:
            m_LinearHardening[m_index.data()[i]].increment();
            break;
        }
    }
}

template <size_t N>
inline xt::xtensor<size_t, N> Array<N>::type() const
{
    return m_type;
}

template <size_t N>
inline xt::xtensor<size_t, N> Array<N>::isElastic() const
{
    xt::xtensor<size_t, N> ret = xt::where(xt::equal(m_type, Type::Elastic), 1ul, 0ul);
    return ret;
}

template <size_t N>
inline xt::xtensor<size_t, N> Array<N>::isPlastic() const
{
    xt::xtensor<size_t, N> ret = xt::where(xt::not_equal(m_type, Type::Elastic), 1ul, 0ul);
    return ret;
}

template <size_t N>
inline xt::xtensor<size_t, N> Array<N>::isLinearHardening() const
{
    xt::xtensor<size_t, N> ret = xt::where(xt::equal(m_type, Type::LinearHardening), 1ul, 0ul);
    return ret;
}

template <size_t N>
inline void Array<N>::setElastic(const xt::xtensor<size_t, N>& I, double K, double G)
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(m_type, I.shape()));
    GMATELASTOPLASTIC_ASSERT(xt::all(xt::equal(I, 0ul) || xt::equal(I, 1ul)));
    GMATELASTOPLASTIC_ASSERT(
        xt::all(xt::equal(xt::where(xt::equal(I, 1ul), m_type, Type::Unset), Type::Unset)));

    for (size_t i = 0; i < m_size; ++i) {
        if (I.data()[i] == 1ul) {
            m_type.data()[i] = Type::Elastic;
            m_index.data()[i] = m_Elastic.size();
            m_Elastic.push_back(Elastic(K, G));
        }
    }
}

template <size_t N>
inline void Array<N>::setLinearHardening(
    const xt::xtensor<size_t, N>& I,
    double K,
    double G,
    double sigy0,
    double H)
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(m_type, I.shape()));
    GMATELASTOPLASTIC_ASSERT(xt::all(xt::equal(I, 0ul) || xt::equal(I, 1ul)));
    GMATELASTOPLASTIC_ASSERT(
        xt::all(xt::equal(xt::where(xt::equal(I, 1ul), m_type, Type::Unset), Type::Unset)));

    for (size_t i = 0; i < m_size; ++i) {
        if (I.data()[i] == 1ul) {
            m_type.data()[i] = Type::LinearHardening;
            m_index.data()[i] = m_LinearHardening.size();
            m_LinearHardening.push_back(LinearHardening(K, G, sigy0, H));
        }
    }
}

template <size_t N>
inline void Array<N>::setElastic(
    const xt::xtensor<size_t, N>& I,
    const xt::xtensor<size_t, N>& idx,
    const xt::xtensor<double, 1>& K,
    const xt::xtensor<double, 1>& G)
{
    GMATELASTOPLASTIC_ASSERT(xt::amax(idx)[0] == K.size() - 1);
    GMATELASTOPLASTIC_ASSERT(K.size() == G.size());
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(m_type, I.shape()));
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(m_type, idx.shape()));
    GMATELASTOPLASTIC_ASSERT(xt::all(xt::equal(I, 0ul) || xt::equal(I, 1ul)));
    GMATELASTOPLASTIC_ASSERT(
        xt::all(xt::equal(xt::where(xt::equal(I, 1ul), m_type, Type::Unset), Type::Unset)));

    for (size_t i = 0; i < m_size; ++i) {
        if (I.data()[i] == 1ul) {
            size_t j = idx.data()[i];
            m_type.data()[i] = Type::Elastic;
            m_index.data()[i] = m_Elastic.size();
            m_Elastic.push_back(Elastic(K(j), G(j)));
        }
    }
}

template <size_t N>
inline void Array<N>::setLinearHardening(
    const xt::xtensor<size_t, N>& I,
    const xt::xtensor<size_t, N>& idx,
    const xt::xtensor<double, 1>& K,
    const xt::xtensor<double, 1>& G,
    const xt::xtensor<double, 1>& sigy0,
    const xt::xtensor<double, 1>& H)
{
    GMATELASTOPLASTIC_ASSERT(xt::amax(idx)[0] == K.size() - 1);
    GMATELASTOPLASTIC_ASSERT(K.size() == G.size());
    GMATELASTOPLASTIC_ASSERT(K.size() == sigy0.size());
    GMATELASTOPLASTIC_ASSERT(K.size() == H.size());
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(m_type, I.shape()));
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(m_type, idx.shape()));
    GMATELASTOPLASTIC_ASSERT(xt::all(xt::equal(I, 0ul) || xt::equal(I, 1ul)));
    GMATELASTOPLASTIC_ASSERT(
        xt::all(xt::equal(xt::where(xt::equal(I, 1ul), m_type, Type::Unset), Type::Unset)));

    for (size_t i = 0; i < m_size; ++i) {
        if (I.data()[i] == 1ul) {
            size_t j = idx.data()[i];
            m_type.data()[i] = Type::LinearHardening;
            m_index.data()[i] = m_LinearHardening.size();
            m_LinearHardening.push_back(LinearHardening(K(j), G(j), sigy0(j), H(j)));
        }
    }
}

template <size_t N>
inline void Array<N>::setStrain(const xt::xtensor<double, N + 2>& arg, bool tangent)
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(arg, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            break;
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].setStrainPtr(&arg.data()[i * m_stride_tensor2]);
            break;
        case Type::LinearHardening:
            m_LinearHardening[m_index.data()[i]].setStrainPtr(&arg.data()[i * m_stride_tensor2], tangent);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::strain(xt::xtensor<double, N + 2>& ret) const
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(ret, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            GMatTensor::Cartesian3d::pointer::O2(&ret.data()[i * m_stride_tensor2]);
            break;
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].strainPtr(&ret.data()[i * m_stride_tensor2]);
            break;
        case Type::LinearHardening:
            m_LinearHardening[m_index.data()[i]].strainPtr(&ret.data()[i * m_stride_tensor2]);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::stress(xt::xtensor<double, N + 2>& ret) const
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(ret, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            GMatTensor::Cartesian3d::pointer::O2(&ret.data()[i * m_stride_tensor2]);
            break;
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].stressPtr(&ret.data()[i * m_stride_tensor2]);
            break;
        case Type::LinearHardening:
            m_LinearHardening[m_index.data()[i]].stressPtr(&ret.data()[i * m_stride_tensor2]);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::tangent(xt::xtensor<double, N + 4>& ret) const
{
    GMATELASTOPLASTIC_ASSERT(xt::has_shape(ret, m_shape_tensor4));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            GMatTensor::Cartesian3d::pointer::O4(&ret.data()[i * m_stride_tensor4]);
            break;
        case Type::Elastic:
            m_Elastic[m_index.data()[i]].tangentPtr(&ret.data()[i * m_stride_tensor4]);
            break;
        case Type::LinearHardening:
            m_LinearHardening[m_index.data()[i]].tangentPtr(&ret.data()[i * m_stride_tensor4]);
            break;
        }
    }
}

template <size_t N>
inline xt::xtensor<double, N + 2> Array<N>::Strain() const
{
    xt::xtensor<double, N + 2> ret = xt::empty<double>(m_shape_tensor2);
    this->strain(ret);
    return ret;
}

template <size_t N>
inline xt::xtensor<double, N + 2> Array<N>::Stress() const
{
    xt::xtensor<double, N + 2> ret = xt::empty<double>(m_shape_tensor2);
    this->stress(ret);
    return ret;
}

template <size_t N>
inline xt::xtensor<double, N + 4> Array<N>::Tangent() const
{
    xt::xtensor<double, N + 4> ret = xt::empty<double>(m_shape_tensor4);
    this->tangent(ret);
    return ret;
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::Epsp() const
{
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);
    this->epsp(ret);
    return ret;
}

template <size_t N>
inline auto Array<N>::getElastic(const std::array<size_t, N>& index) const
{
    GMATELASTOPLASTIC_ASSERT(m_type[index] == Type::Elastic);
    return m_Elastic[m_index[index]];
}

template <size_t N>
inline auto Array<N>::getLinearHardening(const std::array<size_t, N>& index) const
{
    GMATELASTOPLASTIC_ASSERT(m_type[index] == Type::LinearHardening);
    return m_LinearHardening[m_index[index]];
}

template <size_t N>
inline auto* Array<N>::refElastic(const std::array<size_t, N>& index)
{
    GMATELASTOPLASTIC_ASSERT(m_type[index] == Type::Elastic);
    return &m_Elastic[m_index[index]];
}

template <size_t N>
inline auto* Array<N>::refLinearHardening(const std::array<size_t, N>& index)
{
    GMATELASTOPLASTIC_ASSERT(m_type[index] == Type::LinearHardening);
    return &m_LinearHardening[m_index[index]];
}

} // namespace Cartesian3d
} // namespace GMatElastoPlastic

#endif
