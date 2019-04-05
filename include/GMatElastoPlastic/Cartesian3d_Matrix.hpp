/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

================================================================================================= */

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_MATRIX_HPP
#define GMATELASTOPLASTIC_CARTESIAN3D_MATRIX_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Matrix::Matrix(size_t nelem, size_t nip) : m_nelem(nelem), m_nip(nip)
{
  m_type   = xt::ones <size_t>({m_nelem, m_nip}) * Type::Unset;
  m_index  = xt::empty<size_t>({m_nelem, m_nip});
  m_allSet = false;
}

// -------------------------------------------------------------------------------------------------

inline size_t Matrix::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

inline size_t Matrix::nelem() const
{
  return m_nelem;
}

// -------------------------------------------------------------------------------------------------

inline size_t Matrix::nip() const
{
  return m_nip;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> Matrix::type() const
{
  return m_type;
};

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Matrix::K() const
{
  GMATELASTOPLASTIC_ASSERT(m_allSet);

  xt::xtensor<double,2> out = xt::empty<double>({m_nelem, m_nip});
  #pragma omp parallel for
  for (size_t e = 0; e < m_nelem; ++e) {
    for (size_t q = 0; q < m_nip; ++q) {
      switch (m_type(e,q)) {
        case Type::Elastic:         out(e,q) = m_Elastic        [m_index(e,q)].K(); break;
        case Type::LinearHardening: out(e,q) = m_LinearHardening[m_index(e,q)].K(); break;
      }
    }
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Matrix::G() const
{
  GMATELASTOPLASTIC_ASSERT(m_allSet);

  xt::xtensor<double,2> out = xt::empty<double>({m_nelem, m_nip});
  #pragma omp parallel for
  for (size_t e = 0; e < m_nelem; ++e) {
    for (size_t q = 0; q < m_nip; ++q) {
      switch (m_type(e,q)) {
        case Type::Elastic:         out(e,q) = m_Elastic        [m_index(e,q)].G(); break;
        case Type::LinearHardening: out(e,q) = m_LinearHardening[m_index(e,q)].G(); break;
      }
    }
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Matrix::I() const
{
  xt::xtensor<double,4> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});
  #pragma omp parallel
  {
    T2 unit = Cartesian3d::I();
    #pragma omp for
    for (size_t e = 0; e < m_nelem; ++e) {
      for (size_t q = 0; q < m_nip; ++q) {
        auto view = xt::adapt(&out(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        xt::noalias(view) = unit;
      }
    }
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::II() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});
  #pragma omp parallel
  {
    T4 unit = Cartesian3d::II();
    #pragma omp for
    for (size_t e = 0; e < m_nelem; ++e) {
      for (size_t q = 0; q < m_nip; ++q) {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        xt::noalias(view) = unit;
      }
    }
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});
  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4();
    #pragma omp for
    for (size_t e = 0; e < m_nelem; ++e) {
      for (size_t q = 0; q < m_nip; ++q) {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        xt::noalias(view) = unit;
      }
    }
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4rt() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});
  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4rt();
    #pragma omp for
    for (size_t e = 0; e < m_nelem; ++e) {
      for (size_t q = 0; q < m_nip; ++q) {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        xt::noalias(view) = unit;
      }
    }
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4s() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});
  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4s();
    #pragma omp for
    for (size_t e = 0; e < m_nelem; ++e) {
      for (size_t q = 0; q < m_nip; ++q) {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        xt::noalias(view) = unit;
      }
    }
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4d() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});
  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4d();
    #pragma omp for
    for (size_t e = 0; e < m_nelem; ++e) {
      for (size_t q = 0; q < m_nip; ++q) {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        xt::noalias(view) = unit;
      }
    }
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::check() const
{
  if (xt::any(xt::equal(m_type, Type::Unset)))
    throw std::runtime_error("Points without material found");
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::checkAllSet()
{
  if (xt::any(xt::equal(m_type, Type::Unset)))
    m_allSet = false;
  else
    m_allSet = true;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::setElastic(
  const xt::xtensor<size_t,2>& I,
  double K,
  double G)
{
  GMATELASTOPLASTIC_ASSERT(m_type.shape() == I.shape());
  GMATELASTOPLASTIC_ASSERT(xt::all(xt::equal(I,0ul) || xt::equal(I,1ul)));
  GMATELASTOPLASTIC_ASSERT(\
    xt::all(xt::equal(xt::where(xt::equal(I,1ul), m_type, Type::Unset), Type::Unset)));

  m_type = xt::where(xt::equal(I, 1ul), Type::Elastic, m_type);
  m_index = xt::where(xt::equal(I, 1ul), m_Elastic.size(), m_index);
  this->checkAllSet();
  m_Elastic.push_back(Elastic(K, G));
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::setLinearHardening(
  const xt::xtensor<size_t,2>& I,
  double K,
  double G,
  double sigy0,
  double H)
{
  GMATELASTOPLASTIC_ASSERT(m_type.shape() == I.shape());
  GMATELASTOPLASTIC_ASSERT(xt::all(xt::equal(I,0ul) || xt::equal(I,1ul)));
  GMATELASTOPLASTIC_ASSERT(\
    xt::all(xt::equal(xt::where(xt::equal(I,1ul), m_type, Type::Unset), Type::Unset)));

  for (size_t e = 0; e < m_nelem; ++e) {
    for (size_t q = 0; q < m_nip; ++q) {
      if (I(e,q) == 1ul) {
        m_type(e,q) = Type::LinearHardening;
        m_index(e,q) = m_LinearHardening.size();
        m_LinearHardening.push_back(LinearHardening(K, G, sigy0, H));
      }
    }
  }
  this->checkAllSet();
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::setElastic(
  const xt::xtensor<size_t,2>& I,
  const xt::xtensor<size_t,2>& idx,
  const xt::xtensor<double,1>& K,
  const xt::xtensor<double,1>& G)
{
  GMATELASTOPLASTIC_ASSERT(xt::amax(idx)[0] == K.size()-1);
  GMATELASTOPLASTIC_ASSERT(K.size() == G.size());
  GMATELASTOPLASTIC_ASSERT(m_type.shape() == idx.shape());
  GMATELASTOPLASTIC_ASSERT(m_type.shape() == I.shape());
  GMATELASTOPLASTIC_ASSERT(xt::all(xt::equal(I,0ul) || xt::equal(I,1ul)));
  GMATELASTOPLASTIC_ASSERT(\
    xt::all(xt::equal(xt::where(xt::equal(I,1ul), m_type, Type::Unset), Type::Unset)));

  m_type = xt::where(xt::equal(I, 1ul), Type::Elastic, m_type);
  m_index = xt::where(xt::equal(I, 1ul), m_Elastic.size() + idx, m_index);
  this->checkAllSet();
  for (size_t i = 0; i < K.size(); ++i)
    m_Elastic.push_back(Elastic(K(i), G(i)));
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::setLinearHardening(
  const xt::xtensor<size_t,2>& I,
  const xt::xtensor<size_t,2>& idx,
  const xt::xtensor<double,1>& K,
  const xt::xtensor<double,1>& G,
  const xt::xtensor<double,1>& sigy0,
  const xt::xtensor<double,1>& H)
{
  GMATELASTOPLASTIC_ASSERT(xt::amax(idx)[0] == K.size()-1);
  GMATELASTOPLASTIC_ASSERT(K.size() == G.size());
  GMATELASTOPLASTIC_ASSERT(K.size() == H.size());
  GMATELASTOPLASTIC_ASSERT(m_type.shape() == idx.shape());
  GMATELASTOPLASTIC_ASSERT(m_type.shape() == I.shape());
  GMATELASTOPLASTIC_ASSERT(xt::all(xt::equal(I,0ul) || xt::equal(I,1ul)));
  GMATELASTOPLASTIC_ASSERT(\
    xt::all(xt::equal(xt::where(xt::equal(I,1ul), m_type, Type::Unset), Type::Unset)));

  for (size_t e = 0; e < m_nelem; ++e) {
    for (size_t q = 0; q < m_nip; ++q) {
      if (I(e,q) == 1ul) {
        m_type(e,q) = Type::LinearHardening;
        m_index(e,q) = m_LinearHardening.size();
        m_LinearHardening.push_back(
          LinearHardening(K(idx(e,q)), G(idx(e,q)), sigy0(idx(e,q)), H(idx(e,q))));
      }
    }
  }
  this->checkAllSet();
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::stress(const xt::xtensor<double,4>& a_Eps, xt::xtensor<double,4>& a_Sig)
{
  GMATELASTOPLASTIC_ASSERT(m_allSet);
  GMATELASTOPLASTIC_ASSERT(a_Eps.shape() == \
    std::decay_t<decltype(a_Eps)>::shape_type({m_nelem, m_nip, m_ndim, m_ndim}));
  GMATELASTOPLASTIC_ASSERT(a_Eps.shape() == a_Sig.shape());

  #pragma omp parallel for
  for (size_t e = 0; e < m_nelem; ++e) {
    for (size_t q = 0; q < m_nip; ++q) {
      auto Eps = xt::adapt(&a_Eps(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
      auto Sig = xt::adapt(&a_Sig(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
      switch (m_type(e,q)) {
        case Type::Elastic:         m_Elastic        [m_index(e,q)].stress(Eps, Sig); break;
        case Type::LinearHardening: m_LinearHardening[m_index(e,q)].stress(Eps, Sig); break;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::tangent(
  const xt::xtensor<double,4>& a_Eps,
        xt::xtensor<double,4>& a_Sig,
        xt::xtensor<double,6>& a_Tangent)
{
  GMATELASTOPLASTIC_ASSERT(m_allSet);
  GMATELASTOPLASTIC_ASSERT(a_Eps.shape() == \
    std::decay_t<decltype(a_Eps)>::shape_type({m_nelem, m_nip, m_ndim, m_ndim}));
  GMATELASTOPLASTIC_ASSERT(a_Eps.shape() == a_Sig.shape());
  GMATELASTOPLASTIC_ASSERT(a_Tangent.shape() == \
    std::decay_t<decltype(a_Tangent)>::shape_type({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim}));

  #pragma omp parallel for
  for (size_t e = 0; e < m_nelem; ++e) {
    for (size_t q = 0; q < m_nip; ++q) {
      auto Eps = xt::adapt(&a_Eps(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
      auto Sig = xt::adapt(&a_Sig(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
      auto C = xt::adapt(&a_Tangent(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
      switch (m_type(e,q)) {
        case Type::Elastic:         m_Elastic        [m_index(e,q)].tangent(Eps, Sig, C); break;
        case Type::LinearHardening: m_LinearHardening[m_index(e,q)].tangent(Eps, Sig, C); break;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::epsp(xt::xtensor<double,2>& epsp) const
{
  GMATELASTOPLASTIC_ASSERT(m_allSet);
  GMATELASTOPLASTIC_ASSERT(epsp.shape() == \
    std::decay_t<decltype(epsp)>::shape_type({m_nelem, m_nip}));

  #pragma omp parallel for
  for (size_t e = 0; e < m_nelem; ++e) {
    for (size_t q = 0; q < m_nip; ++q) {
      switch (m_type(e,q)) {
        case Type::Elastic:         epsp(e,q) = m_Elastic        [m_index(e,q)].epsp(); break;
        case Type::LinearHardening: epsp(e,q) = m_LinearHardening[m_index(e,q)].epsp(); break;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::increment()
{
  GMATELASTOPLASTIC_ASSERT(m_allSet);

  #pragma omp parallel for
  for (size_t e = 0; e < m_nelem; ++e) {
    for (size_t q = 0; q < m_nip; ++q) {
      switch (m_type(e,q)) {
        case Type::Elastic:         m_Elastic        [m_index(e,q)].increment(); break;
        case Type::LinearHardening: m_LinearHardening[m_index(e,q)].increment(); break;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Matrix::Stress(const xt::xtensor<double,4>& Eps)
{
  xt::xtensor<double,4> Sig = xt::empty<double>(Eps.shape());
  this->stress(Eps, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<xt::xtensor<double,4>,xt::xtensor<double,6>> Matrix::Tangent(
  const xt::xtensor<double,4>& Eps)
{
  xt::xtensor<double,4> Sig = xt::empty<double>({m_nelem,m_nip,m_ndim,m_ndim});
  xt::xtensor<double,6> C = xt::empty<double>({m_nelem,m_nip,m_ndim,m_ndim,m_ndim,m_ndim});
  this->tangent(Eps, Sig, C);
  return std::make_tuple(Sig, C);
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
