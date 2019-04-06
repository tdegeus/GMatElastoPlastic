/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

================================================================================================= */

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_HPP
#define GMATELASTOPLASTIC_CARTESIAN3D_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Tensor2 I()
{
  return Tensor2({{1.0, 0.0, 0.0},
                  {0.0, 1.0, 0.0},
                  {0.0, 0.0, 1.0}});
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 II()
{
  Tensor4 out;
  out.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          if (i == j and k == l)
            out(i,j,k,l) = 1.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 I4()
{
  Tensor4 out;
  out.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          if (i == l and j == k)
            out(i,j,k,l) = 1.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 I4rt()
{
  Tensor4 out;
  out.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          if (i == k and j == l)
            out(i,j,k,l) = 1.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 I4s()
{
  return 0.5 * ( I4() + I4rt() );
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 I4d()
{
  return I4s() - II() / 3.0;
}

// -------------------------------------------------------------------------------------------------

inline double Hydrostatic(const Tensor2& A)
{
  return trace(A) / 3.0;
}

// -------------------------------------------------------------------------------------------------

inline Tensor2 Deviatoric(const Tensor2& A)
{
  return A - trace(A) / 3.0 * I();
}

// -------------------------------------------------------------------------------------------------

inline double Epseq(const Tensor2& Eps)
{
  Tensor2 Epsd = Eps - trace(Eps) / 3.0 * I();
  return std::sqrt(2.0/3.0 * A2_ddot_B2(Epsd,Epsd));
}

// -------------------------------------------------------------------------------------------------

inline double Sigeq(const Tensor2& Sig)
{
  Tensor2 Sigd = Sig - trace(Sig) / 3.0 * I();
  return std::sqrt(1.5 * A2_ddot_B2(Sigd,Sigd));
}

// -------------------------------------------------------------------------------------------------

inline void hydrostatic(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Am)
{
  GMATELASTOPLASTIC_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Am.shape()[0], Am.shape()[1], 3, 3}));

  #pragma omp parallel
  {
    #pragma omp for
    for (size_t e = 0; e < A.shape()[0]; ++e) {
      for (size_t q = 0; q < A.shape()[1]; ++q) {
        auto Ai = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        Am(e,q) = trace(Ai) / 3.0;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void deviatoric(const xt::xtensor<double,4>& A, xt::xtensor<double,4>& Ad)
{
  GMATELASTOPLASTIC_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Ad.shape()[0], Ad.shape()[1], 3, 3}));

  #pragma omp parallel
  {
    Tensor2 I = Cartesian3d::I();
    #pragma omp for
    for (size_t e = 0; e < A.shape()[0]; ++e) {
      for (size_t q = 0; q < A.shape()[1]; ++q) {
        auto Ai  = xt::adapt(&A (e,q,0,0), xt::xshape<3,3>());
        auto Aid = xt::adapt(&Ad(e,q,0,0), xt::xshape<3,3>());
        xt::noalias(Aid) = Ai - trace(Ai) / 3.0 * I;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void epseq(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Aeq)
{
  GMATELASTOPLASTIC_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Aeq.shape()[0], Aeq.shape()[1], 3, 3}));

  #pragma omp parallel
  {
    Tensor2 I = Cartesian3d::I();
    #pragma omp for
    for (size_t e = 0; e < A.shape()[0]; ++e) {
      for (size_t q = 0; q < A.shape()[1]; ++q) {
        auto Ai  = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        auto Aid = Ai - trace(Ai) / 3.0 * I;
        Aeq(e,q) = std::sqrt(2.0/3.0 * A2_ddot_B2(Aid,Aid));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void sigeq(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Aeq)
{
  GMATELASTOPLASTIC_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Aeq.shape()[0], Aeq.shape()[1], 3, 3}));

  #pragma omp parallel
  {
    Tensor2 I = Cartesian3d::I();
    #pragma omp for
    for (size_t e = 0; e < A.shape()[0]; ++e) {
      for (size_t q = 0; q < A.shape()[1]; ++q) {
        auto Ai  = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        auto Aid = Ai - trace(Ai) / 3.0 * I;
        Aeq(e,q) = std::sqrt(1.5 * A2_ddot_B2(Aid,Aid));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Hydrostatic(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Am = xt::empty<double>({A.shape()[0], A.shape()[1]});
  Cartesian3d::hydrostatic(A, Am);
  return Am;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Deviatoric(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,4> Ad = xt::empty<double>(A.shape());
  Cartesian3d::deviatoric(A, Ad);
  return Ad;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Epseq(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Aeq = xt::empty<double>({A.shape()[0], A.shape()[1]});
  Cartesian3d::epseq(A, Aeq);
  return Aeq;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Sigeq(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Aeq = xt::empty<double>({A.shape()[0], A.shape()[1]});
  Cartesian3d::sigeq(A, Aeq);
  return Aeq;
}

// -------------------------------------------------------------------------------------------------

template <class U>
inline double trace(const U& A)
{
  return A(0,0) + A(1,1) + A(2,2);
}

// -------------------------------------------------------------------------------------------------

template <class U, class V, class W>
inline void A2_dyadic_B2(const U& A, const V& B, W& C)
{
  C.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          C(i,j,k,l) += A(i,j) * B(k,l);
}

// -------------------------------------------------------------------------------------------------

template <class U, class V>
inline double A2_ddot_B2(const U& A, const V& B)
{
  return A(0,0) * B(0,0) + 2.0 * A(0,1) * B(0,1) + 2.0 * A(0,2) * B(0,2) +
         A(1,1) * B(1,1) + 2.0 * A(1,2) * B(1,2) +
         A(2,2) * B(2,2);
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
