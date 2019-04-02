/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

================================================================================================= */

#ifndef GMATELASTOPLASTIC_CARTESIAN3D_H
#define GMATELASTOPLASTIC_CARTESIAN3D_H

#include "config.h"

namespace GMatElastoPlastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

// Alias

using T2 = xt::xtensor_fixed<double, xt::xshape<3,3>>;
using T4 = xt::xtensor_fixed<double, xt::xshape<3,3,3,3>>;

// Tensor operations

template <class T> inline double trace (const T& A);
template <class T> inline double ddot22(const T& A, const T& B);

template <class T2, class T4>
inline void dyadic22(const T2& A, const T2& B, T4& C);

// Unit tensors

inline T2 I();
inline T4 II();
inline T4 I4();
inline T4 I4rt();
inline T4 I4s();
inline T4 I4d();

// -------------------------------------------------------------------------------------------------

// Hydrostatic stress/strain

inline double Hydrostatic(const T2& A);

// Deviatoric part of a tensor

inline T2 Deviatoric(const T2& A);

// Equivalent deviatoric stress/stress

inline double Sigeq(const T2& Sig);
inline double Epseq(const T2& Eps);

// Matrix version of the functions above (no allocation)

inline void hydrostatic(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Am);
inline void deviatoric(const xt::xtensor<double,4>& A, xt::xtensor<double,4>& Ad);
inline void sigeq(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Aeq);
inline void epseq(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Aeq);

// Auto-allocation allocation of the functions above

inline xt::xtensor<double,2> Hydrostatic(const xt::xtensor<double,4>& A);
inline xt::xtensor<double,4> Deviatoric(const xt::xtensor<double,4>& A);
inline xt::xtensor<double,2> Sigeq(const xt::xtensor<double,4>& Sig);
inline xt::xtensor<double,2> Epseq(const xt::xtensor<double,4>& Eps);

// -------------------------------------------------------------------------------------------------

class Elastic
{
public:

  // Constructors
  Elastic() = default;
  Elastic(double K, double G);

  // Parameters
  double K() const;
  double G() const;

  // Stress (no allocation, overwrites "Sig")
  template <class T>
  void stress(const T2& Eps, T&& Sig) const;

  // Stress (auto allocation)
  T2 Stress(const T2& Eps) const;

  // Stress & Tangent (no allocation, overwrites "Sig" and "C")
  template <class T, class S>
  void tangent(const T2& Eps, T&& Sig, S&& C) const;

  // Stress & Tangent (auto allocation)
  std::tuple<T2,T4> Tangent(const T2& Eps) const;

  // Increment (does nothing)
  void increment() const;

  // Plastic strain (always zero)
  double epsp() const;

private:

  double m_K; // bulk modulus
  double m_G; // shear modulus
};

// -------------------------------------------------------------------------------------------------

class LinearHardening
{
public:

  // Constructors
  LinearHardening() = default;
  LinearHardening(double K, double G, double sigy0, double H);

  // Parameters
  double K() const;
  double G() const;
  double sigy0() const;
  double H() const;

  // Stress (no allocation, overwrites "Sig")
  template <class T>
  void stress(const T2& Eps, T&& Sig);

  // Stress (auto allocation)
  T2 Stress(const T2& Eps);

  // Stress & Tangent (no allocation, overwrites "Sig" and "C")
  template <class T, class S>
  void tangent(const T2& Eps, T&& Sig, S&& C);

  // Stress & Tangent (auto allocation)
  std::tuple<T2,T4> Tangent(const T2& Eps);

  // Increment (update plastic strain)
  void increment();

  // Plastic strain
  double epsp() const;

private:

  double m_K;     // bulk modulus
  double m_G;     // shear modulus
  double m_sigy0; // initial shear stress
  double m_H;     // hardening modulus

  // plastic strain
  double m_epsp=0.0;
  double m_epsp_t=0.0;

  // elastic strain tensor
  T2 m_Epse=T2({{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}});
  T2 m_Epse_t=T2({{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}});

  // strain tensor
  T2 m_Eps=T2({{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}});
  T2 m_Eps_t=T2({{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}});
};

// -------------------------------------------------------------------------------------------------

struct Type {
  enum Value {
    Unset,
    Elastic,
    LinearHardening,
  };
};

// -------------------------------------------------------------------------------------------------

class Matrix
{
public:

  // Constructors

  Matrix() = default;
  Matrix(size_t nelem, size_t nip);

  // Shape

  size_t nelem() const;
  size_t nip() const;

  // Type

  xt::xtensor<size_t,2> type() const;
  xt::xtensor<size_t,2> isPlastic() const;

  // Parameters

  xt::xtensor<double,2> K() const;
  xt::xtensor<double,2> G() const;

  // Matrix of unit tensors

  xt::xtensor<double,4> I() const;
  xt::xtensor<double,6> II() const;
  xt::xtensor<double,6> I4() const;
  xt::xtensor<double,6> I4rt() const;
  xt::xtensor<double,6> I4s() const;
  xt::xtensor<double,6> I4d() const;

  // Check that a type has been set everywhere (throws if unset points are found)

  void check() const;

  // Set parameters for a batch of points

  void setElastic(
    const xt::xtensor<size_t,2>& I,
    double K,
    double G);

  void setLinearHardening(
    const xt::xtensor<size_t,2>& I,
    double K,
    double G,
    double sigy0,
    double H);

  // Set material definition for a batch of points

  void setElastic(
    const xt::xtensor<size_t,2>& I,
    const xt::xtensor<size_t,2>& idx,
    const xt::xtensor<double,1>& K,
    const xt::xtensor<double,1>& G);

  void setLinearHardening(
    const xt::xtensor<size_t,2>& I,
    const xt::xtensor<size_t,2>& idx,
    const xt::xtensor<double,1>& K,
    const xt::xtensor<double,1>& G,
    const xt::xtensor<double,1>& sigy0,
    const xt::xtensor<double,1>& H);

  // Compute (no allocation, overwrites last argument)

  void stress(
    const xt::xtensor<double,4>& Eps,
          xt::xtensor<double,4>& Sig);

  void tangent(
    const xt::xtensor<double,4>& Eps,
          xt::xtensor<double,4>& Sig,
          xt::xtensor<double,6>& C);

  void epsp(
          xt::xtensor<double,2>& epsp) const;

  // Auto-allocation of the functions above

  xt::xtensor<double,4> Stress(
    const xt::xtensor<double,4>& Eps);

  std::tuple<xt::xtensor<double,4>,xt::xtensor<double,6>> Tangent(
    const xt::xtensor<double,4>& Eps);

  xt::xtensor<double,4> Epsp() const;

private:

  // Material vectors
  std::vector<Elastic> m_Elastic;
  std::vector<LinearHardening> m_LinearHardening;

  // Identifiers for each matrix entry
  xt::xtensor<size_t,2> m_type;  // type (e.g. "Type::Elastic")
  xt::xtensor<size_t,2> m_index; // index from the relevant material vector (e.g. "m_Elastic")

  // Shape
  size_t m_nelem;
  size_t m_nip;
  static const size_t m_ndim=3;

  // Internal check
  bool m_allSet=false;
  void checkAllSet();

};

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#include "Cartesian3d.hpp"
#include "Cartesian3d_Elastic.hpp"
#include "Cartesian3d_LinearHardening.hpp"
#include "Cartesian3d_Matrix.hpp"

// -------------------------------------------------------------------------------------------------

#endif
