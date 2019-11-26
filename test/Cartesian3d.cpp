
#include <catch2/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

#include <GMatElastoPlastic/Cartesian3d.h>

namespace GM = GMatElastoPlastic::Cartesian3d;

// =================================================================================================

TEST_CASE("GMatElastoPlastic::Cartesian3d", "Cartesian3d.h")
{

// =================================================================================================

SECTION("Elastic")
{
  // material model
  // - parameters
  double K = 12.3;
  double G = 45.6;
  // - model
  GM::Elastic mat(K, G);

  // allocate tensors
  GM::Tensor2 Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm = 0.12;
  // - strain
  Eps.fill(0.0);
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm;
  Eps(0,1) = Eps(1,0) = gamma;
  // - stress
  Sig = mat.Stress(Eps);
  // - analytical solution
  EQ(Sig(0,0), 3.0 * K * epsm);
  EQ(Sig(1,1), 3.0 * K * epsm);
  EQ(Sig(2,2), 3.0 * K * epsm);
  EQ(Sig(0,1), 2.0 * G * gamma);
  EQ(Sig(1,0), 2.0 * G * gamma);
  EQ(Sig(0,2), 0);
  EQ(Sig(1,2), 0);
  EQ(Sig(2,0), 0);
  EQ(Sig(2,1), 0);
}

// =================================================================================================

SECTION("Matrix")
{
  // parameters
  double K = 12.3;
  double G = 45.6;
  size_t nelem = 3;
  size_t nip = 2;
  size_t d = 3;

  // allocate matrix
  GM::Matrix mat({nelem, nip});

  // all rows elastic
  {
    xt::xtensor<size_t,2> I = xt::ones<size_t>({nelem,nip});
    mat.setElastic(I, K, G);
  }

  // allocate tensors
  GM::Tensor2 Eps;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm = 0.12;
  // - strain
  Eps.fill(0.0);
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm;
  Eps(0,1) = Eps(1,0) = gamma;
  // - strain/stress matrices
  xt::xtensor<double,4> eps = xt::empty<double>({nelem, nip, d, d});
  xt::xtensor<double,4> sig;
  // - set strain
  for (size_t e = 0; e < nelem; ++e)
    for (size_t q = 0; q < nip; ++q)
      xt::view(eps, e, q) = Eps;
  // - stress
  sig = mat.Stress(eps);

  // - analytical solution
  EQ(sig(0,0,0,0), 3.0 * K * epsm);  EQ(sig(0,1,0,0), 3.0 * K * epsm);
  EQ(sig(0,0,1,1), 3.0 * K * epsm);  EQ(sig(0,1,1,1), 3.0 * K * epsm);
  EQ(sig(0,0,2,2), 3.0 * K * epsm);  EQ(sig(0,1,2,2), 3.0 * K * epsm);
  EQ(sig(0,0,0,1), 2.0 * G * gamma); EQ(sig(0,1,0,1), 2.0 * G * gamma);
  EQ(sig(0,0,1,0), 2.0 * G * gamma); EQ(sig(0,1,1,0), 2.0 * G * gamma);
  EQ(sig(1,0,0,0), 3.0 * K * epsm);  EQ(sig(1,1,0,0), 3.0 * K * epsm);
  EQ(sig(1,0,1,1), 3.0 * K * epsm);  EQ(sig(1,1,1,1), 3.0 * K * epsm);
  EQ(sig(1,0,2,2), 3.0 * K * epsm);  EQ(sig(1,1,2,2), 3.0 * K * epsm);
  EQ(sig(1,0,0,1), 2.0 * G * gamma); EQ(sig(1,1,0,1), 2.0 * G * gamma);
  EQ(sig(1,0,1,0), 2.0 * G * gamma); EQ(sig(1,1,1,0), 2.0 * G * gamma);
  REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::keep(0), xt::keep(2)), 0));
  REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::keep(1), xt::keep(2)), 0));
  REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::keep(2), xt::keep(0)), 0));
  REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::keep(2), xt::keep(1)), 0));
}

// =================================================================================================

}
