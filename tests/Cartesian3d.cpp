#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <GMatElastoPlastic/Cartesian3d.h>
#include <GMatTensor/Cartesian3d.h>

namespace GM = GMatElastoPlastic::Cartesian3d;
namespace GT = GMatTensor::Cartesian3d;

TEST_CASE("GMatElastoPlastic::Cartesian3d", "Cartesian3d.h")
{
    SECTION("Elastic")
    {
        double K = 12.3;
        double G = 45.6;
        double gamma = 0.02;
        double epsm = 0.12;

        xt::xtensor<double, 2> Eps = {
            {epsm, gamma, 0.0},
            {gamma, epsm, 0.0},
            {0.0, 0.0, epsm}};

        xt::xtensor<double, 2> Sig = {
            {3.0 * K * epsm, 2.0 * G * gamma, 0.0},
            {2.0 * G * gamma, 3.0 * K * epsm, 0.0},
            {0.0, 0.0, 3.0 * K * epsm}};

        GM::Elastic mat(K, G);
        mat.setStrain(Eps);

        REQUIRE(xt::allclose(mat.Stress(), Sig));
    }

    SECTION("Tangent (purely elastic response only) - Elastic")
    {
        double K = 12.3;
        double G = 45.6;

        xt::xtensor<double, 2> Eps = xt::random::randn<double>({3, 3});
        xt::xtensor<double, 4> Is = GM::I4s();
        Eps = GT::A4_ddot_B2(Is, Eps);

        GM::Elastic mat(K, G);
        mat.setStrain(Eps);
        auto Sig = mat.Stress();
        auto C = mat.Tangent();
        REQUIRE(xt::allclose(GT::A4_ddot_B2(C, Eps), Sig));
    }

    SECTION("Array")
    {
        double K = 12.3;
        double G = 45.6;
        double gamma = 0.02;
        double epsm = 0.12;

        xt::xtensor<double, 2> Eps = {
            {epsm, gamma, 0.0},
            {gamma, epsm, 0.0},
            {0.0, 0.0, epsm}};

        xt::xtensor<double, 2> Sig = {
            {3.0 * K * epsm, 2.0 * G * gamma, 0.0},
            {2.0 * G * gamma, 3.0 * K * epsm, 0.0},
            {0.0, 0.0, 3.0 * K * epsm}};

        size_t nelem = 3;
        size_t nip = 2;
        size_t ndim = 3;

        GM::Array<2> mat({nelem, nip});

        {
            xt::xtensor<size_t,2> I = xt::ones<size_t>({nelem, nip});
            mat.setElastic(I, K, G);
        }

        xt::xtensor<double, 4> eps = xt::empty<double>({nelem, nip, ndim, ndim});
        xt::xtensor<double, 4> sig = xt::empty<double>({nelem, nip, ndim, ndim});

        for (size_t e = 0; e < nelem; ++e) {
            for (size_t q = 0; q < nip; ++q) {
                double fac = static_cast<double>((e + 1) * nip + (q + 1));
                xt::view(eps, e, q) = fac * Eps;
                xt::view(sig, e, q) = fac * Sig;
            }
        }

        mat.setStrain(eps);

        REQUIRE(xt::allclose(mat.Stress(), sig));
    }

    SECTION("Array - Model")
    {
        double K = 12.3;
        double G = 45.6;
        double gamma = 0.02;
        double epsm = 0.12;

        xt::xtensor<double, 2> Eps = {
            {epsm, gamma, 0.0},
            {gamma, epsm, 0.0},
            {0.0, 0.0, epsm}};

        xt::xtensor<double, 2> Sig = {
            {3.0 * K * epsm, 2.0 * G * gamma, 0.0},
            {2.0 * G * gamma, 3.0 * K * epsm, 0.0},
            {0.0, 0.0, 3.0 * K * epsm}};

        size_t nelem = 3;
        size_t nip = 2;

        GM::Array<2> mat({nelem, nip});

        {
            xt::xtensor<size_t,2> I = xt::ones<size_t>({nelem, nip});
            mat.setElastic(I, K, G);
        }

        for (size_t e = 0; e < nelem; ++e) {
            for (size_t q = 0; q < nip; ++q) {
                double fac = static_cast<double>((e + 1) * nip + (q + 1));
                auto model = mat.getElastic({e, q});
                model.setStrain(xt::eval(fac * Eps));
                REQUIRE(xt::allclose(model.Stress(), fac * Sig));
            }
        }
    }
}
