/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

*/

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

// Enable basic assertions on matrix shape
// (doesn't cost a lot of time, but avoids segmentation faults)
#define GMATELASTOPLASTIC_ENABLE_ASSERT

#include <GMatElastoPlastic/Cartesian3d.h>

namespace py = pybind11;

PYBIND11_MODULE(GMatElastoPlastic, m)
{

m.doc() = "Elasto-plastic material model";

// ---------------------------------------------
// GMatElastoPlasticFiniteStrainSimo.Cartesian3d
// ---------------------------------------------

py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

namespace SM = GMatElastoPlastic::Cartesian3d;

// Unit tensors

sm.def("I2", &SM::I2, "Second order unit tensor.");

sm.def("II", &SM::II, "Fourth order tensor with the result of the dyadic product II.");

sm.def("I4", &SM::I4, "Fourth order unit tensor.");

sm.def("I4rt", &SM::I4rt, "Fourth right-transposed order unit tensor.");

sm.def("I4s", &SM::I4s, "Fourth order symmetric projection tensor.");

sm.def("I4d", &SM::I4d, "Fourth order deviatoric projection tensor.");

// Tensor algebra

sm.def("Hydrostatic",
    py::overload_cast<const SM::Tensor2&>(&SM::Hydrostatic),
    "Hydrostatic part of a 2nd-order tensor. Returns scalar.",
    py::arg("A"));

sm.def("Hydrostatic",
    py::overload_cast<const xt::xtensor<double,3>&>(&SM::Hydrostatic),
    "Hydrostatic part of a 2nd-order tensor. Returns list of scalars.",
    py::arg("A"));

sm.def("Hydrostatic",
    py::overload_cast<const xt::xtensor<double,4>&>(&SM::Hydrostatic),
    "Hydrostatic part of a 2nd-order tensor. Returns matrix of scalars.",
    py::arg("A"));

sm.def("Deviatoric",
    py::overload_cast<const SM::Tensor2&>(&SM::Deviatoric),
    "Deviatoric part of a 2nd-order tensor. Returns 2nd-order tensor.",
    py::arg("A"));

sm.def("Deviatoric",
    py::overload_cast<const xt::xtensor<double,3>&>(&SM::Deviatoric),
    "Deviatoric part of a 2nd-order tensor. Returns list 2nd-order tensors.",
    py::arg("A"));

sm.def("Deviatoric",
    py::overload_cast<const xt::xtensor<double,4>&>(&SM::Deviatoric),
    "Deviatoric part of a 2nd-order tensor. Returns matrix 2nd-order tensors.",
    py::arg("A"));

sm.def("Epseq",
    py::overload_cast<const SM::Tensor2&>(&SM::Epseq),
    "Equivalent strain deviator. Returns scalar.",
    py::arg("Eps"));

sm.def("Epseq",
    py::overload_cast<const xt::xtensor<double,3>&>(&SM::Epseq),
    "Equivalent strain deviator. Returns list of scalars.",
    py::arg("Eps"));

sm.def("Epseq",
    py::overload_cast<const xt::xtensor<double,4>&>(&SM::Epseq),
    "Equivalent strain deviator. Returns matrix of scalars.",
    py::arg("Eps"));

sm.def("Sigeq",
    py::overload_cast<const SM::Tensor2&>(&SM::Sigeq),
    "Equivalent stress deviator. Returns scalar.",
    py::arg("Sig"));

sm.def("Sigeq",
    py::overload_cast<const xt::xtensor<double,3>&>(&SM::Sigeq),
    "Equivalent stress deviator. Returns list of scalars.",
    py::arg("Sig"));

sm.def("Sigeq",
    py::overload_cast<const xt::xtensor<double,4>&>(&SM::Sigeq),
    "Equivalent stress deviator. Returns matrix of scalars.",
    py::arg("Sig"));

// Material point: Elastic

py::class_<SM::Elastic>(sm, "Elastic")

    .def(py::init<double, double>(), "Linear elastic material point.", py::arg("K"), py::arg("G"))

    .def("K", &SM::Elastic::K, "Returns the bulk modulus.")

    .def("G", &SM::Elastic::G, "Returns the shear modulus.")

    .def("Stress",
        &SM::Elastic::Stress,
        "Returns stress tensor, for a given strain tensor.",
        py::arg("Eps"))

    .def("Tangent",
        &SM::Elastic::Tangent,
        "Returns stress and tangent stiffness tensors, for a given strain tensor.",
        py::arg("Eps"))

    .def("__repr__", [](const SM::Elastic&) { return "<GMatElastoPlastic.Cartesian3d.Elastic>"; });

// Material point: LinearHardening

py::class_<SM::LinearHardening>(sm, "LinearHardening")

    .def(py::init<double, double, double, double>(),
        "Elasto-plastic material point.",
        py::arg("K"),
        py::arg("G"),
        py::arg("epsy0"),
        py::arg("H"))

    .def("K", &SM::LinearHardening::K, "Returns the bulk modulus.")

    .def("G", &SM::LinearHardening::G, "Returns the shear modulus.")

    .def("sigy0", &SM::LinearHardening::sigy0, "Returns the initial yield stress.")

    .def("H", &SM::LinearHardening::H, "Returns the hardening modulus.")

    .def("Stress",
        &SM::LinearHardening::Stress,
        "Returns stress tensor, for a given strain tensor.",
        py::arg("Eps"))

    .def("Tangent",
        &SM::LinearHardening::Tangent,
        "Returns stress and tangent stiffness tensors, for a given strain tensor",
        py::arg("Eps"))

    .def("increment",
        &SM::LinearHardening::increment,
        "Apply time-increment: updates the internal history variables.")

    .def("epsp", &SM::LinearHardening::epsp, "Returns the current equivalent plastic strain.")

    .def("__repr__", [](const SM::LinearHardening&) {
        return "<GMatElastoPlastic.Cartesian3d.LinearHardening>";
    });

// Material identifier

py::module smm = sm.def_submodule("Type", "Type enumerator");

py::enum_<SM::Type::Value>(smm, "Type")
    .value("Unset", SM::Type::Unset)
    .value("Elastic", SM::Type::Elastic)
    .value("LinearHardening", SM::Type::LinearHardening)
    .export_values();

// Matrix

py::class_<SM::Matrix>(sm, "Matrix")

    .def(py::init<size_t, size_t>(),
        "Matrix of linear elastic material points",
        py::arg("nelem"),
        py::arg("nip"))

    .def("ndim", &SM::Matrix::ndim, "Return number of (tensor) dimensions.")

    .def("nelem", &SM::Matrix::nelem, "Return number of elements (matrix rows).")

    .def("nip", &SM::Matrix::nip, "Return number of integration points (matrix columns).")

    .def("K", &SM::Matrix::K, "Return matrix with bulk moduli.")

    .def("G", &SM::Matrix::G, "Return matrix with shear moduli.")

    .def("I2", &SM::Matrix::I2, "Return matrix with second order unit tensors.")

    .def("II",
        &SM::Matrix::II,
        "Return matrix with fourth order tensors with the result of the dyadic product II.")

    .def("I4", &SM::Matrix::I4, "Return matrix with fourth order unit tensors.")

    .def("I4rt",
        &SM::Matrix::I4rt,
        "Return matrix with fourth right-transposed order unit tensors.")

    .def("I4s",
        &SM::Matrix::I4s,
        "Return matrix with fourth order symmetric projection tensors.")

    .def("I4d",
        &SM::Matrix::I4d,
        "Return matrix with fourth order deviatoric projection tensors.")

    .def("type", &SM::Matrix::type, "Return matrix with material types.")

    .def("isPlastic",
        &SM::Matrix::isPlastic,
        "Return matrix with boolean: elastic (0) or plastic (1).")

    .def("check",
        &SM::Matrix::check,
        "Check that all matrix entries are set. Throws if any unset point is found.")

    .def("setElastic",
        py::overload_cast<const xt::xtensor<size_t,2>&, double, double>(
            &SM::Matrix::setElastic),
        "Set specific entries 'Elastic'.",
        py::arg("I"),
        py::arg("K"),
        py::arg("G"))

    .def("setLinearHardening",
        py::overload_cast<const xt::xtensor<size_t,2>&, double, double, double, double>(
            &SM::Matrix::setLinearHardening),
        "Set specific entries 'LinearHardening'.",
        py::arg("I"),
        py::arg("K"),
        py::arg("G"),
        py::arg("sigy0"),
        py::arg("H"))

    .def("setElastic",
        py::overload_cast<
            const xt::xtensor<size_t,2>&,
            const xt::xtensor<size_t,2>&,
            const xt::xtensor<double,1>&,
            const xt::xtensor<double,1>&>(&SM::Matrix::setElastic),
        "Set specific entries 'Elastic'.",
        py::arg("I"),
        py::arg("idx"),
        py::arg("K"),
        py::arg("G"))

    .def("setLinearHardening",
        py::overload_cast<
            const xt::xtensor<size_t,2>&,
            const xt::xtensor<size_t,2>&,
            const xt::xtensor<double,1>&,
            const xt::xtensor<double,1>&,
            const xt::xtensor<double,1>&,
            const xt::xtensor<double,1>&>(&SM::Matrix::setLinearHardening),
        "Set specific entries 'LinearHardening'.",
        py::arg("I"),
        py::arg("idx"),
        py::arg("K"),
        py::arg("G"),
        py::arg("sigy0"),
        py::arg("H"))

    .def("increment",
        &SM::Matrix::increment,
        "Apply time-increment: updates the internal history variables.")

    .def("Stress",
        &SM::Matrix::Stress,
        "Returns matrix of stress tensors, for a given matrix of strain tensors.",
        py::arg("Eps"))

    .def("Tangent",
        &SM::Matrix::Tangent,
        "Returns matrices of stress tangent stiffness tensors, "
        "for a given matrix of strain tensors.",
        py::arg("Eps"))

    .def("Epsp", &SM::Matrix::Epsp, "Returns matrix of current equivalent plastic strains.")

    .def("__repr__", [](const SM::Matrix&) { return "<GMatElastoPlastic.Cartesian3d.Matrix>"; });
}
