#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <vector>

#include "GridUniqueness.hpp"
#include <nta_logging.hpp>

using std::pair;
using std::vector;

namespace py = pybind11;

static vector<vector<vector<double>>>
copyArray3D(py::buffer arr)
{
    py::buffer_info info = arr.request();
    NTA_CHECK(info.ndim == 3);

    // We generally call np.asarray in the Python wrapper, so the user can't 
    // mess this up.
    NTA_ASSERT(info.itemsize == sizeof(double));
    NTA_ASSERT(info.format == py::format_descriptor<double>::format());

    vector<vector<vector<double>>> m;
    for (int i = 0; i < info.shape[0]; i++)
    {
        vector<vector<double>> module;
        char *pi = (char *)info.ptr + i*info.strides[0];
        for (int j = 0; j < info.shape[1]; j++)
        {
            vector<double> row;
            char *pj = pi + j*info.strides[1];
            for (int k = 0; k < info.shape[2]; k++)
            {
                char *pk = pj + k*info.strides[2];
                row.push_back(*(double *)pk);
            }
            module.push_back(row);
        }
        m.push_back(module);
    }

    return m;
}

static pair<double, vector<double>>
computeGridUniquenessHypercube(
    py::buffer domainToPlaneByModule,
    py::buffer latticeBasisByModule,
    double phaseResolution,
    double ignoredCenterDiameter,
    double pingInterval = 10.0)
{
    return nupic::experimental::grid_uniqueness::computeGridUniquenessHypercube(
        copyArray3D(domainToPlaneByModule), copyArray3D(latticeBasisByModule),
        phaseResolution, ignoredCenterDiameter, pingInterval);
}

static double
computeBinSidelength(
    py::buffer domainToPlaneByModule,
    double readoutResolution,
    double resultPrecision,
    double upperBound = 1000.0,
    double timeout = -1.0)
{
    return nupic::experimental::grid_uniqueness::computeBinSidelength(
        copyArray3D(domainToPlaneByModule), readoutResolution, resultPrecision,
        upperBound, timeout);
}

static vector<double>
computeBinRectangle(
    py::buffer domainToPlaneByModule,
    double readoutResolution,
    double resultPrecision,
    double upperBound = 1000.0,
    double timeout = -1.0)
{
    return nupic::experimental::grid_uniqueness::computeBinRectangle(
        copyArray3D(domainToPlaneByModule), readoutResolution, resultPrecision,
        upperBound, timeout);
}

PYBIND11_MODULE(_gridcodingrange, m)
{
    m.def("computeGridUniquenessHypercube", &computeGridUniquenessHypercube);
    m.def("computeBinSidelength", &computeBinSidelength);
    m.def("computeBinRectangle", &computeBinRectangle);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
