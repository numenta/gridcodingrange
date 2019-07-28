/* ---------------------------------------------------------------------
 * Numenta Platform for Intelligent Computing (NuPIC)
 * Copyright (C) 2019, Numenta, Inc.  Unless you have an agreement
 * with Numenta, Inc., for a separate license for this software code, the
 * following terms and conditions apply:
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero Public License version 3 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero Public License for more details.
 *
 * You should have received a copy of the GNU Affero Public License
 * along with this program.  If not, see http://www.gnu.org/licenses.
 *
 * http://numenta.org/licenses/
 * ----------------------------------------------------------------------
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <vector>

#include "grid_coding_range.hpp"
#include <nta_logging.hpp>

using std::pair;
using std::vector;

namespace py = pybind11;

static vector<double>
copyArray1D(py::buffer arr)
{
  py::buffer_info info = arr.request();
  NTA_CHECK(info.ndim == 1);

  // We generally call np.asarray in the Python wrapper, so the user can't mess
  // this up.
  NTA_ASSERT(info.itemsize == sizeof(double));
  NTA_ASSERT(info.format == py::format_descriptor<double>::format());

  vector<double> v;
  for (int i = 0; i < info.shape[0]; i++)
  {
    char *pi = (char *)info.ptr + i*info.strides[0];
    v.push_back(*(double*)pi);
  }

  return v;
}

static vector<vector<vector<double>>>
copyArray3D(py::buffer arr)
{
  py::buffer_info info = arr.request();
  NTA_CHECK(info.ndim == 3);

  // We generally call np.asarray in the Python wrapper, so the user can't mess
  // this up.
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
computeCodingRange(
  py::buffer domainToPlaneByModule,
  py::buffer latticeBasisByModule,
  py::buffer scaledbox,
  py::buffer ignorebox,
  double phaseResolution,
  double pingInterval)
{
  return gridcodingrange::computeCodingRange(
    copyArray3D(domainToPlaneByModule), copyArray3D(latticeBasisByModule),
    copyArray1D(scaledbox), copyArray1D(ignorebox), phaseResolution,
    pingInterval);
}

static pair<double, vector<double>>
computeGridUniquenessHypercube(
  py::buffer domainToPlaneByModule,
  py::buffer latticeBasisByModule,
  double phaseResolution,
  double ignoredCenterDiameter,
  double pingInterval)
{
  return gridcodingrange::computeGridUniquenessHypercube(
    copyArray3D(domainToPlaneByModule), copyArray3D(latticeBasisByModule),
    phaseResolution, ignoredCenterDiameter, pingInterval);
}

static double
computeBinSidelength(
  py::buffer domainToPlaneByModule,
  double readoutResolution,
  double resultPrecision,
  double upperBound,
  double timeout)
{
  return gridcodingrange::computeBinSidelength(
    copyArray3D(domainToPlaneByModule), readoutResolution, resultPrecision,
    upperBound, timeout);
}

static vector<double>
computeBinRectangle(
  py::buffer domainToPlaneByModule,
  double readoutResolution,
  double resultPrecision,
  double upperBound,
  double timeout)
{
  return gridcodingrange::computeBinRectangle(
    copyArray3D(domainToPlaneByModule), readoutResolution, resultPrecision,
    upperBound, timeout);
}

PYBIND11_MODULE(_gridcodingrange, m)
{
  m.def("computeCodingRange", &computeCodingRange);
  m.def("computeGridUniquenessHypercube", &computeGridUniquenessHypercube);
  m.def("computeBinSidelength", &computeBinSidelength);
  m.def("computeBinRectangle", &computeBinRectangle);
  m.def("resetCheckPolygonThreshold", &gridcodingrange::resetCheckPolygonThreshold);
  m.def("setCheckPolygonThreshold", &gridcodingrange::setCheckPolygonThreshold);

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
