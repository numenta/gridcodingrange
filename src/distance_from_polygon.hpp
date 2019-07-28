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

#ifndef NTA_DISTANCE_FROM_POLYGON_HPP
#define NTA_DISTANCE_FROM_POLYGON_HPP

#include <utility>
#include <vector>

struct LineSegmentInfo2D {
  LineSegmentInfo2D(
    const std::pair<double,double> &start,
    const std::pair<double,double> &end);
  std::pair<double,double> start;
  std::pair<double,double> end;
  std::pair<double,double> unitvector;
  double length;
};

struct HalfplaneInfo
{
  std::pair<double,double> normalvector;
  double top;
};

struct PolygonInfo {
  PolygonInfo() {}
  PolygonInfo(
    const std::vector<std::pair<double, double>> &vertices);

  bool is_valid_polygon;
  std::pair<double,double> centroid;
  std::vector<std::pair<double,double>> vertices;
  std::vector<double> thetas;
  std::vector<HalfplaneInfo> halfplanes;
  std::vector<LineSegmentInfo2D> edges;
};


double distToSegmentSquared(
  std::pair<double, double> point,
  const LineSegmentInfo2D &segment);

double distToSegmentSquared(
  std::pair<double, double> point,
  std::pair<double, double> start,
  std::pair<double, double> end);


/**
 * Measure the squared distance from a point to a convex polygon.
 *
 * This function is optimized for the case of testing the same polygon with many
 * different points. The PolygonInfo contains cached information that would
 * normally be too expensive to justify computing for a single operation, but if
 * this polygon is going to be tested many times then it is worth it.
 */
double distToConvexPolygonSquared(
  std::pair<double, double> point,
  const PolygonInfo& polygon);


/**
 * The slow version of distToConvexPolygonSquared. Instantiate your own
 * PolygonInfo and cache it somewhere to get the performance benefits.
 */
double distToConvexPolygonSquared(
  std::pair<double, double> point,
  const std::vector<std::pair<double, double>> &vertices);

#endif // NTA_DISTANCE_FROM_POLYGON_HPP
