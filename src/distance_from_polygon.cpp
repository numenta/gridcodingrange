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

#include "distance_from_polygon.hpp"
#include <nta_logging.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>

using std::vector;
using std::pair;

LineSegmentInfo2D::LineSegmentInfo2D(
  const pair<double,double> &_start,
  const pair<double,double> &_end)
  : start(_start), end(_end)
{
  const pair<double,double> v = {end.first - start.first,
                                 end.second - start.second};
  this->length = sqrt(pow(v.first, 2) + pow(v.second, 2));
  this->unitvector = {v.first / length, v.second / length};
}

double distToSegmentSquared(
  pair<double, double> p,
  const LineSegmentInfo2D &segment)
{
  const double distance_along_line =
    (segment.unitvector.first*(p.first - segment.start.first) +
     segment.unitvector.second*(p.second - segment.start.second));

  pair<double,double> nearest_point_on_segment;
  if (distance_along_line <= 0)
  {
    nearest_point_on_segment = segment.start;
  }
  else if (distance_along_line < segment.length)
  {
    nearest_point_on_segment = {
      segment.start.first + segment.unitvector.first * distance_along_line,
      segment.start.second + segment.unitvector.second * distance_along_line,
    };
  }
  else
  {
    nearest_point_on_segment = segment.end;
  }

  return (pow(p.first - nearest_point_on_segment.first, 2) +
          pow(p.second - nearest_point_on_segment.second, 2));
}

double distToSegmentSquared(
  pair<double, double> p,
  pair<double, double> start,
  pair<double, double> end)
{
  return distToSegmentSquared(p, LineSegmentInfo2D(start, end));
}

/**
 * Computing theta can be expensive. Get a number that increases monotonically
 * with theta.
 */
static double getThetaIndex(double dx, double dy)
{
  // // Simple slow version.
  // double theta = atan2(dy, dx);
  // if (theta < 0)
  // {
  //   theta += 2*M_PI;
  // }
  // return theta;

  // Faster version.
  const bool dx_sign = dx >= 0;
  const bool dy_sign = dy >= 0;

  if (dx_sign && dy_sign)
  {
    if (dy < dx)
    {
      return dy / dx;
    }
    else
    {
      return 2.0 - (dx / dy);
    }
  }
  else if (!dx_sign && dy_sign)
  {
    dx = std::abs(dx);
    if (dy > dx)
    {
      return 2.0 + (dx / dy);
    }
    else
    {
      return 4.0 - (dy / dx);
    }
  }
  else if (!dx_sign && !dy_sign)
  {
    dx = std::abs(dx);
    dy = std::abs(dy);
    if (dy < dx)
    {
      return 4.0 + (dy / dx);
    }
    else
    {
      return 6.0 - (dx / dy);
    }
  }
  else
  {
    dy = std::abs(dy);
    if (dy > dx)
    {
      return 6.0 + (dx / dy);
    }
    else
    {
      return 8.0 - (dy / dx);
    }
  }
}

static vector<LineSegmentInfo2D> computeEdges(
  const vector<pair<double,double>> &vertices)
{
  vector<LineSegmentInfo2D> edges;
  edges.reserve(vertices.size());

  for (size_t i = 0; i < vertices.size(); ++i)
  {
    const pair<double,double> &v1 = vertices[i];
    const pair<double,double> &v2 =
      vertices[(i == vertices.size() - 1) ? 0 : i + 1];
    edges.push_back(LineSegmentInfo2D(v1, v2));
  }

  return edges;
}

PolygonInfo::PolygonInfo(
  const vector<pair<double, double>> &vertices)
{
  // Compute the polygon's area times two.
  double acc = 0;
  for (size_t i = 0; i < vertices.size(); ++i)
  {
    const pair<double,double> &v1 = vertices[i];
    const pair<double,double> &v2 =
      vertices[(i == vertices.size() - 1) ? 0 : i + 1];

    acc += v1.first * v2.second;
    acc -= v2.first * v1.second;
  }

  this->is_valid_polygon = (acc != 0);

  if (this->is_valid_polygon)
  {
    // Compute centroid.
    double xsum = 0;
    double ysum = 0;
    for (const pair<double,double> &v : vertices)
    {
      xsum += v.first;
      ysum += v.second;
    }
    this->centroid = {xsum / vertices.size(),
                      ysum / vertices.size()};

    // Compute thetas.
    vector<double> thetas;
    thetas.reserve(vertices.size());
    for (const pair<double,double> &v : vertices)
    {
      thetas.push_back(getThetaIndex(v.first - centroid.first,
                                     v.second - centroid.second));
    }

    // Sort by theta.
    this->vertices.reserve(vertices.size());
    this->thetas.reserve(vertices.size());
    {
      vector<size_t> indices(vertices.size());
      std::iota(indices.begin(), indices.end(), 0);
      std::sort(indices.begin(), indices.end(),
                [&](size_t a, size_t b) {
                  return thetas[a] < thetas[b];
                });
      for (size_t idx : indices)
      {
        this->vertices.push_back(vertices[idx]);
        this->thetas.push_back(thetas[idx]);
      }
    }

    // Compute half planes.
    this->halfplanes.reserve(this->vertices.size());
    for (size_t i = 0; i < this->vertices.size(); ++i)
    {
      const pair<double,double> &v1 = this->vertices[i];
      const pair<double,double> &v2 =
        this->vertices[(i == this->vertices.size() - 1) ? 0 : i + 1];
      const pair<double,double> normalvector = {v2.second - v1.second,
                                                -(v2.first - v1.first)};
      const double top = (normalvector.first*v1.first +
                          normalvector.second*v1.second);
      this->halfplanes.push_back({normalvector, top});
    }

    this->edges = computeEdges(this->vertices);
  }
  else
  {
    bool use_x = false;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      if (vertices[i].first != vertices[i+1].first)
      {
        use_x = true;
        break;
      }
    }

    this->vertices.reserve(2);
    if (use_x)
    {
      auto compare =
        [](const pair<double,double> &a, const pair<double,double> &b)
        {
          return a.first < b.first;
        };

      this->vertices.push_back(
        *std::min_element(vertices.begin(), vertices.end(), compare));
      this->vertices.push_back(
        *std::max_element(vertices.begin(), vertices.end(), compare));
    }
    else
    {
      auto compare =
        [](const pair<double,double> &a, const pair<double,double> &b)
        {
          return a.second < b.second;
        };

      this->vertices.push_back(
        *std::min_element(vertices.begin(), vertices.end(), compare));
      this->vertices.push_back(
        *std::max_element(vertices.begin(), vertices.end(), compare));
    }

    this->centroid = {INFINITY,INFINITY};
    this->thetas = {};
    this->halfplanes = {};
    this->edges = computeEdges(this->vertices);
  }
}

double distToConvexPolygonSquared(
  pair<double, double> point,
  const PolygonInfo &polygon)
{
  if (polygon.is_valid_polygon)
  {
    // Figure out which edge to check.
    const double theta_index = getThetaIndex(
      point.first - polygon.centroid.first,
      point.second - polygon.centroid.second);

    const vector<double>::const_iterator it = (polygon.vertices.size() <= 8)
      ? std::find_if(polygon.thetas.begin(), polygon.thetas.end(),
                     [&](double d)
                     {
                       return theta_index < d;
                     })
      : std::lower_bound(polygon.thetas.begin(), polygon.thetas.end(),
                         theta_index);

    const size_t i = (it == polygon.thetas.begin())
      ? polygon.thetas.size() - 1
      : std::distance(polygon.thetas.begin(), it) - 1;

    // Check whether the lattice point is contained within the polygon.
    const pair<double,double> &normalvector =
      polygon.halfplanes[i].normalvector;
    if (normalvector.first*point.first +
        normalvector.second*point.second
        <= polygon.halfplanes[i].top)
    {
      // The point is contained.
      return 0.0;
    }
  }

  double d = std::numeric_limits<double>::max();
  for (const LineSegmentInfo2D &edge : polygon.edges)
  {
    d = std::min(d, distToSegmentSquared(point, edge));
  }

  return d;
}

double distToConvexPolygonSquared(
  pair<double, double> point,
  const vector<pair<double,double>> &vertices)
{
  return distToConvexPolygonSquared(point, PolygonInfo(vertices));
}
