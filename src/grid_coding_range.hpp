/* ---------------------------------------------------------------------
 * Numenta Platform for Intelligent Computing (NuPIC)
 * Copyright (C) 2018-2019, Numenta, Inc.  Unless you have an agreement
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

/** @file
 * Functions that analyze the uniqueness of locations in grid cell spaces
 */

#ifndef NTA_GRIDCODINGRANGE
#define NTA_GRIDCODINGRANGE

#include <vector>
#include <utility>

namespace gridcodingrange
{
  /**
   * Determine whether any points in a k-dimensional rectangle have a grid code
   * equal to the grid code at the origin.
   *
   * If this function says there's a grid code zero at point (42.0, 42.0), that
   * implies that every location has the same grid code as the location that is
   * 42.0 units up and 42.0 units right.
   *
   * This function uses a recursive "divide and conquer" algorithm. It first
   * quickly samples a few points to see if any of them have grid code zero.
   * Then it tries to prove that grid code zero can't possibly happen in this
   * hyperrectangle by showing that at least one module never has grid code zero
   * for any point in this hyperrectangle. Finally, if neither attempt succeeds,
   * it divides the hyperrectangle in half, and then it tries again on both
   * halves.
   *
   * @param x0
   * The lowest corner of the k-dimensional rectangle that will be searched.
   *
   * @param dims
   * The dimensions of the k-dimensional rectangle that will be searched.
   *
   * @param pointWithGridCodeZero
   * Output parameter. A point with grid code zero in this hyperrectangle. Only
   * populated if the function returns true.
   *
   * @return
   * true if grid code zero is found, false otherwise.
   */
  bool findGridCodeZero(
      const std::vector<std::vector<std::vector<double>>> &domainToPlaneByModule,
      const std::vector<std::vector<std::vector<double>>> &latticeBasisByModule,
      const std::vector<double> &x0,
      const std::vector<double> &dims,
      double readoutResolution,
      std::vector<double> *pointWithGridCodeZero = nullptr);

  /**
   * Given a set of grid cell module parameters, scale a k-dimensional box until
   * it reaches a point with the same grid cell representation as the origin.
   *
   * This function expands scaledboxes outward from the origin, one in each
   * high-dimensional quadrant, by multiplying it with a scaling factor and
   * iteratively building out the larger box out of smaller boxes, relying on
   * findGridCodeZero to analyze each box. When the expanding box collides with
   * a point that has the same grid cell representation as the origin, the
   * scaling factor from before the collision is returned.
   *
   * Each grid cell module is specified by a pair of matrices. The first one
   * projects a k-dimensional point to a 2D plane, and the second matrix
   * specifies the basis of the grid lattice on this plane, specifying which
   * points on the plane have equivalent grid codes. The distance between two
   * grid codes is equal to the shortest distance between them on the plane. The
   * readout resolution parameter is measured in units of this distance.
   *
   * There's not a strictly "correct" way to configure these matrices and the
   * readout resolution, but a typical way is to use unit vectors as lattice
   * basis vectors, use a plane projection that normalizes distances so that 1
   * unit is 1 "scale", and then set the readout resolution to some number in
   * the range (0, 1) so that it is measured in units of "scales".
   *
   * @param domainToPlaneByModule
   * A list of 2*k matrices, one per module. The matrix converts from a point in
   * the domain to a point on a plane, normalizing for grid cell scale.
   *
   * @param latticeBasisByModule
   * A list of m 2*2 matrices, one per module. This matrix contains the basis
   * vectors for a lattice, specifying which points on the plane have equivalent
   * location representations in this module.
   *
   * @param scaledbox
   * The box that is scaled outward in every high-dimensional quadrant.
   *
   * @param ignorebox
   * The box that should be ignored because it contains the *actual* grid code
   * zero. This hypercube starts at the origin and extends in the positive
   * direction. It is applied in every high-dimensional quadrant.
   *
   * @param readoutResolution
   * The precision of readout of this grid code, measured in distance on the
   * plane. For example, if this is 0.2, then all points on the plane are
   * indistinguishable from those in their surrounding +- 0.1 range.
   *
   * @param pingInterval
   * How often, in seconds, the function should print its current status. If <=
   * 0, no printing will occur.
   *
   * @return
   * - The largest tested scaling factor of the scaledbox that contains no
       collisions.
   * - A point just outside this scaled scaledbox that collides with the origin.
   */
  std::pair<double, std::vector<double>> computeCodingRange(
      const std::vector<std::vector<std::vector<double>>> &domainToPlaneByModule,
      const std::vector<std::vector<std::vector<double>>> &latticeBasisByModule,
      const std::vector<double> &scaledbox,
      const std::vector<double> &ignorebox,
      double readoutResolution,
      double pingInterval = 10.0);

  /**
   * Calls computeCodingRange with a unit cube scaledBox and cube ignore box.
   *
   * @param ignoredCenterDiameter
   * The sidelength of the hypercube which should be ignored because it contains
   * the *actual* grid code zero. This hypercube starts at the origin and
   * extends in the positive direction. It is applied in every expansion
   * direction (i.e. in every quadrant/octant/orthant).
   *
   * @param pingInterval
   * How often, in seconds, the function should print its current status. If <=
   * 0, no printing will occur.
   *
   * @return
   * - The diameter of the hypercube that contains no collisions.
   * - A point just outside this hypercube that collides with the origin.
   */
  std::pair<double, std::vector<double>> computeGridUniquenessHypercube(
      const std::vector<std::vector<std::vector<double>>> &domainToPlaneByModule,
      const std::vector<std::vector<std::vector<double>>> &latticeBasisByModule,
      double readoutResolution,
      double ignoredCenterDiameter,
      double pingInterval = 10.0);

  /**
   * Compute the sidelength of the smallest hypercube that encloses the
   * intersection of all of the modules' firing fields centered at the origin.
   * This sidelength is like a resolution of the code.
   *
   * @param domainToPlaneByModule
   * A list of 2*k matrices, one per module. The matrix converts from a point in
   * the domain to a point on a plane, normalizing for grid cell scale.
   *
   * @param readoutResolution
   * The precision of readout of this grid code, measured in distance on the
   * plane. For example, if this is 0.2, then all points on the plane are
   * indistinguishable from those in their surrounding +- 0.1 range.
   *
   * @param resultPrecision
   * The precision level for this sidelength. This algorithm will binary-search
   * for the smallest hypercube it can place around zero, stopping at the point
   * when it's within 'resultPrecision' of the actual smallest cube.
   *
   * @param upperBound
   * Instructs the algorithm when it should give up. If the provided matrix
   * never moves away from grid code zero, the algorithm could search forever.
   * This parameter tells it when it should stop searching.
   *
   * @param timeout
   * Specifies how long to try. This function will give up after 'timeout'
   * seconds. If <= 0, the function will not time out. On timeout, the function
   * throws an exception with message "timeout". In Python this exception is of
   * type RuntimeError.
   *
   * @return
   * The sidelength of this hypercube. Returns -1.0 if a surface can't be found
   * (i.e. if upperBound is reached.)
   */
  double computeBinSidelength(
      const std::vector<std::vector<std::vector<double>>> &domainToPlaneByModule,
      double readoutResolution,
      double resultPrecision,
      double upperBound = 2048.0,
      double timeout = -1.0);

  /**
   * Like computeBinSidelength, but it computes a hyperrectangle rather than a
   * hypercube.
   *
   * @param domainToPlaneByModule
   * A list of 2*k matrices, one per module. The matrix converts from a point in
   * the domain to a point on a plane, normalizing for grid cell scale.
   *
   * @param readoutResolution
   * The precision of readout of this grid code, measured in distance on the
   * plane. For example, if this is 0.2, then all points on the plane are
   * indistinguishable from those in their surrounding +- 0.1 range.
   *
   * @param resultPrecision
   * The precision level for this sidelength. This algorithm will binary-search
   * for the smallest hyperrectangle it can place around zero, stopping at the
   * point when each side is within 'resultPrecision' of the actual smallest
   * hyperrectangle.
   *
   * @param upperBound
   * Instructs the algorithm when it should give up. If the provided matrix
   * never moves away from grid code zero, the algorithm could search forever.
   * This parameter tells it when it should stop searching.
   *
   * @param timeout
   * Specifies how long to try. This function will give up after 'timeout'
   * seconds. If <= 0, the function will not time out. On timeout, the function
   * throws an exception with message "timeout". In Python this exception is of
   * type RuntimeError.
   *
   * @return
   * The dimensions of this hyperrectangle. Returns an empty vector if a surface
   * can't be found (i.e. if upperBound is reached.)
   */
  std::vector<double> computeBinRectangle(
      const std::vector<std::vector<std::vector<double>>> &domainToPlaneByModule,
      double readoutResolution,
      double resultPrecision,
      double upperBound = 2048.0,
      double timeout = -1.0);

} // end namespace gridcodingrange

#endif // NTA_GRIDCODINGRANGE
