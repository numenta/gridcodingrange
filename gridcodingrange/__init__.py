import _gridcodingrange

import numpy as np

def computeCodingRange(domainToPlaneByModule, latticeBasisByModule,
                       boxToScale, ignoreBox, phaseResolution,
                       pingInterval=10.0):
    '''
    Given a set of grid cell module parameters, scale a k-dimensional box until
    it reaches a point with the same grid cell representation as the origin.

    This function expands scaledboxes outward from the origin, one in each
    high-dimensional quadrant, by multiplying it with a scaling factor and
    iteratively building out the larger box out of smaller boxes and scanning
    each small box using a divide-and-conquer algorithm. When the expanding box
    collides with a point that has the same grid cell representation as the
    origin, the scaling factor from before the collision is returned.

    Each grid cell module is specified by a pair of matrices. The first one
    projects a k-dimensional point to a 2D plane, and the second matrix
    specifies the basis of the grid lattice on this plane, specifying which
    points on the plane have equivalent grid codes. The distance between two
    grid codes is equal to the shortest distance between them on the plane. The
    readout resolution parameter is measured in units of this distance.

    There's not a strictly "correct" way to configure these matrices and the
    readout resolution, but a typical way is to use unit vectors as lattice
    basis vectors, use a plane projection that normalizes distances so that 1
    unit is 1 "scale", and then set the readout resolution to some number in the
    range (0, 1) so that it is measured in units of "scales".

    @param domainToPlaneByModule (3D numpy array)
    A list of 2*k matrices, one per module. The matrix converts from a point in
    the domain to a point on a plane, normalizing for grid cell scale.

    @param latticeBasisByModule (3D numpy array)
    A list of m 2*2 matrices, one per module. This matrix contains the basis
    vectors for a lattice, specifying which points on the plane have equivalent
    location representations in this module.

    @param scaledbox (1D numpy array)
    The box that is scaled outward in every high-dimensional quadrant.

    @param ignorebox (1D numpy array)
    The box that should be ignored because it contains the *actual* grid code
    zero. This hypercube starts at the origin and extends in the positive
    direction. It is applied in every high-dimensional quadrant.

    @param readoutResolution (float)
    The precision of readout of this grid code, measured in distance on the
    plane. For example, if this is 0.2, then all points on the plane are
    indistinguishable from those in their surrounding +- 0.1 range.

    @param pingInterval (float)
    How often, in seconds, the function should print its current status. If <=
    0, no printing will occur.

    @return
    - The largest tested scaling factor of the scaledbox that contains no
      collisions.
    - A point just outside this scaled scaledbox that collides with the origin.
    '''
    domainToPlaneByModule = np.asarray(
        domainToPlaneByModule, dtype='float64')
    latticeBasisByModule = np.asarray(
        latticeBasisByModule, dtype='float64')
    boxToScale = np.asarray(
        boxToScale, dtype='float64')
    ignoreBox = np.asarray(
        ignoreBox, dtype='float64')

    return _gridcodingrange.computeCodingRange(
        domainToPlaneByModule, latticeBasisByModule, boxToScale,
        ignoreBox, phaseResolution, pingInterval)


def computeGridUniquenessHypercube(domainToPlaneByModule, latticeBasisByModule,
                                   phaseResolution, ignoredCenterDiameter,
                                   pingInterval=10.0):
    '''
    Calls computeCodingRange with a unit cube scaledBox and cube ignore box.

    @param domainToPlaneByModule (3D numpy array)
    A list of 2*k matrices, one per module. The matrix converts from a point in
    the domain to a point on a plane, normalizing for grid cell scale.

    @param latticeBasisByModule (3D numpy array)
    A list of m 2*2 matrices, one per module. This matrix contains the basis
    vectors for a lattice, specifying which points on the plane have equivalent
    location representations in this module.

    @param ignoredCenterDiameter (float)
    The sidelength of the hypercube which should be ignored because it contains
    the *actual* grid code zero. This hypercube starts at the origin and extends
    in the positive direction. It is applied in every expansion direction (i.e.
    in every quadrant/octant/orthant).

    @param pingInterval
    How often, in seconds, the function should print its current status. If <=
    0, no printing will occur.

    @return
    - The diameter of the hypercube that contains no collisions.
    - A point just outside this hypercube that collides with the origin.
    '''
    domainToPlaneByModule = np.asarray(
        domainToPlaneByModule, dtype='float64')
    latticeBasisByModule = np.asarray(
        latticeBasisByModule, dtype='float64')

    return _gridcodingrange.computeGridUniquenessHypercube(
        domainToPlaneByModule, latticeBasisByModule, phaseResolution,
        ignoredCenterDiameter, pingInterval)


def computeBinSidelength(domainToPlaneByModule, phaseResolution,
                         resultPrecision, upperBound=1000.0, timeout=-1.0):
    '''
    Compute the sidelength of the smallest hypercube that encloses the
    intersection of all of the modules' firing fields centered at the origin.
    This sidelength is like a resolution of the code.

    @param domainToPlaneByModule
    A list of 2*k matrices, one per module. The matrix converts from a point in
    the domain to a point on a plane, normalizing for grid cell scale.

    @param readoutResolution
    The precision of readout of this grid code, measured in distance on the
    plane. For example, if this is 0.2, then all points on the plane are
    indistinguishable from those in their surrounding +- 0.1 range.

    @param resultPrecision
    The precision level for this sidelength. This algorithm will binary-search
    for the smallest hypercube it can place around zero, stopping at the point
    when it's within 'resultPrecision' of the actual smallest cube.

    @param upperBound
    Instructs the algorithm when it should give up. If the provided matrix
    never moves away from grid code zero, the algorithm could search forever.
    This parameter tells it when it should stop searching.

    @param timeout
    Specifies how long to try. This function will give up after 'timeout'
    seconds. If <= 0, the function will not time out. On timeout, the function
    throws an exception with message "timeout". In Python this exception is of
    type RuntimeError.

    @return
    The sidelength of this hypercube. Returns -1.0 if a surface can't be found
    (i.e. if upperBound is reached.)
    '''
    domainToPlaneByModule = np.asarray(
        domainToPlaneByModule, dtype='float64')

    return _gridcodingrange.computeBinSidelength(
        domainToPlaneByModule, phaseResolution, resultPrecision, upperBound, timeout)


def computeBinRectangle(domainToPlaneByModule, phaseResolution,
                        resultPrecision, upperBound=1000.0, timeout=-1.0):
    '''
    Like computeBinSidelength, but it computes a hyperrectangle rather than a
    hypercube.

    @param domainToPlaneByModule
    A list of 2*k matrices, one per module. The matrix converts from a point in
    the domain to a point on a plane, normalizing for grid cell scale.

    @param readoutResolution
    The precision of readout of this grid code, measured in distance on the
    plane. For example, if this is 0.2, then all points on the plane are
    indistinguishable from those in their surrounding +- 0.1 range.

    @param resultPrecision
    The precision level for this sidelength. This algorithm will binary-search
    for the smallest hyperrectangle it can place around zero, stopping at the
    point when each side is within 'resultPrecision' of the actual smallest
    hyperrectangle.

    @param upperBound
    Instructs the algorithm when it should give up. If the provided matrix
    never moves away from grid code zero, the algorithm could search forever.
    This parameter tells it when it should stop searching.

    @param timeout
    Specifies how long to try. This function will give up after 'timeout'
    seconds. If <= 0, the function will not time out. On timeout, the function
    throws an exception with message "timeout". In Python this exception is of
    type RuntimeError.

    @return
    The dimensions of this hyperrectangle. Returns an empty vector if a surface
    can't be found (i.e. if upperBound is reached.)
    '''
    domainToPlaneByModule = np.asarray(
        domainToPlaneByModule, dtype='float64')

    return _gridcodingrange.computeBinRectangle(
        domainToPlaneByModule, phaseResolution, resultPrecision, upperBound, timeout)
