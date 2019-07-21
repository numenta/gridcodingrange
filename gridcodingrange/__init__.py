import _gridcodingrange

import numpy as np


def computeGridUniquenessHypercube(domainToPlaneByModule, latticeBasisByModule,
                                   phaseResolution, ignoredCenterDiameter,
                                   pingInterval=10.0):
    domainToPlaneByModule = np.asarray(
        domainToPlaneByModule, dtype="float64")
    latticeBasisByModule = np.asarray(
        latticeBasisByModule, dtype="float64")

    return _gridcodingrange.computeGridUniquenessHypercube(
        domainToPlaneByModule, latticeBasisByModule, phaseResolution,
        ignoredCenterDiameter, pingInterval)


def computeBinSidelength(domainToPlaneByModule, phaseResolution,
                         resultPrecision, upperBound=1000.0, timeout=-1.0):
    domainToPlaneByModule = np.asarray(
        domainToPlaneByModule, dtype="float64")

    return _gridcodingrange.computeBinSidelength(
        domainToPlaneByModule, phaseResolution, resultPrecision, upperBound, timeout)


def computeBinRectangle(domainToPlaneByModule, phaseResolution,
                        resultPrecision, upperBound=1000.0, timeout=-1.0):
    domainToPlaneByModule = np.asarray(
        domainToPlaneByModule, dtype="float64")

    return _gridcodingrange.computeBinRectangle(
        domainToPlaneByModule, phaseResolution, resultPrecision, upperBound, timeout)
