# ----------------------------------------------------------------------
# Numenta Platform for Intelligent Computing (NuPIC)
# Copyright (C) 2019, Numenta, Inc.  Unless you have an agreement
# with Numenta, Inc., for a separate license for this software code, the
# following terms and conditions apply:
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero Public License version 3 as
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero Public License for more details.
#
# You should have received a copy of the GNU Affero Public License
# along with this program.  If not, see http://www.gnu.org/licenses.
#
# http://numenta.org/licenses/
# ----------------------------------------------------------------------

import argparse
import math
import multiprocessing
import os
import pickle
import sys
import threading

import numpy as np
from scipy.stats import ortho_group

from gridcodingrange import computeBinSidelength, computeBinRectangle


def create_bases(k, s):
    assert(k>1)
    B = np.zeros((k,k))
    B[:2,:2] = np.array([
        [1.,0.],
        [0.,1.]
    ])
    for col in range(2,k):
        B[ :, col] = np.random.randn(k)
        B[ :, col] = B[ :, col]/np.linalg.norm(B[ :, col])

    Q = ortho_group.rvs(k)
    return np.dot(Q,s*B)


def create_orthogonal_projection_base(k, s):
    assert k > 1
    B = np.eye(k)
    Q = ortho_group.rvs(k)

    return np.dot(Q,s*B)


def random_point_on_circle():
    r = np.random.sample()*2.*np.pi
    return np.array([np.cos(r), np.sin(r)])


def create_params(m, k, orthogonal, normalizeScales=True):
    B = np.zeros((m,k,k))
    A = np.zeros((m,2,k))
    S = 1 + np.random.normal(size=m, scale=0.2)

    if normalizeScales:
        S /= np.mean(S)

    for m_ in range(m):
        if k==1:
            B[m_,0,0] = S[m_]*np.random.choice([-1.,1.])
            A[m_] = np.dot(random_point_on_circle().reshape((2,1)), np.linalg.inv(B[m_]))
        else:
            if orthogonal:
                B[m_] = create_orthogonal_projection_base(k, S[m_])
            else:
                B[m_] = create_bases(k, S[m_])
            A[m_] = np.linalg.inv(B[m_])[:2]

    return {
        "A": A,
        "S": S,
    }


def processCubeQuery1(query):
    phr, m, k, forceOrthogonal, normalizeScales, max_binsidelength = query

    upperBound = 4.0
    timeout = 60.0 * 10.0 # 10 minutes
    resultResolution = 0.01

    numDiscardedTooBig = 0
    numDiscardedTimeout = 0

    while True:
        expDict = create_params(m, int(math.ceil(k)), forceOrthogonal,
                                normalizeScales)

        # Rearrange to make the algorithms faster.
        sortOrder = np.argsort(expDict["S"])[::-1]
        expDict["S"] = expDict["S"][sortOrder]
        expDict["A"] = expDict["A"][sortOrder,:,:]

        try:
            binSidelength = computeBinSidelength(expDict["A"], phr,
                                                 resultResolution,
                                                 upperBound, timeout)

            if (binSidelength == -1.0 or
                (max_binsidelength is not None and
                 binSidelength >= max_binsidelength)):
                numDiscardedTooBig += 1
                continue

            expDict["bin_sidelength"] = binSidelength

            return (expDict, numDiscardedTooBig, numDiscardedTimeout)

        except RuntimeError as e:
            if e.message == "timeout":
                print("Timed out on query {}".format(expDict["A"]))

                numDiscardedTimeout += 1
                continue
            else:
                raise


class UniqueBasesScheduler(object):
    def __init__(self, folderpath, numTrials, ms, ks, phaseResolutions,
                 measureRectangle, allowOblique, normalizeScales, filtered):

        assert measureRectangle == False # Not supported (yet)

        self.folderpath = folderpath
        self.numTrials = numTrials
        self.ms = ms
        self.ks = ks
        self.phaseResolutions = phaseResolutions

        self.failureCounter = 0
        self.successCounter = 0

        self.pool = multiprocessing.Pool()
        self.finishedEvent = threading.Event()

        max_binsidelength = (1.0 if filtered else None)
        forceOrthogonal = not allowOblique
        self.param_combinations = [(phr, m, k, forceOrthogonal, normalizeScales,
                                    max_binsidelength)
                                   for phr in phaseResolutions
                                   for m in ms
                                   for k in ks
                                   if 2*m >= k]

        # Keep running tasks on all CPUs until we have generated enough results,
        # then kill the remaining workers.
        for _ in range(multiprocessing.cpu_count()):
            self.queueNewWorkItem()


    def join(self):
        try:
            if sys.version_info >= (3, 0):
                self.finishedEvent.wait()
            else:
                # Python 2
                # Interrupts (ctrl+c) have no effect without a timeout.
                self.finishedEvent.wait(9999999999)
            # Kill remaining workers.
            self.pool.terminate()
            self.pool.join()
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            self.pool.terminate()
            self.pool.join()


    def queueNewWorkItem(self):
        self.pool.map_async(processCubeQuery1, self.param_combinations,
                            callback=self.onWorkItemFinished)


    def onWorkItemFinished(self, results):
        if self.successCounter == self.numTrials:
            return

        discardedTooBig = np.full((len(self.phaseResolutions),
                                   len(self.ms),
                                   len(self.ks)),
                                  0, dtype="int")
        discardedTimeout = np.full((len(self.phaseResolutions),
                                    len(self.ms),
                                    len(self.ks)),
                                   0, dtype="int")
        binSidelengths = np.full((len(self.phaseResolutions),
                                  len(self.ms),
                                  len(self.ks)),
                                 np.nan, dtype="float")

        everyA = {}
        everyS = {}

        for params, result in zip(self.param_combinations, results):
            phr, m, k, _, _, _ = params
            expDict, numDiscardedTooBig, numDiscardedTimeout = result

            everyA[(phr, m, k)] = expDict["A"]
            everyS[(phr, m, k)] = expDict["S"]

            idx = (self.phaseResolutions.index(phr), self.ms.index(m),
                   self.ks.index(k))
            binSidelengths[idx] = expDict["bin_sidelength"]
            discardedTooBig[idx] += numDiscardedTooBig
            discardedTimeout[idx] += numDiscardedTimeout

        resultDict = {
            "phase_resolutions": self.phaseResolutions,
            "ms": self.ms,
            "ks": self.ks,
            "discarded_too_big": discardedTooBig,
            "discarded_timeout": discardedTimeout,
            "bin_sidelength": binSidelengths,
            "every_A": everyA,
            "every_S": everyS,
        }

        # Save the dict
        successFolder = os.path.join(self.folderpath, "in")
        if self.successCounter == 0:
            os.makedirs(successFolder)
        filepath = os.path.join(successFolder, "in_{}.p".format(
            self.successCounter))
        self.successCounter += 1
        with open(filepath, "wb") as fout:
            print("Saving {} ({} remaining)".format(
                filepath, self.numTrials - self.successCounter))
            pickle.dump(resultDict, fout)

        if self.successCounter == self.numTrials:
            self.finishedEvent.set()
        else:
            self.queueNewWorkItem()


def create_submatrixing_params(m, k, normalizeScales=True):
    # Create a 3D projection. Then for each module, append k-3 columns with
    # lengths U(0, maxlength) where maxlength is the longest of the first 3
    # columns.

    params = create_params(m, 3, orthogonal=True,
                           normalizeScales=normalizeScales)

    A = np.zeros((m,2,k))
    A[:,:,:3] = params["A"]

    for m_ in range(m):
        maxlength = max(np.linalg.norm(A[m_,:,k_])
                        for k_ in range(3))

        for k_ in range(3, k):
            length = np.random.uniform(0, maxlength)
            A[m_,:,k_] = random_point_on_circle()*length

    return {
        "A": A[:,:,:k],
        "S": params["S"],
    }


def getQuery(A, S, m, k, phase_resolution):
    A_ = A[:m, :, :k]
    sort_order = np.argsort(S[:m])[::-1]
    A_ = A_[sort_order, :, :]

    return (A_, phase_resolution)


def processCubeQuery2(query):
    A, phase_resolution = query
    resultResolution = 0.01
    upperBound = 2048.0
    timeout = 60.0 * 10.0 # 10 minutes

    try:
        result = computeBinSidelength(A, phase_resolution, resultResolution,
                                      upperBound, timeout)
        if result == -1.0:
            print("Couldn't find bin smaller than {} for query {}".format(
                upperBound, A.tolist()))

        return result
    except RuntimeError as e:
        if e.message == "timeout":
            print("Timed out on query {}".format(A.tolist()))
            return None
        else:
            raise


def processRectangleQuery(query):
    A, phase_resolution = query
    resultResolution = 0.01
    upperBound = 2048.0
    timeout = 60.0 * 10.0 # 10 minutes

    try:
        result = computeBinRectangle(A, phase_resolution, resultResolution,
                                     upperBound, timeout)

        if len(result) == 0:
            print("Couldn't find bin smaller than {} for query {}".format(
                upperBound, A.tolist()))

        return result
    except RuntimeError as e:
        if e.message == "timeout":
            print("Timed out on query {}".format(A.tolist()))
            return None
        else:
            raise



class IterableWithLen(object):
    def __init__(self, iterable, length):
        self.iterable = iterable
        self.length = length

    def __len__(self):
        return self.length

    def __iter__(self):
        return self.iterable



class ReuseBasesScheduler(object):
    def __init__(self, folderpath, numTrials, ms, ks, phaseResolutions,
                 measureRectangle, normalizeScales, buildupBases, filtered):
        self.folderpath = folderpath
        self.numTrials = numTrials
        self.ms = ms
        self.ks = ks
        self.phaseResolutions = phaseResolutions
        self.measureRectangle = measureRectangle
        self.normalizeScales = normalizeScales
        self.buildupBases = buildupBases

        self.failureCounter = 0
        self.successCounter = 0

        self.pool = multiprocessing.Pool()
        self.finishedEvent = threading.Event()

        self.param_combinations = [(phr, m, k)
                                   for phr in phaseResolutions
                                   for m in ms
                                   for k in ks
                                   if 2*m >= k]

        self.max_binsidelength = (1.0 if filtered else None)

        # Keep running tasks on all CPUs until we have generated enough results,
        # then kill the remaining workers.
        for _ in range(multiprocessing.cpu_count()):
            self.queueNewWorkItem()


    def join(self):
        try:
            # Interrupts (ctrl+c) have no effect without a timeout.
            if sys.version_info >= (3, 0):
                self.finishedEvent.wait()
            else:
                # Python 2
                # Interrupts (ctrl+c) have no effect without a timeout.
                self.finishedEvent.wait(9999999999)
            # Kill remaining workers.
            self.pool.terminate()
            self.pool.join()
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            self.pool.terminate()
            self.pool.join()


    def queueNewWorkItem(self):
        if self.buildupBases:
            resultDict = create_submatrixing_params(max(self.ms),
                                                    int(math.ceil(max(self.ks))),
                                                    self.normalizeScales)
        else:
            resultDict = create_params(max(self.ms),
                                       int(math.ceil(max(self.ks))),
                                       orthogonal=True,
                                       normalizeScales=self.normalizeScales)
        resultDict["phase_resolutions"] = self.phaseResolutions
        resultDict["ms"] = self.ms
        resultDict["ks"] = self.ks

        A = resultDict["A"]
        S = resultDict["S"]

        queries = (getQuery(A, S, m, int(math.ceil(k)), phr)
                   for phr, m, k in self.param_combinations)
        # map_async will convert this to a list if it can't get the length.
        queries = IterableWithLen(queries, len(self.param_combinations))

        if self.measureRectangle:
            operation = processRectangleQuery
        else:
            operation = processCubeQuery2

        context = ContextForSingleMatrix(self, resultDict, self.max_binsidelength)
        self.pool.map_async(operation, queries, callback=context.onFinished)


    def handleFailure(self, resultDict):
        failureFolder = os.path.join(self.folderpath, "failures")
        if self.failureCounter == 0:
            os.makedirs(failureFolder)

        filename = "failure_{}.p".format(self.failureCounter)
        self.failureCounter += 1

        filepath = os.path.join(failureFolder, filename)

        with open(filepath, "wb") as fout:
            print("Saving {} ({} remaining)".format(
                filepath, self.numTrials - self.successCounter))
            pickle.dump(resultDict, fout)

        self.queueNewWorkItem()


    def handleSuccess(self, resultDict, results):
        if self.successCounter == self.numTrials:
            return

        # Insert results into dict
        if self.measureRectangle:
            rectangles = {}
            for (phr, m, k), result in zip(self.param_combinations, results):
                rectangles[(phr, m, k)] = result

            resultDict["rectangles"] = rectangles
        else:
            bin_sidelengths = np.full((len(self.phaseResolutions),
                                       len(self.ms),
                                       len(self.ks)),
                                      np.nan, dtype="float")
            for (phr, m, k), result in zip(self.param_combinations, results):
                bin_sidelengths[self.phaseResolutions.index(phr), self.ms.index(m),
                                self.ks.index(k)] = result
            resultDict["bin_sidelength"] = bin_sidelengths

        # Save the dict
        successFolder = os.path.join(self.folderpath, "in")
        if self.successCounter == 0:
            os.makedirs(successFolder)
        filepath = os.path.join(successFolder, "in_{}.p".format(
            self.successCounter))
        self.successCounter += 1
        with open(filepath, "wb") as fout:
            print("Saving {} ({} remaining)".format(
                filepath, self.numTrials - self.successCounter))
            pickle.dump(resultDict, fout)

        if self.successCounter == self.numTrials:
            self.finishedEvent.set()
        else:
            self.queueNewWorkItem()


class ContextForSingleMatrix(object):
    def __init__(self, scheduler, resultDict, max_binsidelength):
        self.scheduler = scheduler
        self.resultDict = resultDict
        self.max_binsidelength = max_binsidelength

    def onFinished(self, results):
        failure = any(result is None or result == -1
                      for result in results)

        if self.max_binsidelength is not None:
            failure = failure or any(result >= self.max_binsidelength
                                     for result in results)

        if failure:
            self.scheduler.handleFailure(self.resultDict)
        else:
            self.scheduler.handleSuccess(self.resultDict, results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folderName", type=str)
    parser.add_argument("--numTrials", type=int, default=1)
    parser.add_argument("--m", type=int, required=True, nargs="+")
    parser.add_argument("--k", type=float, required=True, nargs="+")
    parser.add_argument("--phaseResolution", type=float, default=[0.2], nargs="+")
    parser.add_argument("--measureRectangle", action="store_true")
    parser.add_argument("--allowOblique", action="store_true")
    parser.add_argument("--normalizeScales", action="store_true")
    parser.add_argument("--filtered", action="store_true")
    parser.add_argument("--reuseBases", action="store_true")
    parser.add_argument("--buildupBases", action="store_true")

    args = parser.parse_args()

    cwd = os.path.dirname(os.path.realpath(__file__))
    folderpath = os.path.join(cwd, args.folderName)

    SchedulerClass = (ReuseBasesScheduler if args.reuseBases
                      else UniqueBasesScheduler)

    params = {
        "folderpath": folderpath,
        "numTrials": args.numTrials,
        "ms": args.m,
        "ks": args.k,
        "phaseResolutions": args.phaseResolution,
        "measureRectangle": args.measureRectangle,
        "normalizeScales": args.normalizeScales,
        "filtered": args.filtered,
    }

    if args.reuseBases:
        params["buildupBases"] = args.buildupBases
        ReuseBasesScheduler(**params).join()
    else:
        params["allowOblique"] = args.allowOblique
        UniqueBasesScheduler(**params).join()
