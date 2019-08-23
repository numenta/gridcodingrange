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
import os
import pickle
import re

import numpy as np

from gridcodingrange import computeCodingRange


def create_L(m, theta=np.pi/3.):
    L = np.zeros((m,2,2))
    for i in range(m):
      L[i] = np.array([
          [np.cos(0.), np.cos(theta)],
          [np.sin(0.), np.sin(theta)]
        ])
    return L


def extractNumber(filename):
    return int(re.match("in_([0-9]+).p", filename).groups()[0])


def measureSidelengths(folderPath, scaleMinimalBox):
    inPath = os.path.join(folderPath, "in")
    outPath = os.path.join(folderPath, "out")

    if not os.path.exists(outPath):
        os.makedirs(outPath)

    files = sorted((filename
                    for filename in os.listdir(inPath)
                    if re.match("in_[0-9]+.p", filename)),
                   key=extractNumber)

    for filename in files:
        inFilePath = os.path.join(folderPath, "in", filename)

        iFile = int(re.match("in_([0-9]+).p", filename).groups()[0])

        print("Computing result {}".format(iFile))

        outFilePath = os.path.join(folderPath, "out", "res_{}.p".format(iFile))

        with open(inFilePath, "rb") as fin:
            try:
                result_dict = pickle.load(fin)
            except UnicodeDecodeError:
                # Probably loading Python2 pickled file from Python3
                result_dict = pickle.load(fin, encoding='latin1')

        ms = result_dict["ms"]
        ks = result_dict["ks"]
        phase_resolutions = result_dict["phase_resolutions"]
        rectangles = result_dict["rectangles"]

        L = create_L(max(ms))

        param_combinations = [(phr, m, k)
                              for phr in phase_resolutions
                              for m in ms
                              for k in ks
                              if 2*m >= k]

        max_scale_factors = np.full((len(phase_resolutions),
                                     len(ms), len(ks)),
                                    np.nan, dtype="float")

        for phr, m, k in param_combinations:

            if "A" in result_dict:
                A_ = result_dict["A"][:m, :, :int(math.ceil(k))]
                S = result_dict["S"]
                sort_order = np.argsort(S[:m])[::-1]
                A_ = A_[sort_order, :, :]
            else:
                A_ = result_dict["every_A"][(phr, m, k)]

            L_ = L[:m]

            ignorebox = 0.51*rectangles[(phr, m, k)]
            if scaleMinimalBox:
                scaledbox = ignorebox
            else:
                scaledbox = np.ones(int(math.ceil(k)), dtype='float')

            partial_final_dim = k - math.floor(k)
            if partial_final_dim > 0:
                scaledbox[-1] = partial_final_dim

            max_scale_factor, _ = computeCodingRange(
                A_, L_, scaledbox, ignorebox, phr,
                pingInterval=100.0)

            max_scale_factors[phase_resolutions.index(phr),
                              ms.index(m), ks.index(k)] = max_scale_factor

        result_dict["width"] = max_scale_factors

        with open(outFilePath, "wb") as fout:
            print("Saving {}".format(outFilePath))
            pickle.dump(result_dict, fout)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("folderName", type=str)
    parser.add_argument("--scaleMinimalBox", action="store_true")

    args = parser.parse_args()

    cwd = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(cwd, args.folderName)

    measureSidelengths(folderPath, args.scaleMinimalBox)
