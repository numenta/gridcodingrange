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

import unittest

import math
import pickle

import numpy as np
from scipy.stats import ortho_group

from gridcodingrange import computeCodingRange, computeBinSidelength

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


def create_params(m, k, orthogonal):
    B = np.zeros((m,k,k))
    A = np.zeros((m,2,k))
    S = np.ones(m) + np.random.sample(m)

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

def create_L(m, theta=np.pi/3.):
    L = np.zeros((m,2,2))
    for i in range(m):
        L[i] = np.array([
            [np.cos(0.), np.cos(theta)],
            [np.sin(0.), np.sin(theta)]
        ])
    return L


class AlgorithmTweaksTest(unittest.TestCase):
    def testExpandBoxAsSubspace(self):
        for i in range(100):
            m = 4
            max_k = 4
            A = create_params(m, max_k, True)['A']
            phr = 0.2
            ignorebox_width = 0.51*computeBinSidelength(A, 0.2, 0.01, 1000)
            L = create_L(m)

            scaledbox = np.array([1.0, 1.0], dtype='float')
            ignorebox = ignorebox_width*np.ones(2, dtype='float')
            A_ = A[:,:,:2]
            L_ = L[:,:,:2]
            result1 = computeCodingRange(A_, L_, scaledbox, ignorebox, phr)

            scaledbox = np.array([1.0, 1.0, 0.0, 0.0], dtype='float')
            ignorebox = ignorebox_width*np.ones(4, dtype='float')
            result2 = computeCodingRange(A, L, scaledbox, ignorebox, phr)

            self.assertEqual(
                result1[0],
                result2[0],
                "Different results for A: {} L: {}, results {} != {}".format(
                    A.tolist(),
                    L.tolist(),
                    result1,
                    result2))


if __name__ == "__main__":
  unittest.main()
