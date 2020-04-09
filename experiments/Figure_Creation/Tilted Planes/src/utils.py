import numpy as np
from scipy.stats import ortho_group


hex_base = np.array([
    [1., np.cos(np.pi/3.)],
    [0., np.sin(np.pi/3.)]])


rombus = np.array([ [0., 0.],
                    [1., 0.],
                    [1. + np.cos(np.pi/3.), np.sin(np.pi/3.)],
                    [np.cos(np.pi/3.), np.sin(np.pi/3.)]])

def CAN_distance(p,q):
    q_ = (q - p)%1.
    q_ = q_.reshape((-1,2))
    q_ = np.dot(q_, hex_base.T).reshape((-1,1,2))
    d  = np.linalg.norm(q_ - rombus.reshape(1,4,2), axis=2)
    return np.amin(d,axis=1)

def lattice(orientation=0., angle=np.pi/3.):
    L = np.array([
        [np.cos(orientation), np.cos(orientation + angle)],
        [np.sin(orientation), np.sin(orientation + angle)],
    ])
    return L

def get_3dA(s):
    B = np.eye(3)
    B[:2,:2] = lattice(0,np.pi/3.)
    A = np.zeros((2,3))
    A = np.linalg.inv(s*B)[:2]
    return A


def map_to_torus(A, X):
    Phi = np.dot(X, A.T)%1.
    return Phi


def create_firing_field(X, A):
    Phi  = map_to_torus(A, X)
    Zero = np.zeros(Phi.shape)

    D = CAN_distance(Phi,Zero)
    pr = np.exp(-D**2/0.1)
    pr = pr/np.sum(pr)

    return pr



