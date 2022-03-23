#!/usr/bin/env python

import pytest
from scipy.optimize import RootResults

import itp
def test_success():
    # quadratic fucntion with boudaries 0 and 3
    f = lambda x: (x**2 - 1)

    root = itp.itp(f, 0, 3, full_output=True)
    assert root[0] == 1
    assert isinstance(root[1], RootResults)

def test_boundaries():
    # the root must be inside the boundaries
    f = lambda x: (x**2 - 1)

    with pytest.raises(ValueError):
        itp.itp(f, 2, 5)


def test_convergence_failure():
    # this really shouldn`t happen, but
    import math
    f = lambda x: (x*math.exp(x) -1)

    with pytest.raises(RuntimeError):
        itp.itp(f, -1, 2, maxiter=5)

