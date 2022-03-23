#!/usr/bin/env python
""" Implementation of the ITP method for finding roots.

Implementation of the method proposed on the paper:

I. F. D. Oliveira and R. H. C. Takahashi. 2019. An enhancement 
of the bisection method average performance
preserving minmax optimality. ACM Trans. Math. Softw. 0, 0, 
Article 0 ( 2019), 25 pages.
"""

__author__ = "Fernanda Takahashi"


import math
import numpy as np
from scipy.optimize import RootResults


MAX_ITER = 100
TOLERANCE = 0.0000000001 # 10e-10
KAPPA_1 = 0.1
KAPPA_2 = 2.0

def itp(f, a, b, xtol=TOLERANCE, k1=KAPPA_1, k2=KAPPA_2, 
        nrelax=0, maxiter=MAX_ITER, full_output=False,
        *args, **kargs):
    """
    Find root of a function within an interval using ITP.
    Interpolation, Truncation and Projection method that transforms
    the traditional Regular-Falsi into a minmax strategy with superlinear
    convergence.

    Parameters
    ----------
    f : function
        Python function returning a number.
    a : scalar
        One end of the bracketing interval [a,b].
    b : scalar
        The other end of the bracketing interval [a,b].
    xtol : number, optional
		Tolerance
    k1 : number, optional
		Multiplicative factor of the convergence rate
    k2 : number, optional
		Exponetail factor if the convergence rate
    args : arguments
        Extra arguments for the f function
    kargs: arguments
        Extra named arguments for the f function
    Returns
    -------
    x : float
        Zero of `f` between `a` and `b`.

    """

    if not isinstance(args, tuple):
        args = (*args,)
    if not isinstance(kargs, dict):
        kargs = None

    if xtol <= 0:
        raise ValueError("xtol too small (%g <= 0)" % xtol)

    return _itp(f, a, b, xtol, k1, k2, nrelax, maxiter, full_output, *args, **kargs)


def _itp(f, a, b, xtol, k1, k2, nrelax, maxiter, full_output, *args, **kargs):
    fa = f(a, *args, **kargs)
    fb = f(b, *args, **kargs)
    if fa*fb >0 :
        raise ValueError("f(a) and f(b) must have different signs")
    n_bisect = math.ceil(math.log((b-a)/(2*xtol),2))
    n_max = n_bisect + nrelax
    for k in range(maxiter):
        x = x_bisect = (a+b)/2
        if (b-a) <= (2*xtol):
            return (x, RootResults(x, k, k+2, 0)) if full_output else x
        # Interpolation
        x_f = (b*fa-a*fb)/(fa-fb) #regula falsi

        # Truncation
        sigma = np.sign(x_bisect - x_f)
        delta = k1*(abs(b-a))**k2
        x_trunc = x_f + sigma*delta if delta <= abs(x_bisect - x_f) else x_bisect

        # Projection  
        r = xtol*(2**(n_max-k)) - (b-a)/2
        x = x_trunc if abs(x_trunc - x_bisect) < r else x_bisect - sigma*r

        fx = f(x, *args, **kargs)
        if fx*fa > 0:
            a = x
            fa = fx
        elif fx*fb > 0:
            b = x
            fb = fx
        else:
            a = b = x
    raise RuntimeError("Just how?")
