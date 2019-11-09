from sympy import *

import optics as op
"""
default parameters:
r, th, n, L, R, R1, R2, g, P
[NOTE]desn't enter op.main()
"""

P1 = Symbol('P1')
class GRIN_1(op.GRIN):
    """A GRIN lens with one end polished into convex shape"""

    def __init__(self):
        super().__init__(op.P)
        self._M = op.M7*op.M9*op.M1

    @property
    def Mtx(self):
        return self._M.subs(op.P, 0.25)

test = GRIN_1()
print("Pitch = 0.25")
op.Calculate(test, 'FFL')
