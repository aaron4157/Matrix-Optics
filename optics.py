from sympy import *
init_printing()
#abstract base class
from abc import ABCMeta, abstractmethod

#Ray-transfer vector
r, th = symbols('r θ')
X = Matrix([r,th])

#Prpagation
L = Symbol('L')
M2 = Matrix([[1,L],[0,1]])
#plane refraction(denser)
n = Symbol('n')
M1 = Matrix([[1,0],[0,1/n]])
#plane refraction(thiner)
M3 = Matrix([[1,0],[0,n]])
#into convex lens surface
R1 = Symbol('R1')
M5 = Matrix([[1,0],[(1-n)/n/R1,1/n]])
#out convex lens surface
#note R reverse here
R2 = Symbol('R2')
M7 = Matrix([[1,0],[(1-n)/R2,n]])
#reflection
#concave curverture radius
#convex: R->-R
R = Symbol('R')
M6 = Matrix([[1,0],[2/R,1]])
#GRIN media
#√A = g; √A z = P(pitch)
#To make lens, P ~ 1/4
g, P = symbols('g P')
M9 = Matrix([[cos(2*pi*P),sin(2*pi*P)/g],[-g*sin(2*pi*P),cos(2*pi*P)]])


#light ray across a general lens
#X2 = M7*M2*M5*X

#light ray across a c-lens
X2 = M7*M2*M1*X

#print(pretty(X2))
#pprint(X2)
#LaTex
#print(latex(X2))


class Component(object, metaclass=ABCMeta):
    """A general optical component, including lens, mirror, prism, shutter,
    is mostly specified with refractive index and midension.
    transfer matrix could be brought into that of some complex optical linear system.
    """

    def __init__(self, L):
        self._thickness = L

    @property
    def thickness(self):
        return self._thickness

    @property
    def Mtx(self):
        pass

    @abstractmethod
    def reasonable(self):
        pass

class Slab(Component):
    def __init__(self, ns, ds):
        super().__init__(L)
        self._M = M3*M2*M1
        self._index = ns
        self._thickness = ds

    @property
    def index(self):
        return self._index

    @property
    def Mtx(self):
        return self._M.subs([(n, self._index),(L, self._thickness)])

    def reasonable(self):
        try:
            if self._thickness <= 0 or self._index < 1:
                return False
            elif self._index == 1:
                return 'Void.'
            else:
                return True

        except TypeError as e:
            print('Not numerically assigned yet.')
            return True

class lens(Component):
    """A general thick lens. 
    parameters: refractive index nl,
    front surface curvature Ra,
    back surface curvature Rb
    thickness L
    """

    def __init__(self, nl, Ra, Rb, L):
        super().__init__(L)
        self._M = M7*M2*M5
        self._index = nl
        self._R1 = Ra
        self._R2 = Rb

    @property
    def Mtx(self):
        return self._M.subs([(n, self._index),(R1, self._R1),(R2, self._R2),(L, self._thickness)])

    def reasonable(self):
        try:
            if self._thickness <= 0 or self._index < 1:
                return False
            elif self._index == 1:
                return 'Void.'
            else:
                return True

        except TypeError as e:
            print('Not numerically assigned yet.')
            return True

class mirror(Component):
    """A reflection mirror
    parameter: curvature (> 0 for collimating reflective surface)
    """

    def __init__(self, Rc):
        super().__init__(L)
        self._M = M6
        self._curvature = Rc

    @property
    def Mtx(self):
        return self._M.subs(R,self._curvature)

    def  reasonable():
        try:
            if self._curvature == 0:
                return False
            else:
                return True
        except TypeError as e:
            print("Not numerically assigned...")
            return True

class GRIN(Component):
    """A GRIN media. pitch = 0.25 for collimating lens.
    parameters: pitch
    """

    def __init__(self, P1):
        super().__init__(L)
        self._M = M3*M9*M1
        self._P = P1
    @property
    def Mtx(self):
        return self._M.subs(P, self._P)

    @property
    def Pitch(self):
        return self._P

    def reasonable(self):
        try:
            return True if self._P > 0 else False
        except Exception as e:
            print('Not numerically assigned yet.')
            return True

def Calculate(Component, cmd):
    r_e = cmd.lower()
    output = Component.Mtx*X
    collimated = solve(output[1], th)
    focusing = output.subs(th,0)
    imaging = output.subs(r,0)

    if r_e == 'shift':
        print("Shift:")
        pprint((Component.Mtx*X - X)[0])

    elif r_e == 'optical thickness':
        #obtain negative value, if rays are focused outside
        print("image depth:")
        pprint(imaging[0]/imaging[1])

    elif r_e == 'ffl':
        #terminate to front surface
        print("Front focal length:")
        pprint(r/collimated[0])

    elif r_e == 'bfl':
        #originate from back surface
        print("Back focal length:")
        pprint(-focusing[0]/focusing[1])

    elif r_e == 'na':
        print("Radius = a.")
        a = Symbol('a')
        print("Numerical aperture:")
        pprint(collimated[0].subs(r, a))

    elif r_e == 'output':
        print("Output vector:")
        pprint(output)

    else:
        print("Request not supported...")


def main():
    n1, t1 = symbols('n1 t1')
    prism1 = Slab(n1, t1)
    #print(prism1.reasonable())
    print("Prism slab...")

    pprint(prism1.Mtx)
    Calculate(prism1,'Optical thickness')

    prism2 = Slab(2.0, 0.5)

    class combination(Component):
        """Describes two attached slabs, prism1 & prism2.
        param: total physical thickness
        """

        def __init__(self, L):
            super().__init__(L)
            self._M = prism2.Mtx*M2*prism1.Mtx

        @property
        def Mtx(self):
            return self._M

        def reasonable(arg):
            pass

    sets0 = combination('t1 + 0.5')
    Calculate(sets0, 'Shift')

    print("C-lens...")
    R = Symbol('R')
    lens1 =  lens(n1, oo, R, t1)
    Calculate(lens1, 'BFL')

if __name__ == '__main__':
    main()
