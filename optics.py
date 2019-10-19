from sympy import *
init_printing()

#Jones Matrix
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
#GRIN lens
#√A = g; √A z = p(pitch)
g, P = symbols('g P')
M9 = Matrix([[cos(P),sin(P)/g],[-g*sin(P),cos(P)]])

#light ray across a slab
#X2 = M3*M2*M1*X

#light ray across a general lens
#X2 = M7*M2*M5*X

#light ray across a c-lens
X2 = M7*M2*M1*X

#print(pretty(X2))
#pprint(X2)
#LaTex
#print(latex(X2))

#front focal length
#terminate to front surface
sol = solve(X2[1], th)
print('FFL:')
pprint(r/sol[0])

#back focal length
#originate from back surface
x2B = X2.subs(th,0)
print('BFL:')
pprint(-x2B[0]/x2B[1])
