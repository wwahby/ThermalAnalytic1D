import scipy as sp
import numpy as np
import scipy.linalg as la
import sympy as sym




[a, x, b] = sym.symbols(['a','x','b'])

phi = sym.cos(a*x/sym.sqrt(b))
psi = sym.sin(a*x/sym.sqrt(b))

dphi = sym.diff(phi, x)
dpsi = sym.diff(psi, x)

N1 =  [ -dphi