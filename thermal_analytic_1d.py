import scipy as sp
import numpy as np
import scipy.linalg as la
import sympy as sym



k_si = 130
k_ox = 1.38
k_cu = 385

alpha_si = k_si/2329/700
alpha_ox = k_ox/2203/703
alpha_cu = 1.11e-4

h_air = 1.8e4
h_water = 4.6e4
h_package = 5

alpha = np.array( [alpha_si, alpha_ox, alpha_si] )
k_actual = np.array( [k_si, k_ox, k_si] )
thickness_actual = np.array( [50, 5, 50] ) * 1e-6
pdens_cm2 = np.array( [0, 100, 0] )
h_actual = np.array([h_air, h_water] )

dTR = 1


x_actual = np.cumsum(thickness_actual)

dx = np.diff(x_actual)
dx = np.insert(dx, x_actual[0], 0)

# Construct useful constants
pdens_m2 = pdens_cm2 * 1e4
pdens_m3 = pdens_m2 / dx
hM = h_actual[-1] * x_actual[0]/k_actual[-1]
h1 = h_actual[0] * x_actual[0]/k_actual[0]

num_layers = len(alpha)

# Construct important vectors
k_vec = k_actual[1:]/k_actual[0:-1]
eta_vec = x_actual/x_actual[0]
b_vec = alpha / alpha[0]
g_vec = b_vec * x_actual[0]**2 * pdens_m3 / k_actual / dTR



# Construct basis functions
[a, x, b] = sym.symbols(['a','x','b'])

# Static symbolic basis functions
phi_s = sym.cos(a*x/sym.sqrt(b))
psi_s = sym.sin(a*x/sym.sqrt(b))
dphi_s = sym.diff(phi_s, x)
dpsi_s = sym.diff(psi_s, x)

# Actual basis functions
phi = lambda xx, aa, bb: phi_s.subs( [(a, aa), (b, bb), (x, xx)] )
psi = lambda xx, aa, bb: psi_s.subs( [(a, aa), (b, bb), (x, xx)] )
dphi = lambda xx, aa, bb: dphi_s.subs( [(a, aa), (b, bb), (x, xx)] )
dpsi = lambda xx, aa, bb: dpsi_s.subs( [(a, aa), (b, bb), (x, xx)] )

N1 = np.array([ -dphi(0, a, b_vec[0])/h1 + phi(0, a, b_vec[0]), -dpsi(0, a, b_vec[0])/h1 + psi(0, a, b_vec[0]) ])
N1 = np.append( N1, np.zeros( (1, 2*(num_layers-1) ) )  )

P_core = lambda n, eta, a, b_vec: np.array([phi(eta[n-1], a, b_vec[n-1]), psi(eta[n-1], a, b_vec[n-1]), -phi( eta[n], a, b_vec[n-1] ), -psi( eta[n], a, b_vec[n-1] )])
P = lambda n, eta, a, b_vec: np.append( np.append( np.zeros((1, 2*(n-1))), P_core(n, eta, a, b_vec)), np.zeros(( 1, 2*(num_layers - n - 1)) ) )

print(np.size(N1))




