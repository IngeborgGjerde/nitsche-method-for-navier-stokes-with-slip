
import sympy as sym
import numpy as np
from fenics import *
import math


def potential_flow(radius, center, mesh):
    '''
    Set up potential flow solution for flow past a cylinder
    Pressure is rescaled (using mesh) so that int p dx = 0

    Args: 
        radius (float): cylinder radius
        center (list): [x,y] coords of cylinder center
        mesh (fenis mesh): domain mesh (for computing the rescaling)

    Returns:
        p_a: potential flow pressure
        u_a: potential flow velocity
        f: Right-hand size, i.e. just the zero-vector
        up: fenics tuple with the entire solution
    '''
    
    
    import sympy as sym
    a = radius
    x, y, z = sym.symbols('x[0] x[1] x[2]')

    r2 = (x-center[0])**2.0 + (y-center[1])**2.0

    phi = a*x + x/r2

    u0 = sym.diff(phi, x); u1 = sym.diff(phi, y); u2 = sym.diff(phi, z)

    norm_u2 = u0**2.0 + u1**2.0 + u2**2.0
    p = - 0.5*norm_u2 
    
    p_a = Expression(sym.printing.ccode(p), degree=2)
    u_a = Expression((sym.printing.ccode(u0),sym.printing.ccode(u1)), degree=2)
    
    # Rescale pressure by a constant so that p*dx = 0 over domain
    Q = FunctionSpace(mesh, 'CG', 1)
    pdx = assemble(interpolate(p_a, Q)*dx)
    size_domain = assemble(interpolate(Constant(1.0), Q)*dx)
    p -= pdx/size_domain
    p_a = Expression(sym.printing.ccode(p), degree=2)

    f = Constant((0.0, 0.0))
    
    up = Expression((sym.printing.ccode(u0),sym.printing.ccode(u1),sym.printing.ccode(p), '0.0'), degree=2)
    
    return p_a, u_a, f, up


def polynomial_flow(nu, beta, mesh):
    '''
    Set up tangential flow past a cylinder, (so a linear function),
    and compute the corresponding rhs

    Args: 
        radius (float): cylinder radius
        mesh (fenis mesh): domain mesh (for computing the rescaling)

    Returns:
        p_a: potential flow pressure
        u_a: potential flow velocity
        f: Right-hand size, i.e. just the zero-vector

    '''
    
    import sympy as sym
    x, y = sym.symbols('x[0] x[1]')

    u0, u1 = -y, x # then div u = 0
    p = 0.5*(x + y)

    # manufacture forcing term on gamma, f_gamma = beta dot(u, tau) + nu*n^t*D(u)*tau
    fgamma = beta*u0*y + beta*u1*(-x) # beta dot(u,tau)

    # manufacture forcing term for momentum equation
    f0 = sym.diff(p, x) - nu * sym.diff(sym.diff(u0, x), x) - nu * sym.diff(sym.diff(u0, y), y)
    f0 += u0 * sym.diff(u0, x) + u1 * sym.diff(u0, y)

    f1 = sym.diff(p, y) - nu * sym.diff(sym.diff(u1, x), x) - nu * sym.diff(sym.diff(u1, y), y)
    f1 += u0 * sym.diff(u1, x) + u1 * sym.diff(u1, y)

    # Make fenics expressions and return
    p_a = Expression(sym.printing.ccode(p), degree=2)
    u_a = Expression((sym.printing.ccode(u0),sym.printing.ccode(u1)), degree=2)
    f = Expression((sym.printing.ccode(f0),sym.printing.ccode(f1)), degree=2)
    f_gamma = Expression(sym.printing.ccode(fgamma), degree=2)

    up_a = Expression((sym.printing.ccode(u0), sym.printing.ccode(u1), sym.printing.ccode(p), '0.0'), degree=2)

    # Rescale pressure by a constant so that p*dx = 0 over domain
    Q = FunctionSpace(mesh, 'CG', 1)
    pdx = assemble(interpolate(p_a, Q) * dx)
    size_domain = assemble(interpolate(Constant(1.0), Q) * dx)
    p -= pdx / size_domain

    p_a = Expression(sym.printing.ccode(p), degree=2)


    return p_a, u_a, f, f_gamma, up_a