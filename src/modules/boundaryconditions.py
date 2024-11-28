import numpy as np
from modules.mesh import *
from modules.conserved import *
from modules.thermodynamics import GAMMA

def bc_connect(W, mesh):
    """This function applies the connect boundary 
    condition (between jmin and jmax) to the conserved variables"""
    W.rho [:, 0] = W.rho [:, -(mesh.n_ghosts + 1)]
    W.rhou[:, 0] = W.rhou[:, -(mesh.n_ghosts + 1)]
    W.rhov[:, 0] = W.rhov[:, -(mesh.n_ghosts + 1)]
    W.rhoE[:, 0] = W.rhoE[:, -(mesh.n_ghosts + 1)]

    W.rho [:, 1] = W.rho [:, -(mesh.n_ghosts + 1)]
    W.rhou[:, 1] = W.rhou[:, -(mesh.n_ghosts + 1)]
    W.rhov[:, 1] = W.rhov[:, -(mesh.n_ghosts + 1)]
    W.rhoE[:, 1] = W.rhoE[:, -(mesh.n_ghosts + 1)]

    W.rho [:, -1] = W.rho [:, mesh.n_ghosts + 1]
    W.rhou[:, -1] = W.rhou[:, mesh.n_ghosts + 1]
    W.rhov[:, -1] = W.rhov[:, mesh.n_ghosts + 1]
    W.rhoE[:, -1] = W.rhoE[:, mesh.n_ghosts + 1]

    W.rho [:, -2] = W.rho [:, mesh.n_ghosts]
    W.rhou[:, -2] = W.rhou[:, mesh.n_ghosts]
    W.rhov[:, -2] = W.rhov[:, mesh.n_ghosts]
    W.rhoE[:, -2] = W.rhoE[:, mesh.n_ghosts]

    return W


def bc_wall(W, mesh):
    """This function applies the wall condition (no penetration) at cells 
    at i_phys = 0. It updates the ghost cells to have a no penetration condition at the wall"""
    i_wall = mesh.n_ghosts # Index i of the wall
    for j in range(mesh.n_ghosts, mesh.nj+mesh.n_ghosts):
        rho  = W.rho [i_wall][j]
        u    = W.rhou[i_wall][j]/rho
        v    = W.rhov[i_wall][j]/rho
        dj_x = mesh.Delta_jx[i_wall][j]
        dj_y = mesh.Delta_jy[i_wall][j]
        normal = np.array([dj_y, -dj_x])/np.sqrt(dj_x**2 + dj_y**2) # normal to constant i faces are oriented from i to i+1
        Vn   = u*normal[0] + v*normal[1]
        u    = u - 2*Vn*normal[0]
        v    = v - 2*Vn*normal[1]

        W.rho [0, j] = rho
        W.rhou[0, j] = rho*u
        W.rhov[0, j] = rho*v
        W.rhoE[0, j] = W.rhoE[i_wall, j]

        W.rho [1, j] = rho
        W.rhou[1, j] = rho*u
        W.rhov[1, j] = rho*v
        W.rhoE[1, j] = W.rhoE[i_wall, j]
    return W
    
def bc_farfield(W, mesh, M_inf, AOA):
    """This function implements the farfield boundary condition
    For the moment, only the supersonic case is implemented (M_inf > 1)
    """
    # Compute initial state (for inflow in supersonic cases)
    P_0 = 1                     # Stagnation pressure of free stream flow
    rho_0 = 1                   # Stagnation density of free stream flow
    c_0 = np.sqrt(P_0/rho_0)    # Velocity of sound 

    pressure_ratio = (1+0.5*(GAMMA-1)*M_inf**2)**(-GAMMA/(GAMMA-1))  # isentropic relation associated with Mach number for pressure
    rho_ratio = pressure_ratio**(1/GAMMA)                            # isentropic relation associated with Mach number for density
                                 
    P = pressure_ratio*P_0
    rho = rho_ratio*rho_0
    c = np.sqrt(GAMMA*P/rho)/c_0
    u = c*M_inf*np.cos(AOA*(np.pi)/180)
    v = c*M_inf*np.sin(AOA*(np.pi)/180)
    rhoE = (1/(GAMMA-1))*P+0.5*rho*(u**2+v**2)

    # Get farfield (last physical cell at i=imax) state
    rho_farfield = W.rho [-3, :]
    u_farfield   = W.rhou[-3, :]/rho_farfield
    v_farfield   = W.rhov[-3, :]/rho_farfield

    # Verify if inflow or outflow (sign(V.n))
    dj_x = mesh.Delta_jx[-3][:]
    dj_y = mesh.Delta_jy[-3][:]

    Vn     = u_farfield*dj_y - v_farfield*dj_x   # normal is not unit length but ok since we look at the sign
    for j in range(mesh.nj + 2*mesh.n_ghosts):
        if Vn[j] >=0 : # inflow
            W.rho [-2, j]  = rho
            W.rhou[-2, j] = rho*u
            W.rhov[-2, j] = rho*v
            W.rhoE[-2, j] = rhoE
        else: # outflow
            W.rho [-2, j] = W.rho [-3, j]
            W.rhou[-2, j] = W.rhou[-3, j]
            W.rhov[-2, j] = W.rhov[-3, j]
            W.rhoE[-2, j] = W.rhoE[-3, j]
    W.rho [-1, :] = W.rho [-2, :]
    W.rhou[-1, :] = W.rhou[-2, :]
    W.rhov[-1, :] = W.rhov[-2, :]
    W.rhoE[-1, :] = W.rhoE[-2, :]
    return W

def apply_boundary_conditions(W, mesh, M_inf, AOA):
    W = bc_connect (W, mesh)
    W = bc_wall    (W, mesh)
    W = bc_farfield(W, mesh, M_inf, AOA)
    return W