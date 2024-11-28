from os import makedirs
from time import time
import argparse

from modules.mesh import *
from modules.conserved import *
from modules.boundaryconditions import *
from modules.monitor import *
from modules.writer import *
import modules.numerics as Sol
import modules.post as POST
import modules.utils as UTILS

def main(mesh=None, iterations=None, AOAdeg=None, M_inf=None, CFL=None):
    startTime = time()
    # Mesh file 
    # Read imput parameters
    if mesh == None:
        mesh_filepath = "NACA0012grids/33x33.x"
    else:
        mesh_filepath = mesh

    if iterations == None:
        itermax = 100000
    else:
        itermax = int(iterations)
        
    if AOAdeg == None:
        AOAdeg = 0.0
    else:
        AOAdeg = float(AOAdeg)
    AOA = AOAdeg*np.pi/180

    if M_inf == None:
        AOAdeg = 1.5
    else:
        M_inf = float(M_inf)
    
    if CFL == None:
        CFL = 1.0
    else:
        CFL = float(CFL)

    # Numerical parameters
    n_ghosts = 2
    
    # Beginning of the computation
    time_start = time()

    # Construct mesh
    mesh = Mesh(filepath=mesh_filepath, n_ghosts=n_ghosts)

    # Output paths -> defines where to store the results
    folder = "output/"+ str(mesh.ni) + "_" + str(AOAdeg) + "_" + str(M_inf)+ "_" + str(CFL) + "/"

    # Create the folder if it does not exist
    makedirs(folder, exist_ok=True)

    # Filenames
    output_path_results = folder + "Flow.dat"
    output_path_res_fig = folder + "Residuals.png"
    output_path_res_txt = folder + "Residuals.txt"
    output_path_Cp_fig  = folder + "Cp.png"
    output_path_Cp_txt  = folder + "Cp.txt"
    output_path_CL_CD   = folder + "CL_CD.txt"
    output_comp_time    = folder + "computational_time.txt"

    # Construct conserved variables and init the flow field
    W = ConservedVariables(ni=mesh.ni, nj=mesh.nj, n_ghosts=n_ghosts)
    W.init_flow(M_inf, AOA)

    # Apply boundary conditions
    W = apply_boundary_conditions(W, mesh, M_inf, AOA)

    centralscheme = Sol.CentralSchemeWithArtificialDissipation(mesh, k2=1/2, k4=1/64)
    RK2 = Sol.RK2(CFL, mesh, local_time=True)
    monitor = Monitor(centralscheme, RK2, itermax=itermax, eps=10**(-10))

    # Iteration loop
    monitor.iterate(mesh, W, M_inf, AOA)

    # Get the computational time
    time_end = time()
    computational_time = time_end-time_start

    monitor.plot_residuals(output_path_res_fig, output_path_res_txt)

    # Write the results file
    write4Tecplot(output_path_results, mesh, W)

    # Post-processing
    POST.plot_save_CP(W, mesh, M_inf, AOA, output_path_Cp_fig, output_path_Cp_txt)
    POST.compute_save_CL_CD(M_inf, AOA, mesh, W, output_path_CL_CD)
    POST. write_computational_time(computational_time, output_comp_time)

    endTime = time()
    print("STATS:")
    print("-"*6)
    print(f"Program took {(endTime - startTime):10.03f}s to execute.")


if __name__ == "__main__":

    # Initialize the parser
    parser = argparse.ArgumentParser(description="PyEuler 2D")
    parser.add_argument('-i', '--iterations', metavar="", help="Number of iterations.")
    parser.add_argument('-C', '--CFL', metavar="", help="Number of iterations.")
    parser.add_argument('-m', '--mesh', metavar="", help="Pass the location of the mesh file after this flag.")
    parser.add_argument('-A', '--AOA', metavar="", help="Angle of Attack.")
    parser.add_argument('-M', '--MACH', metavar="", help="Mach number.")
    arguments = parser.parse_args()

    main(mesh=arguments.mesh,
        iterations=arguments.iterations,
        AOAdeg=arguments.AOA,
        M_inf=arguments.MACH,
        CFL = arguments.CFL
        )