## Run convergence test for Nitsches method for approximation potential flow past a cylinder

from fenics import *
import solve as solvers
import testproblems as testproblems
import utils as utils
from mesh_functions import get_meshes, mesh_size_surface
parameters["refinement_algorithm"] = "plaza_with_parent_facets" 
'''
Script for running the convergence tests in the paper
'''

def run_convergence_test(N_i, params):

    from pathlib import Path
    code_path = str(Path.cwd())

    ## Get parameters
    beta, nu, gamma_i, degree = [params[name] for name in ['beta', 'nu', 'gammas', 'degree']]
    gammas = [10 ** i for i in range(-gamma_i, gamma_i+1)]

    ## Get meshes
    Ns = [2 ** (i+2) for i in range(0, N_i)]

    radius, center = 1.0, [0.0, 0.0]
    #meshes, facet_functions, hs, hgammas = get_meshes(Ns, use_refined_mesh=params['cyl_refinement'])
    meshes = [Mesh()]
    with XDMFFile("mesh.xdmf") as xdmf:
        xdmf.read(meshes[0])
        mvc = MeshValueCollection("size_t", meshes[0], 1) 
    with XDMFFile("facet_mesh.xdmf") as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(meshes[0], mvc)
    facet_functions = [mf]
    hs = [meshes[0].hmax()]
    hgammas = [mesh_size_surface(facet_functions[0], 1)]

    for i in range(1, N_i):
        new_mesh = adapt(meshes[-1])
        new_mf = adapt(facet_functions[-1], new_mesh)
        mesh_size_surface
        meshes.append(new_mesh)
        facet_functions.append(new_mf)
        hs.append(meshes[-1].hmax())
        hgammas.append(mesh_size_surface(facet_functions[-1], 1))

    ## Get testproblem
    if beta == -2.0:
        p_a, u_a, f, up_a = testproblems.potential_flow(radius, center, meshes[-1])
        f_gamma = Constant(0.0)
    else:
        p_a, u_a, f, f_gamma, up_a = testproblems.polynomial_flow(nu, beta, meshes[-1])

    u_errors, p_errors, a_errors = [], [], []


    ## Run convergence test
    print('*******************************************************')
    print('h     h_Gamma    ||u_e||_a   ||u_e||_H1    ||p_e||_L2 ')
    print('*******************************************************')
    for ix in range(len(meshes)):

        # Get mesh and boundary markers
        mesh, boundary_markers = meshes[ix], facet_functions[ix]

        V2 = VectorFunctionSpace(mesh, 'CG', degree+1, 2)
        Q2 = FunctionSpace(mesh, 'CG', degree)

        u_a_i, p_a_i = interpolate(u_a, V2), interpolate(p_a, Q2)

        # Solve navier stokes with given parameters and record errors
        u_error, p_error, a_error = [], [], []
        for ix_gamma, gamma in enumerate(gammas):
            
            u, p, lamda, _ = solvers.navier_stokes(mesh, boundary_markers, params, gamma, u_a, up_a, f, f_gamma)

            u_error.append( errornorm(u, u_a_i, 'H1') )
            p_error.append( errornorm(p, p_a_i) )

            a_e = solvers.get_a_norm_error(u_a_i, u, V2, boundary_markers, params, gamma)
            a_error.append( a_e )

        u_errors.append(u_error); p_errors.append(p_error); a_errors.append(a_error)

        # we print out the conv. rates for the largest gamma for the user to see
        print(f'{hs[ix]:1.2f}  {hgammas[ix]:1.2f}       {a_error[-1]:1.7f}    {u_error[-1]:1.7f}    {p_error[-1]:1.7f}')


    ## Write errors and convergence rates to file
    errors = [p_errors, u_errors, a_errors]

    file_loc = code_path + '/Results/Tables/'
    ffile_names = [file_loc +'p_error', file_loc + 'u_error', file_loc + 'a_error']
    variable_names = ['$\Vert p_e \Vert_L^2{(\Omega)}$', '$\Vert \mathbf{u}_e \Vert_{H^1(\Omega)}$', '$\Vert \mathbf{u}_e \Vert_{a}$']

    for ix, error in enumerate(errors):
        utils.write_convergence_table(error, ffile_names[ix], variable_names[ix], params, gammas, meshes, hs, hgammas)

    # Plot solutions to paraview
    var_names = ['u_a', 'p_a', 'u_h', 'p_h']
    file_list = [File(code_path + '/Results/Plots/' + var_name + '.pvd') for var_name in var_names]
    var_list = [u_a_i, p_a_i, u, p]
    for ix, file in enumerate(file_list):
        file << var_list[ix]
    from IPython import embed;embed()
    root = Path(code_path)
    xdmf_file_list = [XDMFFile(str(((root / "Results" / "Plots" / var_name).with_suffix(".xdmf")))) for var_name in var_names]
    for file, name, var in zip(xdmf_file_list, var_names, var_list):
        file.write_checkpoint(var, name, 0.0, append=False)
    for file in xdmf_file_list:
        file.close()

    # Plot solutions using matplotlib
    import matplotlib.pyplot as plt
    for func, name in [[p_a_i, 'p_a'], [p, 'p_h']]:
        fig = plt.figure()
        c = plot(func, mode='color')
        plt.colorbar(c)
        fig.savefig(code_path + '/Results/Plots/' + name + '.png')

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--cyl_refinement', help='whether or not to use a refined mesh', default=1, type=int)

    parser.add_argument('--nrefs', help='Number of refinements', default=4, type=int)

    parser.add_argument('--beta', help='friction parameter', default=-2.0, type=float)

    parser.add_argument('--nu', help='viscosity', default=1.0, type=float)

    parser.add_argument('--degree', help='polynomial degree, k=1 gives the Taylor Hood P2xP1 elements', default=2, type=int)

    parser.add_argument('--ngammas', help='Number of Nitsche stabilization parameters to try', default=1, type=int)

    parser.add_argument('--normal_choice', help='str giving normal type, options are projected_normal and discrete_normal', default='projected_normal', type=str)

    parser.add_argument('--verbose', help='whether to print Newton output and errors', default=0, type=int)

    args = parser.parse_args()
    params = {
        'beta': args.beta,
        'nu': args.nu,
        'gammas': args.ngammas,
        'degree': args.degree,
        'normal_choice': args.normal_choice,
        'cyl_refinement': args.cyl_refinement
    }

    assert params['degree']>1, "Polynomial degree for Taylor-Hood elements needs to be >1"
    assert params['normal_choice'] == 'projected_normal' or params['normal_choice'] == 'discrete_normal', "Normal choices available are 'projected_normal' and 'discrete_normal'"

    N = args.nrefs

    if not args.verbose:
        set_log_level(LogLevel.ERROR) #suppress Newton solver output in fenics

    # Run convergence tests
    run_convergence_test(N, params)




