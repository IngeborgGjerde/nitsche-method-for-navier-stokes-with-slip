from fenics import *
import numpy as np

radius = 1.0

def mesh_size_surface(facet_f, tag):
    hmin = 1E15
    mesh = facet_f.mesh()
    tdim = mesh.topology().dim()
    assert tdim == 2  # So that tdim-1 are edges
    mesh.init(tdim-1, 0)
    
    x = mesh.coordinates()
    for facet in facets(mesh):
        if facet_f[facet] == tag:
            hmin = min(np.linalg.norm(np.diff(x[facet.entities(0)])), hmin)
    return hmin

def get_meshes(Ns, use_refined_mesh):
    from pathlib import Path
    code_path = str(Path.cwd()) + '/Meshing/'

    meshes, facet_functions, hs, hgammas = [], [], [], []
    
    for scale in Ns:  
        if use_refined_mesh:
            path = code_path + 'RefinedMeshes/mesh%g' % scale
        else:
            path = code_path + 'UniformMeshes/mesh%g' % scale


        mesh = Mesh(path + '.' +  'xml')
        facetfunction = MeshFunction("size_t", mesh, path + "_facet_region.xml")
        
        meshes.append( mesh )
        hs.append(mesh.hmax())

        facet_functions.append(facetfunction)
        hgammas.append( mesh_size_surface(facetfunction, 1) )

    File(code_path + 'bmarkers.pvd') << facet_functions[0]

    return meshes, facet_functions, hs, hgammas