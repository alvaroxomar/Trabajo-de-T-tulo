from sfepy import data_dir
from sfepy.discrete import Problem
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem import FEDomain
from sfepy.discrete.fem import Field
from sfepy.terms import Term
from sfepy.discrete.conditions import Conditions, EssentialBC
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.discrete.variables import FieldVariable
import numpy as np

# Crear la malla rectangular en el dominio [0, 1] x [0, 1]
mesh = Mesh.from_file(data_dir + '/meshes/2d/square_unit_tri.mesh')
domain = FEDomain('domain', mesh)

# Definir el campo escalar sobre la malla, usando elementos de orden 1 (P1)
field = Field.from_args('phi', np.float64, 'scalar', domain, approx_order=1)

# Definir la variable dependiente (solución) y la variable de prueba
phi = FieldVariable('phi', 'unknown', field)
v = FieldVariable('v', 'test', field, primary_var_name='phi')

# Definir la ecuación de Poisson en forma débil: -laplace(phi) = f
# Con el término de fuente f = -6
f = -6.0
t = Term.new('dw_laplace(v, phi)', integral, domain, v=v, u=phi) + Term.new('dw_volume_lvf(v, val)', integral, domain, v=v, val=f)

# Definir condiciones de borde:
# phi(x, 0) = 0 en el borde inferior
bc_inferior = EssentialBC('bc_inferior', domain.get_boundaries()['Bottom'], {'phi.0': 0.0})

# phi(0.5, 0.5) = 1 en el nodo central
bc_central = EssentialBC('bc_central', domain.get_vertices()[np.where(np.allclose(mesh.coors, [0.5, 0.5], atol=1e-2))], {'phi.0': 1.0})

# Agrupar las condiciones de borde
bcs = Conditions([bc_inferior, bc_central])

# Configurar el problema y resolver
problem = Problem('poisson', equations={'balance': t}, nls=Newton(), ls=ScipyDirect())
problem.set_bcs(ebcs=bcs)
state = problem.solve()

# Graficar la solución usando Matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# Extraer la solución y las coordenadas para graficar
phi_values = state.get_state_parts()['phi']
x, y = mesh.coors[:, 0], mesh.coors[:, 1]
triangulation = tri.Triangulation(x, y, mesh.get_conn())

plt.figure()
plt.tripcolor(triangulation, phi_values, shading='gouraud')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Solución de la ecuación de Poisson usando SfePy")
plt.show()
