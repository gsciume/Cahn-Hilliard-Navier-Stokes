#------------------------------------------------------------------------------
#------------------------------------------------------------------------------ 
# Simulation of pipette aspiration with a Cahn-Hilliard-Navier-Stokes model
# Partial Wetting condition (theta = 150)
#
# Reference paper:
# "A bi-component model to assess the rheology of soft cellular aggregates 
# probed using the micropipette aspiration technique"
# Authored by Giuseppe Sciume, Karine Guevorkian, Pierre Nassoy
#
# Author of the code: GS
# email: giuseppe.sciume@u-bordeaux.fr
# Version : 28/08/2024
# 
# Repository : https://github.com/gsciume/Cahn-Hilliard-Navier-Stokes
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

from fenics import *
from mshr   import *
from ufl import *
from dolfin import *
from fenics import *
from math import floor, ceil
import random as random
import numpy as np
import time
import csv
import ufl
import os
import shutil

#import matplotlib.pyplot as py
from mshr import *


#----------------------------------------------------------------------
# General setting
#----------------------------------------------------------------------

# Messages
Messages = 1

# Export paraview results (1 = YES, 0 = NO)
Export_PVD = 1

# Viscoelastic rheology (1 = Viscoelastic, 0 = Newtonian)
Viscoelastic = 1

# Stabilized (Stabilized = 1)
Stabilized = 1

# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True
parameters["allow_extrapolation"]           = True

# Definition of files for outputs
file_mesh     = File('mesh.pvd')
file_boundary = File('boundary.pvd')
file0         = File("potential.pvd"  ,  "compressed")
file1         = File("cells.pvd"      ,  "compressed")
file2         = File("pressure.pvd"   ,  "compressed")
file3         = File("velocity.pvd"   ,  "compressed")
file4         = File("strain_rate.pvd",  "compressed")
file5         = File("stress.pvd"     ,  "compressed")
file6         = File("stresso.pvd"    ,  "compressed")


#----------------------------------------------------------------------
# Procedures 
#----------------------------------------------------------------------
class InitialConditionsNE(UserExpression):
	def __init__(self, xxc, yyc, r, we, **kwargs):
#		random.seed(5 + MPI.rank(MPI.comm_world))
		self.r   = r
		self.we = we
		self.xc  = xxc
		self.yc  = yyc
#		super().__init__(**kwargs)
	def eval(self, values, x):
		if sqrt((x[0]-self.xc)*(x[0]-self.xc)+(x[1]-self.yc)*(x[1]-self.yc)) < self.r :
			values[0] = 0.
			values[1] = self.we
			values[2] = 0.
			values[3] = 0.
			values[4] = 0.
		else:
			values[0] = 0.
			values[1] = 0.
			values[2] = 0.
			values[3] = 0.
			values[4] = 0.
	def value_shape(self):
		return(5,)

# u_init = SphericalInitialConditions(1.7*H1, 0., RIC, cequ, degree = 1)
# Create intial conditions and interpolate phase-field
class InitialConditionsVE(UserExpression):
	def __init__(self, xxc, yyc, r, we, **kwargs):
		random.seed(14 + MPI.rank(MPI.comm_world))
		self.r   = r
		self.we = we
		self.xc  = xxc
		self.yc  = yyc
		super().__init__(**kwargs)
	def eval(self, values, x):
		if sqrt((x[0]-self.xc)*(x[0]-self.xc)+(x[1]-self.yc)*(x[1]-self.yc)) < self.r :
			values[0] = 0.
			values[1] = self.we
			values[2] = 0.
			values[3] = 0.
			values[4] = 0.
			values[5] = 0.
			values[6] = 0.
			values[7] = 0.
			values[8] = 0.
			values[9] = 0.
			values[10] = 0.
			values[11] = 0.
			values[12] = 0.
			values[13] = 0.
            
		else:
			values[0] = 0.
			values[1] = 0.
			values[2] = 0.
			values[3] = 0.
			values[4] = 0.
			values[5] = 0.
			values[6] = 0.
			values[7] = 0.
			values[8] = 0.
			values[9] = 0.
			values[10] = 0.
			values[11] = 0.
			values[12] = 0.
			values[13] = 0.
	def value_shape(self):
		return(14,)


def Max_height(ch_, mesh):
    chsub1 = ch_.split()[1]
    F = Function(FunctionSpace(mesh, 'CG', 1))
    F.interpolate(Expression('c>0.4?sqrt(x[0]*x[0]):0', degree = 1, c = chsub1))
    ymax = np.max(F.vector().get_local())
    return ymax
    
    
def Output_pvd(sol, t, strain_out, stress_out, Viscoelastic, TENSORE):
    file0 << (sol.split()[0], t/60)
    file1 << (sol.split()[1], t/60)
    file2 << (sol.split()[2], t/60)
    file3 << (sol.split()[3], t/60)
    file4 << (strain_out,     t/60)
    file5 << (stress_out,     t/60)
    
    if Viscoelastic == 1:
        tau_i  = sol.split()[4]
        tau_p = project(tau_i, TENSORE)
        taup_out.assign(tau_p)
        file6 << (taup_out,     t/60)
      
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Data of the experiment (Karine paper)
#----------------------------------------------------------------------

# Rayon pipette
R_PIPE = 35.e-6
T_PIPE = 20.e-6

# Rayon initial agregat cellulaire
RIC = 175.e-6
XC = 1.75*RIC 
YC = 0.

# Other geometrical parameters
L1 = 1800.e-6
H1 = 1.4 * RIC
L2 = L1 - (4.*H1)
H2 = T_PIPE  
L3 = L2 - (0.5*H1)
H3 = H1 - R_PIPE - T_PIPE
RC = 0.5*T_PIPE 
FACT = 1.10

# Pression aspiration
p_suction = Constant(-1180.)


#----------------------------------------------------------------------
# Model parameters
#----------------------------------------------------------------------

# Equilibrium mass fraction cells 
cequ     = 0.80

# Wetting angle (180 correspond to the no-wet condition)
angle = 150  

# Epsilon
epsilon  = 1.4e-6

# Density of the medium species
rho_m = 1000. 

# Density of the cell species
rho_c = 1000. 

# Dynamic viscosity culture medium
eta_medium  = 1600. 

# Identified parameters
eta_cells_maxwell = 102825.
G_maxwell         = 19.28
eta_cells_sliping = 17600.
sigma_cm          = 0.00995
Mc                = 3.796e-16

# In the Newtonian case we assume an infinite stiffness of the spring 
# so the total viscosity of the cells is eta_cells_sliping + eta_cells_maxwell
if Viscoelastic == 0:
    eta_cells_sliping = eta_cells_sliping + eta_cells_maxwell

# Output message 
if Messages == 1:
    print('General setting and parameters') 
   
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Construction of the FE mesh
#----------------------------------------------------------------------
R1 = Rectangle(Point(0., 0.), Point(L1, H1))
R2 = Rectangle(Point(2*H1, R_PIPE), Point((2*H1 + L2), (R_PIPE+H2)))
R3 = Rectangle(Point((2.25*H1), (R_PIPE+T_PIPE)), Point((2.25*H1 + L3), (R_PIPE+T_PIPE+H3)))
R4 = Rectangle(Point(0., 0.), Point(L1, R_PIPE)) 
R5 = Rectangle(Point(0., (R_PIPE + T_PIPE)), Point((2.1*H1), (R_PIPE + 1.1*T_PIPE))) 
R6 = Rectangle(Point((L1 - 2.1*H1), (R_PIPE + T_PIPE)), Point(L1, (R_PIPE + 1.1*T_PIPE)))
R7 = Rectangle(Point(0.5*L1, 0.), Point(L1, H1)) 

C1 = Circle(Point(2*H1,(R_PIPE + 0.5*T_PIPE)), (FACT*RC), 16)
C2 = Circle(Point((2*H1 + L2),(R_PIPE + 0.5*T_PIPE)), (FACT*RC), 16)
domain1 = R1 - R2 - R3 - C1 - C2 - R4 - R5 - R6
domain2 = domain1 + R4 + R5 + R6 - R7

mesh   = generate_mesh(domain2,  40)

# Mesh refinement in the needed areas
ref_t = 15.e-6
markers = MeshFunction("bool", mesh, 2)
markers.set_all(False)
for c in cells(mesh):
    for f in facets(c):
        if (sqrt((f.midpoint()[0] - XC)*(f.midpoint()[0] - XC)+(f.midpoint()[1] - YC)*(f.midpoint()[1] - YC)) < (RIC + ref_t) and \
           sqrt((f.midpoint()[0] - XC)*(f.midpoint()[0] - XC)+(f.midpoint()[1] - YC)*(f.midpoint()[1] - YC)) > (RIC - ref_t)) \
           or (f.midpoint()[0] > (0.22*L1) and f.midpoint()[1] < R_PIPE):
            markers[c] = True

mesh = refine(mesh, markers, redistribute=False)

# Mesh writting in file
file_mesh << mesh

# Output message 
if Messages == 1:
    print('FE meshing') 
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Time discretization of the problem
#----------------------------------------------------------------------
minute  = 60.
TFIN0   = 30. * minute
TFIN1   = TFIN0 + (180. * minute)
TFIN2   = TFIN1 + (180. * minute)

# Time step
dt0 = 10.
dt1 = 10.
dt2 = 10.

# Output message 
if Messages == 1:
    print('Time dicretizaiton') 
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Definition of functional space, trial & test functions
#----------------------------------------------------------------------
# Define function space
CP1  = FiniteElement("CG", mesh.ufl_cell(), 1)
CV2  = VectorElement("CG", mesh.ufl_cell(), 2)

# Newtonian case
CHS  = FunctionSpace(mesh, MixedElement(CP1, CP1, CP1, CV2))

# Viscoelastic case 
if Viscoelastic == 1:
    CT3  = TensorElement("CG", mesh.ufl_cell(), 1, (3,3), symmetry = True)
    CHS  = FunctionSpace(mesh, MixedElement(CP1, CP1, CP1, CV2, CT3))

# Tensor space (for stress and strain)
TENSOR = FunctionSpace(mesh, MixedElement(CP1, CP1, CP1, CP1, CP1, CP1, CP1, CP1, CP1))

# Trial functions and solution vectors (n and n+1) 
dsol  = TrialFunction(CHS)
sol   = Function(CHS)       # Current solution
sol0  = Function(CHS)       # Previous solution

if Viscoelastic == 0:
    m,  c,  p,  u        = split(sol)
    m0, c0, p0, u0       = split(sol0)

if Viscoelastic == 1:
    m,  c,  p,  u,  tau  = split(sol)
    m0, c0, p0, u0, tau0 = split(sol0)
    
# Test functions for the weak form
if Viscoelastic == 0:
    v_m, v_c, q, w     = TestFunctions(CHS)

if Viscoelastic == 1:
    v_m, v_c, q, w, B  = TestFunctions(CHS)

# Definition of the normal vector
n = FacetNormal(mesh)
x = SpatialCoordinate(mesh)

# Output message 
if Messages == 1:
    print('Definition of trial and test function') 
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of relevant bounds
# This time with boundary markers to extract subdomains
#----------------------------------------------------------------------

boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
boundary_markers.set_all(0)
tol  = 1.0e-7
tolc = 5.0e-6

# Boundary_1(Left bound)
class Boundary_1(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], 0, tol)

bound_1 = Boundary_1()
bound_1.mark(boundary_markers, 1)

# Boundary_2U(upper bound big channel)
class Boundary_2U(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], H1 , tol)

bound_2U = Boundary_2U()
bound_2U.mark(boundary_markers, 2)

# Boundary_3U(Left and right mid vertical bound )
class Boundary_3U(SubDomain):
	def inside(self, x, on_boundary):
		return (near(x[0], ((L1 - L3)/2), tol) and (x[1] >= (R_PIPE + T_PIPE ))) \
                or (near(x[0], (L1 - L3)/2 + L3, tol) and (x[1] >= (R_PIPE + T_PIPE )))

bound_3U = Boundary_3U()
bound_3U.mark(boundary_markers, 3)

# Boundary_4U(Left and right upper bound small channel)
class Boundary_4U(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], (R_PIPE + T_PIPE ) , tol) and (x[0] >= (0.5*(L1 - L2))) and ((0.5*(L1 - L2) + L2) >= x[0])

bound_4U = Boundary_4U()
bound_4U.mark(boundary_markers, 4)

# Boundary_5U(Left and right upper circular bound)
xc1 = 2*H1
yc1 = R_PIPE + 0.5*T_PIPE
xc2 = 2*H1 + L2
yc2 = R_PIPE + 0.5*T_PIPE

class Boundary_5U(SubDomain):
	def inside(self, x, on_boundary):
		return (abs(sqrt((x[1]-yc1)*(x[1]-yc1)+(x[0]-xc1)*(x[0]-xc1))-(FACT*RC)) < tolc and on_boundary) \
                or (abs(sqrt((x[1]-yc2)*(x[1]-yc2)+(x[0]-xc2)*(x[0]-xc2))-(FACT*RC)) < tolc and on_boundary)

bound_5U = Boundary_5U()
bound_5U.mark(boundary_markers, 5)

# Boundary_6U(upper bound small channel)
class Boundary_6U(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], R_PIPE, tol) and (x[0] > (0.5*(L1 - L2))) and (x[0] < (0.5*(L1 - L2) + L2))

bound_6U = Boundary_6U()
bound_6U.mark(boundary_markers, 6)

# Boundary_7(Right bound)
x_maxx = 0.5*L1

class Boundary_7(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], x_maxx, tol)

bound_7 = Boundary_7()
bound_7.mark(boundary_markers, 7)

# Boundary_8 (symmetry plan)
class Boundary_8(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 0., tol)

bound_8 = Boundary_8()
bound_8.mark(boundary_markers, 8)

# Boundaries writting in file
file_boundary << boundary_markers

# Redefiniction dx and ds
dx = Measure("dx", domain = mesh)
ds = Measure("ds", domain = mesh, subdomain_data = boundary_markers)

# Output message 
if Messages == 1:
    print('Definition of relevant bounds') 
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Definition of dirichlet boundary conditions
#----------------------------------------------------------------------
v_null = Constant((0.0, 0.0))
zero   = Constant(0)

# No-slip boundary condition for velocity
bcNS2    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 2)
bcNS3    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 3)
bcNS4    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 4)
bcNS5    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 5)
bcNS6    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 6)

# outflow boundary condition for velocity
bcNS7y = DirichletBC(CHS.sub(3).sub(1), zero,   boundary_markers, 7)

# inflow vy = 0
bcNS1y = DirichletBC(CHS.sub(3).sub(1), zero, boundary_markers, 1)

# Symmetry plan vy = 0
bcNS8y = DirichletBC(CHS.sub(3).sub(1), zero, boundary_markers, 8)

# Collect boundary conditions
bc_tot = [bcNS2, bcNS3, bcNS4, bcNS5, bcNS6, bcNS7y, bcNS1y, bcNS8y]

# No-wet condition if angle =180:
if angle == 180:
    bc_ch4  = DirichletBC(CHS.sub(1), zero, boundary_markers, 4) 
    bc_ch5  = DirichletBC(CHS.sub(1), zero, boundary_markers, 5)
    bc_ch6  = DirichletBC(CHS.sub(1), zero, boundary_markers, 6)
    bc_tot = [bcNS2, bcNS3, bcNS4, bcNS5, bcNS6, bcNS7y, bcNS1y, bcNS8y, bc_ch4, bc_ch5, bc_ch6]
    
# Output message 
if Messages == 1:
    print('Definition boundary conditions') 
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Definition of Initial condition
#----------------------------------------------------------------------
if Viscoelastic == 0:
    u_init = InitialConditionsNE(XC, YC, RIC, cequ, degree = 1)

if Viscoelastic == 1:
    u_init = InitialConditionsVE(XC, YC, RIC, cequ, degree = 1)

sol.interpolate(u_init)
sol0.interpolate(u_init)

# Output message 
if Messages == 1:
    print('Definition initial conditions') 
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Definition of divergence and gradient in cylindrycal coordinates
#----------------------------------------------------------------------

import dolfin as df
I = df.Identity(3)

r = abs(x[1])

def grad_cyl(u):
	return as_tensor([[u[0].dx(0), u[0].dx(1), 0.], [u[1].dx(0), u[1].dx(1), 0.], [0., 0., u[1]/r]])

def div_cyl(u):
	return u[1]/r + u[0].dx(0) + u[1].dx(1)

def strain_rate(u): 
	return 0.5*(grad_cyl(u) + grad_cyl(u).T)

def sigma(u,p,eta_i): 
	return -1.*p*I + 2*eta_i*strain_rate(u)

# Output message 
if Messages == 1:
    print('Definition of divergence and gradient in cylindrycal coordinates') 
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of SUPG Coefficient (for stabilized solution)
#----------------------------------------------------------------------
# This correspond to the average mesh size / average velocity magnitude
# Ref. : "Implementing multiphysics models in FEniCS: Viscoelastic flows,
#         poroelasticity, and tumor growth" Birkan Tunç et al. (2023)
SUPG_C = 100.
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Solution of the nonlinear problem (PHASE 0 - for initial solution)
#
# p_right = 0
#
# During this phase of 30 minutes we wait that the aggregate reaches
# the equilibrium condition
#
# NB : during the phase 0 we do not take into account viscoelasticity 
#----------------------------------------------------------------------
#----------------------------------------------------------------------

if Messages == 1:
    print('------------------------------------------------------------------') 
    print('Solution of the nonlinear problem (PHASE 0 - for initial solution)')
    print('------------------------------------------------------------------')    

p_right = Constant(0.)
dti     = dt0

stress_out = Function(TENSOR)
strain_out = Function(TENSOR)
taup_out   = Function(TENSOR)
TENSORE    = TensorFunctionSpace(mesh, 'CG', 1, (3,3))

# Wilson-theta
theta   = 0.5

# Interface regularization parameter
alfa     = (6*sqrt(2))/(cequ**3.)

# Form of the chemical potential
beta = 1./4.
dfdc = 2*beta*c*(cequ-c)**2 - 2*(cequ-c)*beta*(c**2)

# Theta-consistent quantities 
m_mid  = (1.0-theta)*m0 + theta*m
u_mid  = (1.0-theta)*u0 + theta*u
c_mid  = (1.0-theta)*c0 + theta*c

if Viscoelastic == 1:
    tau_mid = (1.0-theta)*tau0 + theta*tau

# eta_newton is the average Newtonian viscosity of the mixture Eqn 
eta_newton = eta_medium + c_mid*(eta_cells_sliping - eta_medium)

# The medium and cells have the same density = 1000 km/m3
rho_mixture = 1000.       

eta_maxw = c_mid*eta_cells_maxwell
lambdac  = eta_maxw/G_maxwell

# Cahn-Hilliard Equations
L0 = c*v_m*r*dx - c0*v_m*r*dx + dti*Mc*dot(grad(m_mid), grad(v_m))*r*dx + dti*dot(u_mid, grad(c_mid))*v_m*r*dx
L1 = (1.0/(epsilon*epsilon))*m*(epsilon/(alfa*sigma_cm))*v_c*r*dx - (1.0/(epsilon*epsilon))*dfdc*v_c*r*dx - dot(grad(c), grad(v_c))*r*dx + (Constant((tan(angle*np.pi/180))**-1.0))*abs(dot(as_vector([-1*n[1], n[0]]),grad(c)))*v_c*r*ds(4) + (Constant((tan(angle*np.pi/180))**-1.0))*abs(dot(as_vector([-1*n[1], n[0]]),grad(c)))*v_c*r*ds(5) + (Constant((tan(angle*np.pi/180))**-1.0))*abs(dot(as_vector([-1*n[1], n[0]]),grad(c)))*v_c*r*ds(6)

# For the no-wet case
if angle == 180:
    L1 = (1.0/(epsilon*epsilon))*m*(epsilon/(alfa*sigma_cm))*v_c*r*dx - (1.0/(epsilon*epsilon))*dfdc*v_c*r*dx - dot(grad(c), grad(v_c))*r*dx 

# Navier-Stokes equations (divergence form)
if Viscoelastic == 0:
    L2 =  rho_mixture*dot(u, w)*(1./dti)*r*dx - rho_mixture*dot(u0, w)*(1./dti)*r*dx + eta_newton*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta_newton*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx - m_mid*dot(grad(c_mid), w)*r*dx +  p_right*dot(n, w)*r*ds(7)

if Viscoelastic == 1:
    L2 =  rho_mixture*dot(u, w)*(1./dti)*r*dx - rho_mixture*dot(u0, w)*(1./dti)*r*dx + eta_newton*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta_newton*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx - m_mid*dot(grad(c_mid), w)*r*dx +  p_right*dot(n, w)*r*ds(7) + inner(tau, grad_cyl(w))*r*dx

L3 = q*div_cyl(u)*r*dx 

# Oldroyd equation
# For the phase 0 we do not take into account viscoelasticity
lambdac0 = 0.*lambdac

if Viscoelastic == 1:
    L4 = inner((tau + lambdac0*((tau - tau0)*(1./dti) + dot(u, nabla_grad(tau)) \
    - dot((grad_cyl(u)).T, tau) - dot(tau, grad_cyl(u))) -                \
    eta_maxw*(grad_cyl(u) + (grad_cyl(u).T))), B)*r*dx 

# Assembling of the system of eqs  in weak form
if Viscoelastic == 0:
    L = L0 + L1 + L2 + L3
    
if Viscoelastic == 1:
    L = L0 + L1 + L2 + L3 + L4
    
# Output message 
if Messages == 1:
    print('Assembling of the system of eqs in weak form') 
#----------------------------------------------------------------------


# Compute directional derivative about u in the direction of du (Jacobian)
a = derivative(L, sol, dsol)

# Create nonlinear problem and Newton solver
problem = NonlinearVariationalProblem(L, sol, bcs=bc_tot, J=a)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "lu"
solver.parameters["newton_solver"]["convergence_criterion"] = "incremental"
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-6
solver.parameters["newton_solver"]['maximum_iterations'] = 30

# Temps and lenght arrays
tps = np.array([]) ; radi = np.array([])

# Step in time
t = 0.
i = 1

# Calculation of strain rate and stress tensors
ciii = sol.split()[1]
piii = sol.split()[2] 
uiii = sol.split()[3]
eta_i = eta_medium + ciii*(eta_cells_sliping - eta_medium)

dddi         = strain_rate(uiii)
sssi         = sigma(uiii,piii,eta_i)
strain_ratei = project(dddi,TENSORE)
stressi      = project(sssi,TENSORE)
stress_out.assign(stressi) 
strain_out.assign(strain_ratei) 

if Export_PVD == 1 :
    Output_pvd(sol, t, strain_out, stress_out, Viscoelastic, TENSORE)


while (t < TFIN0):

	t += dti
	print('t =', t, 'dt = ', dti, 'iteration numero', i)

	sol0.vector()[:] = sol.vector()
	solver.solve()

	# File outputs each 6 time steps (if i%6 == 0:)
	if i%6 == 0:
		piii = sol.split()[2]
		uiii = sol.split()[3]
		ciii = sol.split()[1]
		eta_i = eta_medium + ciii*(eta_cells_sliping - eta_medium)
		
		dddi = strain_rate(uiii)
		sssi = sigma(uiii,piii,eta_i)
		strain_ratei = project(dddi,TENSORE)
		stressi      = project(sssi,TENSORE)
		stress_out.assign(stressi)
		strain_out.assign(strain_ratei)
		
		if Export_PVD == 1 : Output_pvd(sol, t, strain_out, stress_out, Viscoelastic, TENSORE)

	i += 1;



#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Solution of the nonlinear problem (PHASE 1 - aspiration)
#
# p_right = p_suction
#----------------------------------------------------------------------
#----------------------------------------------------------------------

if Messages == 1:
    print('------------------------------------------------------------------') 
    print('Solution of the nonlinear problem (PHASE 1 - aspiration)          ')
    print('------------------------------------------------------------------')   
    
p_right = p_suction
dti     = dt1

if Viscoelastic == 0:
    m0, c0, p0, u0  = sol0.split()
    
if Viscoelastic == 1:
    m0, c0, p0, u0, tau0  = sol0.split()

# Ridefinition of certain weak form (update of time step and pressure on the right side)
# Cahn-Hilliard Equations
L0 = c*v_m*r*dx - c0*v_m*r*dx + dti*Mc*dot(grad(m_mid), grad(v_m))*r*dx + dti*dot(u_mid, grad(c_mid))*v_m*r*dx

# Navier-Stokes equations (divergence form)
if Viscoelastic == 0:
    L2 =  rho_mixture*dot(u, w)*(1./dti)*r*dx - rho_mixture*dot(u0, w)*(1./dti)*r*dx + eta_newton*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta_newton*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx - m_mid*dot(grad(c_mid), w)*r*dx +  p_right*dot(n, w)*r*ds(7)

if Viscoelastic == 1:
    L2 =  rho_mixture*dot(u, w)*(1./dti)*r*dx - rho_mixture*dot(u0, w)*(1./dti)*r*dx + eta_newton*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta_newton*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx - m_mid*dot(grad(c_mid), w)*r*dx +  p_right*dot(n, w)*r*ds(7) + inner(tau, grad_cyl(w))*r*dx

# Oldroyd equation

if Viscoelastic == 1:
    L4 = inner((tau + lambdac*((tau - tau0)*(1./dti) + dot(u, nabla_grad(tau)) \
    - dot((grad_cyl(u)).T, tau) - dot(tau, grad_cyl(u))) -                \
    eta_maxw*(grad_cyl(u) + (grad_cyl(u).T))), B)*r*dx 

if (Viscoelastic == 1) and (Stabilized == 1):
    L4 = inner((tau + lambdac*((tau - tau0)*(1./dti) + dot(u, nabla_grad(tau)) \
    - dot((grad_cyl(u)).T, tau) - dot(tau, grad_cyl(u))) -                \
    eta_maxw*(grad_cyl(u) + (grad_cyl(u).T))), (B + dot(SUPG_C*u, nabla_grad(B))))*r*dx
    print('Stabilized solution') 

# Assembling of the system of eqs  in weak form
if Viscoelastic == 0:
    L   = L0 + L1 + L2 + L3
    
if Viscoelastic == 1:
    L   = L0 + L1 + L2 + L3 + L4
# Output message 
if Messages == 1:
    print('Assembling of the system of eqs in weak form') 
#----------------------------------------------------------------------

# Compute directional derivative about u in the direction of du (Jacobian)
a = derivative(L, sol, dsol)

# Create nonlinear problem and Newton solver
problem = NonlinearVariationalProblem(L, sol, bcs=bc_tot, J=a)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "lu"
solver.parameters["newton_solver"]["convergence_criterion"] = "incremental"
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-6
solver.parameters["newton_solver"]['maximum_iterations'] = 30

# Calculation of strain rate and stress tensors
piii = sol.split()[2] 
uiii = sol.split()[3]
ciii = sol.split()[1]
eta_i = eta_medium + ciii*(eta_cells_sliping - eta_medium)

dddi = strain_rate(uiii)
sssi = sigma(uiii,piii,eta_i)
strain_ratei = project(dddi,TENSORE)
stressi      = project(sssi,TENSORE)
stress_out.assign(stressi) 
strain_out.assign(strain_ratei) 

if Export_PVD == 1 : Output_pvd(sol, t, strain_out, stress_out, Viscoelastic, TENSORE)

while (t < TFIN1):

	t += dti
	print('t =', t, 'dt = ', dti, 'iteration numero', i)

	sol0.vector()[:] = sol.vector()
	solver.solve()

	# File outputs each 6 time steps (if i%6 == 0:)
	if i%6 == 0:
		piii = sol.split()[2]
		uiii = sol.split()[3]
		ciii = sol.split()[1]
		eta_i = eta_medium + ciii*(eta_cells_sliping - eta_medium)
		
		dddi = strain_rate(uiii)
		sssi = sigma(uiii,piii,eta_i)
		strain_ratei = project(dddi,TENSORE)
		stressi      = project(sssi,TENSORE)
		stress_out.assign(stressi)
		strain_out.assign(strain_ratei)
		
		if Export_PVD == 1 : Output_pvd(sol, t, strain_out, stress_out, Viscoelastic, TENSORE)

    # At each time step, store time, aggregate radius and pressure
	tps = np.append(tps, (t/60.)); radi = np.append(radi, (Max_height(sol, mesh)*1.E6))
	np.savetxt("sauvegarde.csv", radi, delimiter=",")
    
	i += 1;

with open('aspiration_phase1.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(['tps', 'height'])
        for i in range (np.size(tps)):
            spamwriter.writerow([tps[i], radi[i]])


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Solution of the nonlinear problem (PHASE 2 - retraction)
#
# p_right = 0
#----------------------------------------------------------------------
#----------------------------------------------------------------------

if Messages == 1:
    print('------------------------------------------------------------------') 
    print('Solution of the nonlinear problem (PHASE 2 - retraction)          ')
    print('------------------------------------------------------------------')   
    
p_right = Constant(0.)
dti     = dt2

if Viscoelastic == 0:
    m0, c0, p0, u0  = sol0.split()
    
if Viscoelastic == 1:
    m0, c0, p0, u0, tau0  = sol0.split()

# Ridefinition of certain weak form (update of time step and pressure on the right side)
# Cahn-Hilliard Equations
L0 = c*v_m*r*dx - c0*v_m*r*dx + dti*Mc*dot(grad(m_mid), grad(v_m))*r*dx + dti*dot(u_mid, grad(c_mid))*v_m*r*dx

# Navier-Stokes equations (divergence form)
if Viscoelastic == 0:
    L2 =  rho_mixture*dot(u, w)*(1./dti)*r*dx - rho_mixture*dot(u0, w)*(1./dti)*r*dx + eta_newton*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta_newton*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx - m_mid*dot(grad(c_mid), w)*r*dx +  p_right*dot(n, w)*r*ds(7)

if Viscoelastic == 1:
    L2 =  rho_mixture*dot(u, w)*(1./dti)*r*dx - rho_mixture*dot(u0, w)*(1./dti)*r*dx + eta_newton*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta_newton*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx - m_mid*dot(grad(c_mid), w)*r*dx +  p_right*dot(n, w)*r*ds(7) + inner(tau, grad_cyl(w))*r*dx

# Oldroyd equation
if Viscoelastic == 1:
    L4 = inner((tau + lambdac*((tau - tau0)*(1./dti) + dot(u, nabla_grad(tau)) \
    - dot((grad_cyl(u)).T, tau) - dot(tau, grad_cyl(u))) -                \
    eta_maxw*(grad_cyl(u) + (grad_cyl(u).T))), B)*r*dx 

if (Viscoelastic == 1) and (Stabilized == 1):
    L4 = inner((tau + lambdac*((tau - tau0)*(1./dti) + dot(u, nabla_grad(tau)) \
    - dot((grad_cyl(u)).T, tau) - dot(tau, grad_cyl(u))) -                \
    eta_maxw*(grad_cyl(u) + (grad_cyl(u).T))), (B + dot(SUPG_C*u, nabla_grad(B))))*r*dx
    print('Stabilized solution') 

# Assembling of the system of eqs  in weak form
if Viscoelastic == 0:
    L   = L0 + L1 + L2 + L3
    
if Viscoelastic == 1:
    L   = L0 + L1 + L2 + L3 + L4
# Output message 
if Messages == 1:
    print('Assembling of the system of eqs in weak form') 
#----------------------------------------------------------------------


# Compute directional derivative about u in the direction of du (Jacobian)
a = derivative(L, sol, dsol)

# Create nonlinear problem and Newton solver
problem = NonlinearVariationalProblem(L, sol, bcs=bc_tot, J=a)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "lu"
solver.parameters["newton_solver"]["convergence_criterion"] = "incremental"
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-6
solver.parameters["newton_solver"]['maximum_iterations'] = 30

# Calculation of strain rate and stress tensors
piii = sol.split()[2] 
uiii = sol.split()[3]
ciii = sol.split()[1]
eta_i = eta_medium + ciii*(eta_cells_sliping - eta_medium)

dddi = strain_rate(uiii)
sssi = sigma(uiii,piii,eta_i)
strain_ratei = project(dddi,TENSORE)
stressi      = project(sssi,TENSORE)
stress_out.assign(stressi) 
strain_out.assign(strain_ratei) 

if Export_PVD == 1 : Output_pvd(sol, t, strain_out, stress_out, Viscoelastic, TENSORE)

while (t < TFIN2):

	t += dti
	print('t =', t, 'dt = ', dti, 'iteration numero', i)

	sol0.vector()[:] = sol.vector()
	solver.solve()

	# File outputs each 6 time steps (if i%6 == 0:)
	if i%6 == 0:
		piii = sol.split()[2]
		uiii = sol.split()[3]
		ciii = sol.split()[1]
		eta_i = eta_medium + ciii*(eta_cells_sliping - eta_medium)
		
		dddi = strain_rate(uiii)
		sssi = sigma(uiii,piii,eta_i)
		strain_ratei = project(dddi,TENSORE)
		stressi      = project(sssi,TENSORE)
		stress_out.assign(stressi)
		strain_out.assign(strain_ratei)
		
		if Export_PVD == 1 : Output_pvd(sol, t, strain_out, stress_out, Viscoelastic, TENSORE)

    # At each time step, store time, aggregate radius and pressure
	tps = np.append(tps, (t/60.)); radi = np.append(radi, (Max_height(sol, mesh)*1.E6))
	np.savetxt("sauvegarde.csv", radi, delimiter=",")
    
	i += 1;


with open('aspiration_phase2.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(['tps', 'height'])
        for i in range (np.size(tps)):
            spamwriter.writerow([tps[i], radi[i]])

