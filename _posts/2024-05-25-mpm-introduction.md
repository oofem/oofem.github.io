---
title: "First steps with new OOFEM multi-physics module (mpm)"
date: 2024-05-25
author: Borek Patzak
categories:
  - blog
tags:
  - tutorial mpm
---

## First steps with new OOFEM multi-physics module (MPM)

What is MPM?
MPM is new multiphysic OOFEM module, developped to make implementation of multiphysics problems more simple. 

What are the design ideas behind MPM?
* MPM defines **Interpolations**, **Variables**, **Terms**, **Integrals** as resuable blocks to create multiphysics formulations
* **Variables** represent unknown fields (or test felds) in a weak form of the problem. The variable has its interpolation, type (scalar, vector) and physical meaning defined.
* **Terms** represent an integrand in weak formulation to be evaluated (integrated). Term definition is independent on underlying element geometry and interpolation. It defines two key methods:
    *  method to evaluate the term value, typycally contributing to RHS of discrete system.
    *  method to evaluate the consistent linearization of the term, so if Term is T(u), depending on unknown u, this term evaluates dT/du, which typically contributes to the LHS. 
* **Interpolations** are representing different FE interpolations.
* **Integral** represents the integral of term in a weak form. It can compute its contributions to the dicrete set of equations. 

The concept allows for parametrization with different element geometries and interpolations. Also, the components (interpolations, terms) can be reused/shared between different formulations.
I will ilustrate the concept using oofem python bindings, which allows for fast prototyping.

In the example below, we reuse already defined terms in oofem, but it is possible to define your own terms, interpolations, etc.

As an example, consider the weak form of equlibrium equations:

$$ \int_\Omega (\partial w)^T \sigma (\partial u)\ d\Omega = \int_\Omega w^T \rho g d\Omega + \int_\Gamma w^T t d\Gamma $$

where
* u,w are variables (fileds), represented by _Variable_ class instances
* $\left[ (\partial w)^T \sigma (\partial u)\right] $ and $\left[ (w)^Tt \right]$ are Terms, parametrized (to be evaluated) by u,w, represented by classes _BTSigmaTerm_ and _NTfTerm_ (derived from parent _Term_ class).

When term is evaluated, the interpolations of the test and uknown fields as well as the element geometry are substituted. In our example, the approximations on the element level are $u^e=\sum N_u r_u$ and $w^e=\sum N_w r_w$. As this in general yields to a nonlinear system of equations, the term on left hand side evaluates
* residual contribution for given element, esentially evaluating itself with all variables known. In our example this corresponds to evaluating $\int_\Omega^e (\partial N_w)^T\sigma(\partial N_u r_u)\ d\Omega^e$
* its linearization, cooresponding in our case to $\int_\Omega^e (\partial N_w)^T \frac{\partial \sigma}{\partial \varepsilon} (\partial N_u)\ d\Omega^e$.


### Simple example
We will use OOFEM python interface to demonstrate the concept on the problem presented above.
First set up simple mesh, consisting of single quad element, defined by four nodes. Plus standard definition of materials, boundary conditions, etc.

The demo plane-stress elasticity problem consist of two elements. One is quad element representing the domain, the second is boundary element (linear segment) representing boundary between nodes 2,3 subjected to distributed loading.
Node 1 is fixed in all directions, node 4 is fixed in x direction. 

```
     y↑ 4[0,5]         3[3,4]
      +---------------+  +→
      |               |  |→
      |               |  |→
      |       1       | 2|→
      |               |  |→
      |               |  |→
      +---------------+  +→ fx=1.0   --->x
      1 [0,0]           2[1,2]
```


```python
# Requires oofempy compiled with __MPM_MODULE ON
import sys
sys.path.extend(['/home/bp/devel/oofem.git/build', '/home/bp/devel/oofem.git/bindings/python'])
import oofempy
import util

 

# Create a new dummy problem (placeholder for our demo) with one domain.
problem = oofempy.dummyProblem(nSteps=1, outFile='test_7.out')
domain = oofempy.domain(1, 1, problem, oofempy.domainType._HeatTransferMode, tstep_all=1, dofman_all=0, element_all=0)
problem.setDomain(1, domain, True)
   
# Define nodes
n1 = oofempy.node(1, domain, coords=(0, 0, 0. ))
n2 = oofempy.node(2, domain, coords=(1., 0.0, 0. ))
n3 = oofempy.node(3, domain, coords=(1., 1., 0. ))
n4 = oofempy.node(4, domain, coords=(0, 1, 0. ))
   
# Defdine elements, note that q1 defines just element geometry.
q1 = oofempy.q1(1, domain, nodes=(1,2,3,4), mat=1, crossSect=1) # quad element #1
l1 = oofempy.l1(2, domain, nodes=(2,3), mat=1, crossSect=1)     # boundary element #2

# Dirichlet Boundary conditions
bc1 = oofempy.boundaryCondition(1, domain, loadTimeFunction=1, dofs=(1,2), values=(0.,0.), set=1)
bc2 = oofempy.boundaryCondition(2, domain, loadTimeFunction=1, dofs=(1,),  values=(0.,),   set=2)
# material and cross section
mat = oofempy.isoLE(1, domain, d=1., e=1., n=0.3, talpha=1.)
cs  = oofempy.simpleCS(1, domain, mat=1, thickness=1.0)
# time functions
ltf1 = oofempy.constantFunction(1, domain, f_t=1.0)
# some sets (groups of nodes and elements) for later use
s1 = oofempy.createSet(1, domain, nodes=(1,))
s2 = oofempy.createSet(2, domain, nodes=(4,))
s3 = oofempy.createSet(3, domain, elements=(1,))
bs1 = oofempy.createSet(4, domain, elements=(2,))
util.setupDomain(domain, nodes=(n1,n2,n3,n4), elems=(q1,l1), css=(cs,), mats=(mat,), bcs=(bc1,bc2), ics=(), ltfs=(ltf1,), sets=(s1,s2,s3,bs1))

```

Now the interesting part begins:
* Define the $u$ and $w$ fields (variables), with linear interpolation, with physical meaning of Displacement vector of size 2, with two degrees of freedom u,v (2D)
* Create instances of BTSigmaTerm, NTfTerm to be evaluated for $u,w$ in plane stress mode.



```python

interpolation = oofempy.linearinterpolation()
w= u = oofempy.Variable(interpolation, oofempy.VariableQuantity.Displacement, oofempy.VariableType.vector, 2, [1,2], None)
mt = oofempy.BTSigmaTerm(w, u, oofempy.MaterialMode._PlaneStress)
lt = oofempy.NTfTerm(w, u, oofempy.MaterialMode._PlaneStress)
tstep = problem.giveNextStep()

```

Now define integrals forming our weak form. Integral `I1` over $\Omega$, defined by set `s3` (containing element 1) and boundary integral `I2` over boundary edge 2-3, defined by set `bs1`.  
The integrals need to be initialized (allocates the needed DOFs, sets up the element integration rules, etc.) as well as the entire problem (number the equations).


```python

I1 = oofempy.Integral(domain, s3, mt)
problem.addIntegral(I1)
I1.initialize()

I2 = oofempy.Integral(domain, bs1, lt)
I2.initialize()

problem.postInitialize()
problem.forceEquationNumbering()

```

    OctreeLocalizer: init
    Spatial localizer init done





    5



We are now approaching the part where we will use all the concepts:
* Create sparse matrix instance (_lhs_) to hold stifness matrix 
* The stifness matrix is assembled by integral _I1_ instance usiing its _assemble\_lhs_ method, that will integrate the linearization of our term over all elements in domain defined by integral domain.


```python

lhs=oofempy.skyline()
lhs.buildInternalStructure(problem, 1, oofempy.EModelDefaultEquationNumbering());
I1.assemble_lhs(lhs, oofempy.EModelDefaultEquationNumbering(), tstep)
lhs.printYourself()
```

    FloatMatrix with dimensions : 5 5
     4.945e-01  -1.786e-01   5.495e-02  -1.374e-02   1.786e-01  
    -1.786e-01   4.945e-01   1.374e-02  -3.022e-01  -2.473e-01  
     5.495e-02   1.374e-02   4.945e-01   1.786e-01  -1.374e-02  
    -1.374e-02  -3.022e-01   1.786e-01   4.945e-01   5.495e-02  
     1.786e-01  -2.473e-01  -1.374e-02   5.495e-02   4.945e-01  


Assemble load vector forming right hand side of the problem


```python

rhs = oofempy.FloatArray(5)

I2.assemble_rhs(rhs, oofempy.EModelDefaultEquationNumbering(), tstep)
rhs.printYourself()

```

    FloatArray of size : 5 
     5.000e-01   0.000e+00   5.000e-01   0.000e+00   0.000e+00  


Finally, use suitable linear solver to solve for unknown displacement field.  


```python
r = oofempy.FloatArray(5)
linsolv = oofempy.ldltfactorization(domain, problem)
linsolv.solve(lhs, rhs, r)
print ("Displacement vector = ", r)
```

    Displacement vector =  <oofempy.FloatArray: {1, -2.77556e-17, 1, -0.3, -0.3, }>


Hope you enjoed!
