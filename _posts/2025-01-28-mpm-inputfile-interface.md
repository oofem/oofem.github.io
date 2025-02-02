---
title: "First steps with symbolic mpm module (input file syntax)"
date: 2025-01-28
author: Borek Patzak
categories:
  - blog
tags: tutorial mpm
---

## First steps with new OOFEM multi-physics module (MPM)

What is MPM?
MPM is new multi-physic OOFEM module, developed to make implementation of multi-physics problems more simple. In this tutorial we will focus on symbolic MPM extension. This extension is somehow different from traditional OOFEM modules.

In traditional OOFEM problem setup, the user discretizes the problem using problem-specific elements, boundary conditions, etc. that provide necessary functionality to solve the problem.
On the other hand, the symbolic MPM module allows to use universal, problem independent elements. These elements just define geometry. What is going to be evaluated is defined using **Terms**, that are evaluated on elements and integrated using **Integrals**. **Terms** can be evaluated on any element using user defined **Interpolations**. In this way, the problem definition is more like defining the weak form of the problem. 
     
The basic design ideas behind MPM are following:
* MPM defines **Interpolations**, **Variables**, **Terms**, **Integrals** as reusable blocks to create multiphysics formulations
* **Variables** represent unknown fields (or test felds) in a weak form of the problem. The variable has its interpolation, type (scalar, vector) and physical meaning defined.
* **Terms** represent an integrand in weak formulation to be evaluated (integrated). Term definition is independent on underlying element geometry and interpolation. It defines two key methods:
    *  method to evaluate the term value, typically contributing to RHS of discrete system.
    *  method to evaluate the consistent linearization of the term, so if Term is T(u), depending on unknown u, this term evaluates dT/du, which typically contributes to the LHS. 
* **Interpolations** are representing different FE interpolations.
* **Integral** represents the integral of term in a weak form. It can compute its contributions to the discrete set of equations. 

The concept allows for parametrization with different element geometries and interpolations. Also, the components (interpolations, terms) can be reused/shared between different formulations.
In the [previous post on mpm module](https://oofem.github.io/blog/mpm-introduction/) I have introduced the Python interface of MPM module. In this post, I will illustrate mpm features on example set-up from traditional OOFEM input deck.  

We will consider the same linear 2D elasticity problem, defined by following weak form of equilibrium equations:

$$ \int_\Omega (\partial w)^T \sigma (\partial u)\ d\Omega = \int_\Omega w^T \rho g d\Omega + \int_\Gamma w^T t d\Gamma $$

where
* u,w are variables (fileds), represented by _Variable_ class instances
* $\left[ (\partial w)^T \sigma (\partial u)\right] $ and $\left[ (w)^Tt \right]$ are Terms, parametrized (to be evaluated) by u,w, represented by classes _BTSigmaTerm_ and _NTfTerm_ (derived from parent _Term_ class).

When term is evaluated, the interpolations of the test and unknown fields as well as the element geometry are substituted. In the following notation, the approximations on the element level are expressed as $u^e=\sum N_u r_u$ and $w^e=\sum N_w r_w$. 


### Simple example

We will use OOFEM input deck to demonstrate the concept of setting up the problem of cantilever beam fixed on the left-hand side and loaded by distributed loading on free, right-hand edge.

```
     y↑ 
      |6              7               8               9                10
      +---------------+---------------+---------------+---------------+  +­­­­­­­­­­↑
      |               |               |               |               |  |↑
      |               |               |               |               |  |↑
      |       1       |       2       |       3       |       4       | 5|↑
      |               |               |               |               |  |↑
      |               |               |               |               |  |↑
      +---------------+---------------+---------------+---------------+  +↑ fy=1.0   --->x
      1               2               3               4                 5
```
In the following I illustrate the setup of the problem using traditional oofem input deck syntax, which has been extended to allow for definition of Variables, Terms, and Integrals. If you are already familiar with OOFEM input syntax you will conceptually understand, the details and syntax can be found in OOFEM input manual. For presentation purposes, the complete input deck has been broken into parts to allow for intermediate comments:
```
demo.out
Demo of symbolic mpm problem; bending of clamped cantilever (l=3, h=0.3) loaded at the free end
#
# Note the new input sections for Terms, Variables and Integrals
# this exampkle requires oofem to be compiled with MPM support (USE_MPM=ON)
#  
test nsteps 1 nvariables 2 nterms 2 nintegrals 2 lhsterms 1 1 rhsterms 1 2 nmodules 1
#vtkxml primvars 1 1 tstep_all
errorcheck
```
Introduce test and unknown vector fields with linear approximation space:
``` 
Variable name "u" interpolation "feilin" type 1 quantity 0 size 2 dofs 2 1 2 # displacement 
Variable name "w" interpolation "feilin" type 1 quantity 0 size 2 dofs 2 1 2 # test function
```
Set up Terms appearing in the weak form
```
BTSigmaTerm 1 variable "u"  testvariable "w" mmode 6
NTfTerm 2 variable "u" testvariable "w" mmode 6 flux 2 0. 1.
```
Set up two integrals (one over the volume, second over the boundary where load is applied)
```
Integral 1 domain 1 set 1 term 1
Integral 2 domain 1 set 2 term 2
```
The input file continues with traditional records. Note that we are using universal elements (quads q1 and lines l1)
```
domain HeatTransfer
outputmanager tstep_all dofman_all element_all
ndofman 10 nelem 5 nbc 2 ncrosssect 1 nic 0 nltf 2 nmat 1 nset 4
# Set up nodes
node 1  coords 3 0.            0.           0.
node 2  coords 3 0.75          0.           0.
node 3  coords 3 1.5           0.           0.
node 4  coords 3 2.25          0.           0.
node 5  coords 3 3.0           0.           0.
node 6  coords 3 0.            0.3          0.
node 7  coords 3 0.75          0.3          0.
node 8  coords 3 1.5           0.3          0.
node 9  coords 3 2.25          0.3          0.
node 10 coords 3 3.0           0.3          0.
# Set up universal (quad) elements
q1 1 nodes 4 1 2 7 6  mat 1 crosssect 1
q1 2 nodes 4 2 3 8 7  mat 1 crosssect 1
q1 3 nodes 4 3 4 9 8  mat 1 crosssect 1
q1 4 nodes 4 4 5 10 9 mat 1 crosssect 1
# boundary element 
l1 5 nodes 2 5 10     mat 1 crosssect 1
# Traditional set up of materials, cross sections boundary conditions and sets
simplecs 1 thick 1.0
isole 1 d 1 e 1 n 0.3 talpha 1.
# x-displacement
boundarycondition 1 loadtimefunction 1 set 3 values 1 0   dofs 1 1 
# y-displacement
boundarycondition 2 loadtimefunction 1 set 4 values 1 0   dofs 1 2 
#
#
constantfunction 1 f(t) 1
PiecewiseLinFunction 2 nPoints 4 t 4 -10. 0. 1. 5. f(t) 4 0. 0. 1.0 1.0 
set 1 elements 4 1 2 3 4
# 
set 2 elements 1 5
# x-bc, note that set 3 is defined using element edge
set 3 elementedges 2 1 4
# y-bc, note that set 3 is defined using element edge
#set 4 nodes 3 1 6 11
set 4 elementedges 2 1 4 
```
Save the above input file into demo.in file. 
```
./oofem -f demo.in
```
The complete input deck can be found in [tests/mpm/mpms03.in](https://raw.githubusercontent.com/oofem/oofem/refs/heads/mpm2/tests/mpm/mpms03.in) file.

### Simple example - quadratic interpolation

Let's switch now to quadratic interpolation. The nice thing here is that it is sufficient to change only two lines in input file to get solved the problem using quadratic interpolation. We just need to locate the records defining test and unknown fields and update the interpolation:
```
Variable name "u" interpolation "feiquad" type 1 quantity 0 size 2 dofs 2 1 2 # displacement 
Variable name "w" interpolation "feiquad" type 1 quantity 0 size 2 dofs 2 1 2 # test function
```
All the magic needed to introduce additional nodes on shared edges, setting up integration rules, etc. will happen automatically for you.
The complete input deck can be found in [tests/mpm/mpms04.in](https://raw.githubusercontent.com/oofem/oofem/refs/heads/mpm2/tests/mpm/mpms04.in) file.

### Simple example - results & conclusions
The analytical solution (beam theory, assuming only the bending moment contribution) is   
$w_{ex}=FL^3/(3EI) = 0.3\times3^3/(3\times1\times0.3^3/12.)=1200$.

The resulting end deflection (node 10), obtained with linear approximation is $w_{lin}=345.03$.
The deflection obtained using quadratic interpolation is $w_q=1178.8$, clearly demonstrating the superior convergence properties of quadratic interpolation over the linear one.

|                |  linear approx. | quadratic approx.  | exact (beam theory) |
|----------------|-----------------|--------------------|---------------------|
| End deflection | 345.03          | 1178.8             | 1200.0              |   


However, the purpose of this post was to illustrate the power of MPM symbolic module and its capabilities. 
With this I conclude today post on mpm module. 

Hope you enjoyed and stay tuned for following updates!
You can leave a comment below to give a feedback.
