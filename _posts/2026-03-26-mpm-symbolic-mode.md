---
title: "Symbolic MPM mode"
date: 2026-03-06
author: Borek Patzak
categories:
  - blog
tags: tutorial mpm
---

## Symbolic MPM: Expression-Based Term Definition

I'm excited to share a major advancement in OOFEM's multiphysics module (mpm). The new **symbolic term feature** enables you to define individual terms directly in input files using mathematical expressions—eliminating the need for hard-coded term libraries.

### Key Innovation: Expression Compiler & Virtual Machine

The implementation includes:
- **Internal Expression Compiler**: Parses symbolic expressions with support for matrix algebra and mathematical operators
- **Virtual Machine Evaluation**: Expressions are compiled once into optimized "code" for efficient runtime evaluation
- **Rich Operator Library**: Predefined operators (functors) for divergence, gradient, field interpolation, and more

This approach combines the flexibility of symbolic definitions with the performance of compiled code execution.

## Example 1: 1D Heat Transport

Let's start with a simple example of 1D stationary heat transport:

$$\lambda\frac{d^2T}{dx^2} + Q = 0$$

$$T=\bar{T}\ \rm{on}\ \Gamma_d$$

$$q\cdot n=\bar{q}\ \rm{on}\ \Gamma_q$$


### Weak Form

The weak form can be written as

$$ \int_\Omega\underbrace{\frac{d\delta T}{dx}\lambda\frac{dT}{dx}}_{T_1}-\int_{\Gamma_q}\underbrace{\delta T\bar{q}}_{T2} + \int_\Omega \underbrace{\delta TQ}_{T3}=0$$

### Variable Definition

Define the variables appearing in the weak form:
```
Variable name "t" interpolation "feilin" type 1 quantity 0 size 1 dofs 1 10 # temperature 
Variable name "dt" interpolation "feilin" type 1 quantity 0 size 1 dofs 1 10 # test function temperature
```
### Symbolic Term Definition

Instead of relying on hard-coded terms, we use the `SymbolicTerm` directive. For term $T_1$:
```
SymbolicTerm 1 variable "t"  testvariable "dt" mmode 25 lexpression "Grad(dt,gp).T*[[1.0]]*Grad(t,gp)" rexpression "Grad(dt,gp).T*[[1.0]]*Grad(t,gp)*ru(t, cell, ts)"
```
**Expression Components:**
- `rexpression`: Defines the weak form residual term
- `lexpression`: Defines the linearization with respect to the trial variable (for Newton iterations)
- `variable`: Trial variable to solve for
- `testvariable`: Test function for weak form

**Built-in Operators Used:**
- `Grad(var, gp)`: Evaluates gradient at Gauss point (e.g., $\nabla N^{e}(x_{gp}) $)
- `ru(var, cell, ts)`: Retrieves nodal unknown values for variable at current cell and time step
- `N(var, gp)`: Interpolation matrix $N^e(x_{gp})$ at Gauss point

**Heat Source Term** $T_3$ (for $Q=1$):

Note that this term requires no linearization since it doesn't depend on the unknown field:
```
SymbolicTerm 2 variable "t"  testvariable "dt" mmode 25 lexpression "0.0" rexpression "N(w,gp).T*[[1.0]]" # constant heat source
```
The complete example is available in the [tests/regression/mpm/mpms06.in](https://raw.githubusercontent.com/oofem/oofem/refs/heads/devel/tests/regression/mpm/mpms06.in) test case.

## Example 2: Multiphysics—Incompressible Elasticity

The symbolic term concept extends naturally to multiphysics problems. Here we illustrate the mixed formulation of incompressible elasticity with pressure and displacement fields:

$$\nabla p + 2\mu\nabla\cdot\rm{dev}[\nabla^s u]+f = 0\;\rm{in}\ \Omega$$

$$\frac{1}{K}p-\nabla\cdot u=0\;\rm{in}\ \Omega$$

$$u=\bar u\;\rm{on}\ \Gamma_u $$

$$\sigma n = \bar t\;\rm{on}\ \Gamma_t$$

The weak form of the above equations can be written as
$$\int_\Omega\overbrace{2\mu\nabla^sw:\rm{dev}[\nabla^su]}^{T_1}\ d\Omega-\int_\Omega\overbrace{ \nabla\cdot w\ p}^{T_2}\ d\Omega=\int_{\Gamma_t}\overbrace{w\cdot \bar t}^{T_5}\ d\Gamma $$

$$-\int_\Omega\underbrace{q\ \nabla\cdot u}_{T_3 = T_2^T}\ d\Omega+\int_\Omega\underbrace{\frac{1}{K}q\ p}_{T_4}\ d\Omega = 0
$$

###  Weak Form of Mixed u-p Formulation

This is a classic mixed finite element formulation for incompressible media.

### Problem Variables

- $u, p$: Unknown displacement and pressure fields
- $w, q$: Corresponding test functions

### Symbolic Input Syntax

Variables are defined in the input deck as:
```
Variable name "u" interpolation "feiquad" type 1 quantity 0 size 2 dofs 2 1 2 
Variable name "w" interpolation "feiquad" type 1 quantity 0 size 2 dofs 2 1 2 
Variable name "p" interpolation "feilin"  type 0 quantity 3 size 1 dofs 1 11  
Variable name "q" interpolation "feilin" type 0 quantity 3 size 1 dofs 1 11 
```

**Field Definitions:**
- Fields ($u, w$): Quadratic interpolation (`feiquad`), vector fields (`type 1`), 2 DOFs (x and y components of displacement vector)
- Pressure fields ($p, q$): Linear interpolation (`feilin`), scalar fields (`type 0`)


### Defining Weak Form Terms

Each component of the weak form is defined as a `SymbolicTerm`:

**Term $T_1$**: Stiffness—displacement test function and displacement unknown, plain strain (_mmode 7_)
```
SymbolicTerm 1 variable "u"  testvariable "w" mmode 7 \
  lexpression "B(w,gp).T*Dm_dev(gp,ts)*B(u,gp)" \
  rexpression "B(w,gp).T*Sig_dev(u,gp,ts)"
```
Where `B(var,gp)` is the strain-displacement matrix and `Dm_dev`/`Sig_dev` are deviatoric material stiffness and stress operators.

**Term $T_2$**: Coupling—displacement test and pressure unknown, plain strain (_mmode 7_)

```
SymbolicTerm 2 variable "p"  testvariable "w" mmode 7 \
  lexpression "Div(w,gp).T*N(p,gp)" \
  rexpression "Div(w,gp).T*N(p,gp)*ru(p, cell, ts)"
```
Where `Div(var,gp)` computes divergence at the Gauss point.

**Terms $T_3$ and $T_4$**: Pressure test—displacement and pressure unknowns:
```
SymbolicTerm 3 variable "u"  testvariable "dp" mmode 7 \
  lexpression "N(dp,gp).T*Div(u,gp)" \
  rexpression "N(dp,gp).T*Div(u,gp)*ru(u, cell, ts)"

SymbolicTerm 4 variable "p"  testvariable "dp" mmode 7 ctype 27 uvmt 1 \
  lexpression "N(dp,gp).T*2.4000000000024e-7*N(p,gp)" \
  rexpression "N(dp,gp).T*2.4000000000024e-7*N(p,gp)*ru(p, cell, ts)"
```

### Integration Definition

Terms are integrated over specific domains ($\Omega,\ \Gamma_t$) using predefined sets:
```
Integral 1 domain 1 set 1 term 1
Integral 2 domain 1 set 1 term 2 factor -1.0
Integral 3 domain 1 set 1 term 3 factor -1.0
Integral 4 domain 1 set 1 term 4 
Integral 5 domain 1 set 2 term 5
```


The complete example is available in the [tests/regression/mpm/mpms_cook2_u2p1.in](https://raw.githubusercontent.com/oofem/oofem/refs/heads/devel/tests/regression/mpm/mpms_cook2_u2p1.in) test case.

## Why Symbolic Terms Matter

This symbolic approach offers significant advantages:

1. **Flexibility**: Define custom physics without modifying core libraries
2. **Readability**: Mathematical expressions directly mirror the weak form on paper
3. **Performance**: One-time compilation produces optimized VM bytecode
4. **Maintainability**: Physics definitions live in input files, not scattered across C++ code
5. **Extensibility**: Easy to add new operators and mathematical constructs

## Conclusion

The new **symbolic term feature with expression compiler and virtual machine** represents a paradigm shift for OOFEM's multiphysics capabilities. By combining mathematical expressiveness with compiled execution efficiency, it enables researchers to rapidly prototype and deploy new coupled formulations without the overhead of C++ development.

Stay tuned for further enhancements to the symbolic framework!

Feel free to share feedback or questions in the comments section below.

