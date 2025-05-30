---
title: "OOFEM integrates with JupiterLite"
date: 2025-05-31
author: Borek Patzak
categories:
  - blog
tags: tutorial python JupyterLite
---

## OOFEM Now Runs in Your Browser with JupyterLite Integration

OOFEM supports Python scripting through its Python bindings, enabling users to write scripts, develop custom elements, material models, and solvers—all in Python. Thanks to recent advancements, OOFEM and its binding code can now be compiled into a WebAssembly (WASM) Python module.

By leveraging [Pyodide](https://github.com/pyodide/pyodide) (a WebAssembly-based Python distribution for browsers and Node.js) and [JupyterLite](https://jupyterlite.readthedocs.io/en/stable/) (a browser-based version of JupyterLab), you can now run and interact with OOFEM entirely within your web browser—no installation needed.

While this might seem like a lightweight demo, it actually unlocks powerful new opportunities for teaching, training, and getting started with finite element analysis. It’s especially useful for educational settings or for newcomers exploring OOFEM for the first time.

Check out the [OOFEM python tutorial](https://oofem.github.io/jupyter-demos/lab/index.html?path=Welcome.ipynb) and try out the browser-based examples directly—no setup required!

The support by [Václav Šmilauer](https://github.com/eudoxos) is highly acknowledged.

<img src="/assets/oofem-jupyterlite-screenshot.png" alt="Screenshot of JupyterLite with oofem" width="50%" height="50%">


____

With this we conclude today post on oofem JupyterLite integration. 
Hope you enjoyed and stay tuned for following updates!

You are welcome to leave a comment below to give a feedback.

