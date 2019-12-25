# Molecular Dynamics Method: Theory and Implementation

<a href="https://github.com/kaityo256/md2019"> <div class="btn-square"><i class="fab fa-github"></i> View on GitHub</div></a>

## About this document

This is a lecture note for an intensive lecture at Kanazawa University (from Dec. 18 to 20, 2019). 

[HTML](https://kaityo256.github.io/md2019/)

[PDF (Written in Japanese)](https://kaityo256.github.io/md2019/md2019.pdf)

## Table of Contents

All notes are written in Japanese and slides are written in English.

* 0\. [Introduction](about/README.md) [[Slides]](https://speakerdeck.com/kaityo256/md2019-introduction)
    * What is MD?
    * Questions about MD
    * Purpose of this lecture
* [1. Classical Mechanics](basic/README.md)
    * 1.1 Euler-Lagrange equation
    * 1.2 Hamiltonian's equation
    * 1.3 Liouville Operator
    * 1.4 Variables and Observables
* [2. Pressure](pressure/README.md)
    * 2.1 Global Pressure
    * 2.2 Local Stress
* [3. Temperature](temperature/README.md)
    * 3.1 Maxwell Distribution
    * 3.2 Canonical Ditribution
    * 3.3 Generarized Virial Theorem
* [4. Numerical Integration](integration/README.md)
    * 4.1 Integration of ODE
    * 4.2 Integration of Equations of Motion
    * 4.3 Symplectic Integrator
* [5. Nose-Hoover Method](nosehoover/README.md)
    * 5.1 Temperature Controll
    * 5.2 Nose-Hoover Method
    * 5.3 Problems on Nose-Hoover method
* [6. Langevin Thermostat](langevin/README.md)
    * 6.1 Langevin Equation
    * 6.2 Euler-Maruyama Method
    * 6.3 H Theorem
* [7. Integration scheme for non-Hamliton systems](respa/README.md)
    * 7.1 Non-Hermiticity of Liouville Operator
    * 7.2 RESPA
    * 7.3 7.3 Time Reversibility
* [8. Generalized Liouville's Theorem of non-Hamiltonian systems](liouville/README.md)
    * 8.1 Jacobi's Formula
    * 8.2 Dynamics of Jacobian
    * 8.3 Generalized Liouvllie's Theorem
* 9\. Implementations and Optimization [[Slides]](https://speakerdeck.com/kaityo256/md2019-implementations-and-optimization)
    * 9.1 Architecture of Computer
    * 9.2 Memory Optimization
    * 9.3 Architecture Dependent Optimization
* 10\. Programming Design [[Slides]](https://speakerdeck.com/kaityo256/md2019-programming-design)
    * 10.1 Module Coupling
    * 10.2 Design for Parallelization

## 謝辞

中村壮伸さんにIrving-Kirkwoodによる微視的圧力定義、およびEuler-Maruyamaの方法を教えていただきました。

## 参考文献

* 大学演習「熱学・統計力学」 久保亮五編 裳華房
* [Tuckerman教授(NYU)の講義ノート](http://www.nyu.edu/classes/tuckerman/stat.mechII/lectures.html)
