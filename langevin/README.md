# Langevin Thermostat

Nose-Hoover法は決定論的な熱浴であったが、次は確率的な熱浴であるランジュバン熱浴を考えてみよう。

## Langevin Equation

以下のようなランジュバン方程式を考える。

$$
\begin{aligned}
\dot{p} &= -\frac{\partial H}{\partial q} - \gamma \frac{\partial H}{\partial p} + \sqrt{2D} \hat{R} \\
\dot{q} &= \frac{\partial H}{\partial p}
\end{aligned}
$$

ただし、$\gamma$は定数、$\hat{R}$は

$$
\left<R(t_1)R(t_2)\right> = \delta(t1-t2)
$$

を満たすランダムな力(White Noise)である。

さて、この方程式が$H$に関するカノニカル分布

$$
f \sim \exp\left(-\beta H\right)
$$

を定常状態に持って欲しい。そこで、Nose-Hoover法と同様に分布関数$f$に関する連続の式、

$$
\frac{\partial f}{\partial t} = - \nabla \cdot \vec{J}
$$

を考える。$\vec{J}$を求めるためにKramers–Moyal展開をしよう。

$$
\frac{\partial f}{\partial t} 
= \sum_{n=0}^\infty
\frac{1}{n!}
(- \nabla)^n
\left(C_n f\right)
$$

ただし、$C_n$は、遷移確率のモーメントであり、

$$
C_n(\vec{\Gamma}) = \int (\vec{\Gamma}'-\vec{\Gamma})^n W(\vec{\Gamma}'| \vec{\Gamma}) d \vec{\Gamma}'
$$

が定義である。

一次のモーメント$C_1$はドリフト項と呼ばれ、ランダム力以外の項が残る。

$$
\nabla(C_1 f)= 
\left[
\frac{\partial}{\partial p}
\left(
-\frac{\partial H}{\partial q} - \gamma \frac{\partial H}{\partial p}\right)
+ 
\frac{\partial}{\partial q}
\left(
\frac{\partial H}{\partial q}
\right)
\right]f
$$