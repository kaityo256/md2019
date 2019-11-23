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

二次のモーメント$C_2$は拡散項と呼ばれ、ランダム力のみが残る。

$$
\frac{1}{2}\nabla^2(C_2f)=
\frac{\partial^2}{\partial p^2} (Df)
$$

これらを連続の式の形に書くと、

$$
\frac{\partial f}{\partial t} =
- \frac{\partial}{\partial p}
\left(
- \frac{\partial H}{\partial q}
\underbrace{- \gamma \frac{\partial H}{\partial p}
- D \frac{\partial}{\partial p}}_{(*)}
\right) f
-\frac{\partial}{\partial q}
\left(\frac{\partial H}{\partial p} \right)f
$$

定常状態$f_\mathrm{eq}$では $\partial f_\mathrm{eq}/\partial = 0$となることと、ハミルトンの運動方程式由来の項がキャンセルすることを使うと、結局(*)の項しか残らず、

$$
\left(- \gamma \frac{\partial H}{\partial p} 
- D\underbrace{\frac{\partial}{\partial p}}_{-\beta \frac{\partial H}{\partial p}}\right)f_\mathrm{eq} = 0
$$

が要請される。ここで、$f_\mathrm{eq} \sim \exp(-\beta H)$であるから、

$$
\gamma =D \beta
$$

つまり、

$$
\beta = \gamma/D
$$

これはEinsteinの関係式に他ならない。

以上から、摩擦係数$\gamma$と、拡散係数$D$の比を適切に設定すれば、指定の温度のカノニカル分布が定常状態となる。

