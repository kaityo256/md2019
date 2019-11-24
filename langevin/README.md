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

これはEinsteinの関係式と呼ばれ、ランジュバン系の温度は、摩擦係数(散逸力)と揺動力の比が決めるということを意味する。

以上から、摩擦係数$\gamma$と、拡散係数$D$の比を適切に設定すれば、指定の温度のカノニカル分布が定常状態となる。

## ランジュバン熱浴の実装

TODO: 書く？

## H Theorem

Nose-Hoover法が保証するのは、「位相空間をボルツマン重みに比例する確率で走る」ということのみであり、さらに運動がエルゴード的であって初めてカノニカル分布が達成される。先ほど、Langevin方程式もカノニカル分布を定常状態に持つことを示したが、Langevin系の場合はNose-Hoover系よりもう少し強いことが言える。すなわち、Langevin系では、いかなる初期条件から初めても、必ずカノニカル分布に収束することを示すことができる。以下、それを見てみよう。

まず、時間に依存する自由エネルギー$F$を以下のように定義しよう。

$$
F = U - TS
$$

ここで$U$は内部エネルギーであり、ハミルトニアン$H$の期待値である。

$$
U = \int  H f dp dq
$$

$S$はエントロピーで、定義は以下の通り。

$$
S = - k_B \int f \ln f dp dq
$$

以上から、自由エネルギーは

$$
F = \int  (Hf+ k_B T f\ln f) dp dq
$$

と表せる。後の便利のために両辺$k_B T$で割っておこう。

$$
\beta F = \int  (\beta Hf + f \ln f) dp dq
$$

さて、この両辺を時間で微分する。

$$
\begin{aligned}
\beta \dot{F} &= \int \left(\beta H \frac{\partial f}{\partial t} + \underbrace{\frac{\partial f}{\partial t}}_{=0} + \ln f  \frac{\partial f}{\partial t}\right)dpdq \\
&= \int \frac{\partial f}{\partial t}\left( \beta H + \ln f \right) dpdq
\end{aligned}
$$

途中で、確率の保存則

$$
\int f dpdq = 1
$$

を用いた。さて、$f$は連続の式、

$$
\frac{\partial f}{\partial t} = -
\frac{\partial }{\partial p}
\left(- \frac{\partial H}{\partial q} - \gamma  \frac{\partial H}{\partial p} - D \frac{\partial}{\partial p}\right)f
- \frac{\partial }{\partial q}\left(\frac{\partial H}{\partial p}\right)f
$$

を満たす。このうち、ハミルトニアン由来の項はエントロピーの増加に寄与しないため、それ以外を代入すると、

$$
\begin{aligned}
\beta \dot{F} &= \int \frac{\partial}{\partial p}\left[\left(\gamma \frac{\partial H}{\partial p} + D\frac{\partial}{\partial p}  \right)f \right] (\beta H + \ln f )dpdq \\
&= - \int \left(\gamma \frac{\partial H}{\partial p}f + D \frac{\partial f}{\partial p} \right)
\left(
\beta \frac{\partial H}{\partial p} + \frac{1}{f} \frac{\partial f}{\partial p}
\right) dpdq \\
&= -\int \frac{D}{f} \left( \underbrace{\frac{\gamma}{D}}_{\beta} \frac{\partial H}{\partial p}f + \frac{\partial f}{\partial f} \right)
\left(
\beta \frac{\partial H}{\partial p}f + \frac{\partial f}{\partial p}
\right)
dqdq \\
&= -\int \frac{D}{f} \left( \beta \frac{\partial H}{\partial p}f + \frac{\partial f}{\partial p}\right)^2 dpdq 
 \leq 0
\end{aligned}
$$

すなわち、自由エネルギーが単調減少することが示された。
