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
\left<R(t_1)R(t_2)\right> = \delta(t_1-t_2)
$$

を満たすランダムな力(White Noise)である。

さて、この方程式がハミルトニアン$H$に関するカノニカル分布

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

定常状態$f_\mathrm{eq}$では $\partial f_\mathrm{eq}/\partial t = 0$となることと、ハミルトンの運動方程式由来の項がキャンセルする(ハミルトニアン由来の項は非圧縮流を作るため、分布関数を変化させない)ことを使うと、結局(*)の項しか残らず、

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

ランジュバン熱浴を実装するためには、シミュレーションで用いている時間刻み$h$の間だけホワイトノイズを積分しただけの力積を計算する必要がある。しかし、ホワイトノイズ(連続的な確率過程)を積分する、という処理を数学的に厳密に扱うのは難しく、真面目にやるならWiener過程を導入して、などとやるのであろうが、以下では厳密さを犠牲にして直感的な導出を試みる。

今、運動量$p$がホワイトノイズ$\sqrt{2D} \hat{R}$にさらされているとしよう。運動方程式は、

$$
\dot{p} = \sqrt{2D} \hat{R}
$$

であり、$\hat{R}$は

$$
\left<R(t_1)R(t_2)\right> = \delta(t1-t2)
$$

を満たす標準揺動力とする。この運動量の確率分布関数$f(p,t)$を考える。先程の運動方程式に対応するFocker-Plank方程式は

$$
\begin{aligned}
\frac{\partial f}{\partial t} &= -\frac{\partial}{\partial p}
\left(-D \frac{\partial}{\partial p} \right)f\\
&= D \frac{\partial^2 f}{\partial p^2}
\end{aligned}
$$

これは拡散係数$D$を持つ一次元拡散方程式に他ならない。今、時刻$t$において運動量が$p_0$であったとすると、時刻$t+h$の分布は

$$
f(p, t+h) = \mathcal{N}(p_0, \underbrace{2Dh}_{\sigma^2})
$$

つまり、平均$p_0$、分散$2Dh$となるガウス分布となる。ここで、Einsteinの関係式から

$$
D = k_B \gamma T 
$$

$k_B=1$とする単位系を取れば、最終的に

$$
f(p, t+h) = \mathcal{N}(p_0, 2 \gamma T h)
$$

となる。実装では、Langevin部分は一次のオイラー法を用いて

$$
p(t+h) = p(t) - \gamma p(t)h + \hat{w}
$$

とすることが多い。ただし、$\hat{w}$は平均$0$、分散$2\gamma T h$のガウス分布に従う乱数であり、例えばC++を使うなら、

```c++
std::normal_distribution<double> nd(0.0, 2.0 * gamma * T * h);
std::mt19937 mt;
double w = nd(mt);
```

などとして生成することができる。なお、この手法をEuler-Maruyamaの方法と呼ぶ。

## H Theorem

Nose-Hoover法が保証するのは、「位相空間をボルツマン重みに比例する確率で走る」ということのみであり、さらに運動がエルゴード的であって初めてカノニカル分布が達成される。先ほど、Langevin方程式もカノニカル分布を定常状態に持つことを示したが、Langevin系の場合はNose-Hoover系よりもう少し強いことが言える。すなわち、Langevin系では、いかなる初期条件から初めても、必ずカノニカル分布に収束することを示すことができる。以下、それを見てみよう。簡単のため、一自由度系を考え、位相空間$(p,q)$上に、ハミルトニアン$H(p,q)$が定義されているとしよう。

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

を満たす。このうち、ハミルトニアン由来の項はエントロピーの増加に寄与しない(しつこいが、非圧縮流れを作るため)ため、それ以外を代入すると、

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

以上から、Langevin系は時間発展に伴って必ず自由エネルギーが単調減少することが示された。自由エネルギーが変化しなくなった場合は、定常状態としてカノニカル分布に収束する。

ここで、自由エネルギーは、現在の分布と、カノニカル分布とのKullback–Leibler (KL)距離になっていることに注意したい。実際、

$$
\begin{aligned}
D_{KL}(f | f_\mathrm{eq} ) &\equiv \int f \ln \frac{f}{f_\mathrm{eq}} dpdq \\
&= \int (f \ln f - f \ln f_\mathrm{eq})dpdq \\
&= \int (f \ln f + \beta H f)dpdq + C \\
&= \beta F + C
\end{aligned}
$$

ただし$C$は定数である。ランジュバン系では、分布関数とカノニカル分布のKL距離が単調に減少し、最終的に距離がゼロ、つまりカノニカル分布が実現した時が定常状態であることがわかる。
