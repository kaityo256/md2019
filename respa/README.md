# Integration scheme for non-Hamliton systems

常微分方程式の数値積分には多数の方法があるが、分子動力学法、すなわちハミルトンの運動方程式の積分には、ほとんどの場合においてシンプレクティック積分が用いられている。その理由は、シンプレクティック積分が軌道には誤差を持ちつつも、エネルギーを厳密に保存するからであった。シンプレクティック積分は

* リュービル演算子がエルミートになっている
* リュービル演算子が、エルミート演算子の和に分解できる
* 分解した演算子それぞれについて、指数関数の肩に乗せた時間発展演算子が厳密に計算できる

という性質を使い、指数分解の公式により数値積分法を構築する手法である。

さて、


## Thermostat Term

簡単のため、一自由度系を考える。また、熱浴の質量も1としておこう。運動方程式は以下の通りである。

$$
\begin{aligned}
\dot{p} &= \underbrace{-\frac{\partial H}{\partial q}}_{iL_K} \underbrace{- p\zeta}_{iL_T} \\
\dot{q} &= \underbrace{\frac{\partial H}{\partial p}}_{iL_V} \\
\dot{\zeta} &= \underbrace{\left(p\frac{\partial H}{\partial p} - \frac{1}{\beta}\right)}_{iL_Z}
\end{aligned}
$$

この方程式は以下のカノニカル分布を定常状態に持つ。

$$
f = Z^{-1} \exp\left[ -\beta\left( H + \frac{\zeta^2}{2}\right) \right]
$$

この系のリュービル演算子$iL$は、以下のように記述できる。

$$
iL = \dot{p} \frac{\partial}{\partial p} +
\dot{q} \frac{\partial}{\partial p} + 
\dot{\zeta} \frac{\partial}{\partial \zeta}
$$

数値積分法を構築するとは、これを指数関数の肩に乗せた演算子

$$
U(h) = \mathrm{e}^{ih L}
$$

をなんらかの手段で近似、評価することである。

シンプレクティック積分の場合と同様に、リュービル演算子を以下のように分解しよう。

$$
iL = iL_K + iL_V + iL_T + iL_Z
$$

これまで、リュービル演算子を分解したとき、指数関数の二次以上の項が消えるために厳密に計算ができることを利用してシンプレクティック積分を構築していた。

たとえば、ハミルトンの運動方程式の運動項由来の部分は

$$
\begin{aligned}
U_K(h) &= \mathrm{e}^{i h L_K} \\
&= \exp\left(-h \frac{\partial H}{\partial q} \frac{\partial}{\partial p} \right)\\
&= \sum_k \frac{1}{k!}\left( -h\frac{\partial H}{\partial q} \frac{\partial}{\partial p}\right)^n
\end{aligned}
$$

となっている。$p$による偏微分があるため、

$$
iL_K q = iL_K \zeta = 0
$$

になるのはすぐにわかる。問題は$p$に演算した場合だが、$p$の偏微分の左側に$p$依存性がないのがポイントである。

つまり、リュービル演算子を一度かけると、$p$依存性が消えてしまう。

$$
\begin{aligned}
iL_K p &= \left(-h \frac{\partial H}{\partial q} \frac{\partial }{\partial p}\right) p \\
&= -h \frac{\partial H}{\partial q}
\end{aligned}
$$

したがって、リュービル演算子を二回かけるとゼロになる。

$$
\begin{aligned}
(iL_K)^n p &= iL_K \left( -h \frac{\partial H}{\partial q} \right) \\
&=0
\end{aligned}
$$

結局、

$$
\begin{aligned}
U_K(h) p &= \mathrm{e}^{i h L_P} p \\
& =(1 + ih L_K) p\\
&= p - h \frac{\partial H}{\partial q}
\end{aligned}
$$

となり、単純な一次のオイラー法に帰着する。

同様に$iL_V, iL_Z$についても、指数関数の肩に乗せると二次の項が消える。問題は摩擦項$iL_T$である。


$$
iL_T = - p\zeta \frac{\partial}{\partial p}
$$

これを指数関数の肩に乗せた演算子は

$$
\begin{aligned}
U_T(h) &= \mathrm{e}^{i h L_T} \\
&= \exp\left(-h p\zeta \frac{\partial}{\partial p} \right)
\end{aligned}
$$

となる。$p$による偏微分を含むため、$q$や$\zeta$にかけると0になる。問題は$p$にかけた場合である。この部分リュービル演算子$iL_T$を$p$にかけてやると

$$
iL_T p = -hp\zeta
$$

となり、$p$依存性が残る。したがって、$iL_T$をもう一度かけてもゼロとはならない。

$$
(iL_T)^2 p = h^2 p \zeta^2
$$

高次項はゼロとはならないのだが、$n$回かけたものが簡単に計算できる。

$$
(iL_T)^n p = (-h\zeta)^n p
$$

これにより、