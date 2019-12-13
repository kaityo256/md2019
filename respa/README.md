# Integration scheme for non-Hamliton systems

常微分方程式の数値積分には多数の方法があるが、分子動力学法、すなわちハミルトンの運動方程式の積分には、ほとんどの場合においてシンプレクティック積分が用いられている。その理由は、シンプレクティック積分が軌道には誤差を持ちつつも、エネルギーを厳密に保存するからであった。シンプレクティック積分は

* リュービル演算子がエルミートになっている
* リュービル演算子が、エルミート演算子の和に分解できる
* 分解した演算子それぞれについて、指数関数の肩に乗せた時間発展演算子が厳密に計算できる

という性質を使い、指数分解の公式により数値積分法を構築する手法である。こうして構築されたシンプレクティック積分は、位相空間の体積を厳密に保存する。

では、温度制御が入った場合はどうだろうか？温度制御が入るということは、本質的に分布関数が揺らぐため、位相空間の体積は保存しない。この時、シンプレクティック積分と同様な数値積分法が構築できるだろうか？

ここでは、シンプレクティック積分と同様に指数分解の方法を使って数値積分法を構築するRESPAと呼ばれる手法と、その性質について紹介する。

## Thermostat Term

簡単のため、一自由度系を考える。また、熱浴の質量も1としておこう。ハミルトニアン$H$にNose-Hoover熱浴をつけた運動方程式は以下の通りである。

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

ハミルトンの運動方程式では、リュービル演算子を分解したとき、指数関数の二次以上の項が消えるために厳密に計算ができることを利用してシンプレクティック積分を構築していた。

たとえば、ハミルトンの運動方程式の運動項由来の部分は

$$
\begin{aligned}
U_K(h) &= \mathrm{e}^{i h L_K} \\
&= \exp\left(-h \frac{\partial H}{\partial q} \frac{\partial}{\partial p} \right)\\
&= \sum_k \frac{1}{k!}\left( -h\frac{\partial H}{\partial q} \frac{\partial}{\partial p}\right)^k
\end{aligned}
$$

となっている。$p$による偏微分があるため、$p$や$\zeta$にかけると、

$$
(iL_K) q = (iL_K) \zeta = 0
$$

と、それぞれ0になるのはすぐにわかる。問題は$p$に演算した場合だが、$p$の偏微分の左側に$p$依存性がないのがポイントで、リュービル演算子を一度かけると、$p$依存性が消えてしまう。

$$
\begin{aligned}
iL_K p &= \left(-h \frac{\partial H}{\partial q} \frac{\partial }{\partial p}\right) p \\
&= -h \frac{\partial H}{\partial q}
\end{aligned}
$$

したがって、リュービル演算子をもう一度かけるとゼロになる。

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

と、一次のオイラー法のような更新が導かれる。

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

しかし、高次項はゼロとはならないのだが、$k$回かけたものが簡単に計算できる。

$$
(iL_T)^k p = (-h\zeta)^k p
$$

これにより、$U_T(h)p$が厳密に計算できてしまう。

$$
\begin{aligned}
U_T(h) p &= \exp(iL_K) p\\
&= \sum_k \frac{1}{k!} 
\underbrace{\left(-hp\zeta \frac{\partial}{\partial p} \right)^k p}_{(-h\zeta)^kp}\\
&= p \underbrace{\sum_k \frac{(-h\zeta)^k}{k!}}_{\exp(-h\zeta)}\\
&= \mathrm{e}^{-h\zeta} p
\end{aligned}
$$

これにより、リュービル演算子を分解した4つの要素を指数関数の肩に乗せたものが全て厳密に評価できた。これを使ってシンプレクティック積分と同様に数値積分法を構築するのがRESPA (reference system propagator algorithm)である。特に時間反転対称にしたものをr-RESPAというが、現在では単にRESPAというと時間反転対称にしたものを指すことが多い。

RESPAの構築方法には複数あるが、通常は二次の対称分解の形として、真ん中に最も計算が重い力の計算を持ってくることが多い。

$$
\tilde{U}(h)_\mathrm{RESPA} = \mathrm{e}^{ihL_K/2} \mathrm{e}^{ihL_T/2} \mathrm{e}^{ihL_V} \mathrm{e}^{ihL_T/2} \mathrm{e}^{ihL_K/2}
$$

この時間発展は

1. 位置を現在の速度で$h/2$だけ進める
2. 運動量を$\zeta$を使って$p \rightarrow p \mathrm{e}^{-h\zeta/2}$とスケールする
3. 現在の位置において力を計算し、運動量を$h$だけ更新する
4. 運動量を$\zeta$を使って$p \rightarrow p \mathrm{e}^{-h\zeta/2}$とスケールする
5. 位置を現在の速度で$h/2$だけ進める

という手続きになり、二次精度で、かつ時間反転対称となっている。ただし、$L_T$がエルミート演算子ではないために、全体としてシンプレクティック変換にはなっていない。

## Time Reversibility

### Linear System

ここで、時間発展演算子の時間反転対称性についてまとめておこう。時間を$h$だけ進める時間発展演算子$U(h)$が時間反転対称であるとは、

$$
U^{-1}(h) = U(-h)
$$

を満たすことを言う。

この式の意味を見るために、調和振動子系を考えて見よう。いま、時刻$t$において$(p,q)$にいた系が、時間発展により時刻$t+h$で$(P,Q)$に移ったとする。すなわち

$$
\begin{pmatrix}
P \\
Q
\end{pmatrix}
= 
U(h)
\begin{pmatrix}
p \\
q
\end{pmatrix}
$$

この時、$U(h)$は行列になる。$(P,Q)$は$(p,q)$の関数であるが、それを逆に解いて、$(p,q)$を$(P,Q)$で表してみよう。線形の場合は逆行列をかけるだけだ。

$$
\begin{pmatrix}
p \\
q
\end{pmatrix}
= 
U^{-1}(h)
\begin{pmatrix}
P \\
Q
\end{pmatrix}
$$

さて、時間反転対称であるとは、時間を$h$だけ進める演算子を演算した結果を逆に解いたら、それは時間を$-h$だけ進める演算子をかけた結果と等しい、ということである。従って、

$$
\begin{pmatrix}
p \\
q
\end{pmatrix}
= 
U(-h)
\begin{pmatrix}
P \\
Q
\end{pmatrix}
$$

さて、よく誤解されているのだが、シンプレクティック積分は時間反転対称性を持つとは限らない。実際、一次のシンプレクティック積分に対応する行列

$$
\tilde{U}_1(h) = 
\begin{pmatrix}
1 -h^2 & -h \\
h  & 1
\end{pmatrix}
$$

の逆行列は

$$
\tilde{U}^h_1(h) = 
\begin{pmatrix}
1  & h \\
-h  & 1 - h^2
\end{pmatrix}
$$

であり、これは明らかに$\tilde{U}(-h)$とは一致しない。

一方、二次のシンプレクティック積分に対応する行列

$$
\tilde{U}_2(h)
= 
\begin{pmatrix}
1 - h^2/2 & h \\
-h + h^3/4 & 1 - h^2/2
\end{pmatrix}
$$

逆行列は

$$
\tilde{U}^{-1}_2(h)
= 
\begin{pmatrix}
1 - h^2/2 & h \\
h - h^3/4 & 1 - h^2/2
\end{pmatrix}
$$

となり、これは$\tilde{U}_2(-h)$に一致する。

### Operator Formulation

時間発展演算子を指数分解の形で書くと、時間反転対称性がわかりやすい。以下のハミルトンダイナミクスを考える。

$$
\begin{aligned}
\dot{p} &= - \frac{\partial H}{\partial q} \\
\dot{q} &= \frac{\partial H}{\partial p}
\end{aligned}
$$

対応するリュービル演算子は以下の通り。

$$
iL = 
\underbrace{-\frac{\partial H}{\partial q} \frac{\partial}{\partial p}}_{iL_V}
+
\underbrace{\frac{\partial H}{\partial p} \frac{\partial}{\partial q}}_{iL_K}
$$

まず、厳密な時間発展演算子を考えよう。

$$
U(h) = \mathrm{e}^{ih L}
$$

この時間発展演算子で時刻$t$において$(p,q)$であった状態が時刻$t+h$において$(P,Q)$に移動したのなら

$$
\begin{aligned}
\begin{pmatrix}
P\\
Q
\end{pmatrix}
&= U(h) 
\begin{pmatrix}
p\\
q
\end{pmatrix} \\
&= 
\mathrm{e}^{ihL}
\begin{pmatrix}
p\\
q
\end{pmatrix} 
\end{aligned}
$$

$(p,q)$を$(P,Q)$で表すには、時間発展演算子を逆側に持ってくれば良いから、

$$
\begin{aligned}
\begin{pmatrix}
p\\
q
\end{pmatrix}
&= 
\mathrm{e}^{-ihL}
\begin{pmatrix}
P\\
Q
\end{pmatrix} \\
&\equiv U^{-1}(h) 
\begin{pmatrix}
p\\
q
\end{pmatrix}
\end{aligned}
$$

ここから、

$$
U^{-1}(h) = U(-h)
$$

が成立するのは自明であろう。要するに、時間発展演算子の時間反転対称性は、指数関数を移項すると負符号がつくことに対応している。

次に、時間発展演算子を指数分解公式で近似する場合を考えよう。一次のシンプレクティック積分に対応する時間発展演算子は

$$
\tilde{U}_1(h) = \mathrm{e}^{ihL_K}\mathrm{e}^{ihL_V}
$$

となる。

$$
\begin{aligned}
\begin{pmatrix}
P\\
Q
\end{pmatrix}
&= \tilde{U}_1(h) 
\begin{pmatrix}
p\\
q
\end{pmatrix} \\
&= 
\mathrm{e}^{ihL_K}\mathrm{e}^{ihL_V}
\begin{pmatrix}
p\\
q
\end{pmatrix} 
\end{aligned}
$$

$(p,q)$について解くと、

$$
\begin{aligned}
\begin{pmatrix}
p\\
q
\end{pmatrix}
&= 
\mathrm{e}^{-ihL_V}\mathrm{e}^{-ihL_K}
\begin{pmatrix}
P\\
Q
\end{pmatrix} \\
&\equiv
\tilde{U}_1^{-1}(h) 
\begin{pmatrix}
P\\
Q
\end{pmatrix} \\
\end{aligned}
$$

と、$\mathrm{e}^{-ihL_V}$と$\mathrm{e}^{-ihL_K}$の位置が入れ替わる。一般に$iL_K$と$iL_V$は交換しないため、$\tilde{U}_1^{-1}(h)$が$\tilde{U}_1(-h)$と等しくないことがわかるであろう。ただし、この演算子$\tilde{U}_1(h)$はユニタリ演算子であるから、この演算子による時間発展はシンプレクティックになる。なぜなら分解した部分リュービル演算子は全てエルミートであり、それに虚数単位$i$をつけて指数関数の肩に乗せたものはユニタリであり、ユニタリ演算子の積はユニタリになるからだ。

次に、二次のシンプレクティック積分を考えよう。時間発展演算子を対称分解すると

$$
\tilde{U}_2(h) = \mathrm{e}^{ihL_K/2}\mathrm{e}^{ihL_V}\mathrm{e}^{ihL_K/2}
$$

これが、時間反転対称であるのは自明であろう。演算子が左右対称な形をしているため、移項しても形が変わらない。従って、

$$
\tilde{U}_2^{-1}(h) = \tilde{U}_2(-h)
$$

が成り立つ。また、近似された時間発展演算子はユニタリ演算子の積で構築されているため、全体としてやはりユニタリになっており、位相空間体積を保存する。すなわちシンプレクティック積分となっている。

さて、Nose-Hoover熱浴がついた系にRESPAを適用すると、

$$
\tilde{U}(h) = \mathrm{e}^{ihL_K/2} \mathrm{e}^{ihL_T/2} \mathrm{e}^{ihL_V} \mathrm{e}^{ihL_T/2} \mathrm{e}^{ihL_K/2}
$$

となるのであった。二次の対称分解の場合と同様に、これが時間反転対称性を持っていることは明らかであろう。しかし、途中で非ユニタリな項$\mathrm{e}^{ihL_T/2}$を含むため、全体としても非ユニタリになっている。

以上見てきた通り、シンプレクティック積分においては、一次の場合は時間反転対称ではなく、二次の対称分解の場合は時間反転対称であった。そして、どちらもシンプレクティックであった。ハミルトンダイナミクスにシンプレクティック積分を適用するとエネルギーが保存するのだが、それと時間反転対称性は関係がない(時間反転対称でなくてもシンプレクティックであり得る)。

二次のシンプレクティック積分と同様に、Nose-Hoover熱浴がついた系にRESPAを適用すると、時間反転対称な時間発展を得ることができる。しかし、温度制御された系において時間発展が時間反転対称性を持つことがどれだけうれしいかは、さほど自明ではない。

## Non-Hermiticity of Liouville Operator

ハミルトンダイナミクスにおいてはリュービル演算子がエルミートになり、さらにエルミート性から「位相空間の流れ」が非圧縮となることを見た。確率流が非圧縮であることから、分布関数が不変になること、すなわち(エルゴード的であれば)ミクロカノニカルであることが結論されるのであった。

では、定常状態としてカノニカル分布を持つような系のリュービル演算子がどのような性質を持つか見てみよう。

簡単のため、位相空間を $\Gamma = \vec{z} = (z_1, z_2, \cdots)$と書く。なんらかの方法により、この空間に運動方程式$\dot{\vec{z}} = (\dot{z_1},\dot{z_2},\cdots)$が導入されたとしよう。この系のリュービル演算子は

$$
-iL = \sum_i \dot{z_i} \frac{\partial }{\partial z_i}
$$

となる(虚数単位$i$と添え字が紛らわしいが、文脈で判別して欲しい)。

この空間に住む分布関数を$f$とすると、確率保存から連続の式

$$
\begin{aligned}
\frac{\partial f}{\partial t} &= 
- \nabla \cdot \left( \dot{z} f\right)\\
&= - \sum_i \frac{\partial}{\partial z_i} \left( \dot{z} f\right)
\end{aligned}
$$

定常状態としてカノニカル分布

$$
f_\mathrm{eq} = Z^{-1} \exp(-\beta H)
$$

を持つならば、

$$
\sum_i \frac{\partial }{\partial z_i} (\dot{z}_i \mathrm{e}^{-\beta H}) = 0
$$

が成り立つ必要がある。従って、

$$
\begin{aligned}
\sum_i\frac{\partial \dot{z}_i}{\partial z_i} &= \beta 
\sum_i \frac{\partial H}{\partial z_i} 
\end{aligned}
$$

が成り立つ必要がある。Nose-HooverでもKinetic-Momentsでも、Nose-Hoover-Chainでも、カノニカル分布を定常状態に持つ決定論的運動方程式は必ずこの関係式を満たしている。

さて、この式の意味を見てみよう。この位相空間に住むスカラー関数$f, g$に対して、内積$(f, g) \in \mathcal{R}$が定義されている時、リュービル演算子がエルミートであるとは、

$$
(Lf, g) = (f, Lg)
$$

が成り立つことであった。それぞれ式で書くと、

$$
\begin{aligned}
(f, Lg) &= \int d\Gamma f^* \left(i \sum_i \dot{z}_i \frac{\partial g}{\partial z_i}\right) \\
(L, g) &= \int d\Gamma \left(i \sum_i \dot{z}_i \frac{\partial f}{\partial z_i}\right)^* g
\end{aligned}
$$

となる。さて、$(f, Lg)$の式を部分積分すると、

$$
\begin{aligned}
(f, Lg) &= \int d\Gamma f^* \left(i \sum_i \dot{z}_i \frac{\partial g}{\partial z_i}\right) \\
&= -i \int d\Gamma g \sum_i  \frac{\partial }{\partial z_i}
\left( \dot{z}_i f^*\right) \\
&= -i \int d\Gamma f^* g \sum_i  \frac{\partial\dot{z}_i }{\partial z_i}
-i \int d\Gamma  g \sum_i  \dot{z}_i \frac{\partial f^* }{\partial z_i} \\
&= -i \beta \int d \Gamma f^* g \sum_i \frac{\partial H}{\partial z_i} + (Lf, g) \\
&= \left(\left[L + i\beta \sum_i\frac{\partial H}{\partial z_i} \right]f, g \right) \\
&\equiv (L^\dagger f, g)
\end{aligned}
$$

ここから、以下の関係式が導かれる。

$$
L^\dagger = L + i\beta \sum_i \frac{\partial H}{\partial z_i}
$$

ハミルトンダイナミクスの場合には、Liouville演算子がエルミート、すなわち

$$
L^\dagger = L
$$

であったことを思い出そう。カノニカル分布を定常状態に保つ運動方程式に付随するLiouville演算子は必ず非エルミートとなり、その非エルミート部分のamplitudeに(逆)温度が現れる。
