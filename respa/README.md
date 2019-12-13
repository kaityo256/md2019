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
&= \sum_k \frac{1}{k!}\left( -h\frac{\partial H}{\partial q} \frac{\partial}{\partial p}\right)^n
\end{aligned}
$$

となっている。$p$による偏微分があるため、

$$
iL_K q = iL_K \zeta = 0
$$

になるのはすぐにわかる。問題は$p$に演算した場合だが、$p$の偏微分の左側に$p$依存性がないのがポイントで、リュービル演算子を一度かけると、$p$依存性が消えてしまう。

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

高次項はゼロとはならないのだが、$k$回かけたものが簡単に計算できる。

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

これにより、リュービル演算子を分解した4つの要素を指数関数の肩に乗せたものが全て厳密に評価できた。これを使ってシンプレクティック積分と同様に数値積分法を構築するのがRESPA (reference system propagator algorithm)である。特に時間反転対称にしたものをr-RESPAというが、現在では単にRESPAというと時間反転対称にしたものを指すと思われる。

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

さて、時間反転対称であるとは、時間を$h$だけ進める演算子を演算子た結果を逆に解いたら、それは時間を$-h$だけ進める演算子をかけた結果と等しい、ということである。従って、

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

これがなぜかは、演算子の形で考えると容易に理解できる。

今、厳密な時間発展演算子$U(h)$を、一次のシンプレクティック積分$\tilde{U}_1(h)$で近似したとしよう。運動項の部分Liouville演算子を$iL_K$、ポテンシャル項を$iL_V$とすると、

$$
\tilde{U}_1(h) = \mathrm{e}^{ihL_K} \mathrm{e}^{ihL_V}
$$

である。従って、この時間発展により$(p,q)$が$(P,Q)$に移ったのなら、

$$
\begin{aligned}
\begin{pmatrix}
P \\
Q
\end{pmatrix}
&= 
\tilde{U}_1(h)
\begin{pmatrix}
p \\
q
\end{pmatrix} \\
&= 
\mathrm{e}^{ihL_K} \mathrm{e}^{ihL_V}
\begin{pmatrix}
p \\
q
\end{pmatrix}
\end{aligned}
$$

これを$(p,q)$について逆に解くと、

$$
\begin{aligned}
\begin{pmatrix}
p \\
q
\end{pmatrix}
&= 
\mathrm{e}^{-ihL_V} \mathrm{e}^{-ihL_K}
\begin{pmatrix}
P \\
Q
\end{pmatrix}\\
&\equiv
\tilde{U}_1^{-1}(h)
\begin{pmatrix}
P \\
Q
\end{pmatrix}\\
\end{aligned}
$$

要するに、指数分解された項の左右が入れ替わってしまうため、$\tilde{U}_1^{-1}(h) = \tilde{U}_1(-h)$が成り立たない。

それに対して、二次の対称分解で構築した時間発展演算子$\tilde{U}_2(h)$は

$$
\tilde{U}_2(h) = \mathrm{e}^{ihL_K/2} \mathrm{e}^{ihL_V} \mathrm{e}^{ihL_K/2}
$$

と左右対称な形をしているため、右辺から左辺に移項しても形が変わらず、それにより$\tilde{U}_2^{-1}(h) = \tilde{U}_2(-h)$が成立していることがわかる。

全く同様にして、温度制御がある系にRESPAを適用して構築した時間発展演算子、

$$
\tilde{U}(h)_\mathrm{RESPA} = \mathrm{e}^{ihL_K/2} \mathrm{e}^{ihL_T/2} \mathrm{e}^{ihL_V} \mathrm{e}^{ihL_T/2} \mathrm{e}^{ihL_K/2}
$$

も左右対称な形になっているため、

$$
\tilde{U}(h)^{-1}_{\mathrm{RESPA}} = \tilde{U}(-h)_\mathrm{RESPA}
$$

が成立し、時間反転対称であることがわかる。

しかし、ハミルトンダイナミクスの場合には、時間反転対称でない一次の公式であっても位相空間の体積が保存し、それに伴ってエネルギーが保存するというメリットがあったが、温度制御された系においては時間反転対称であることが数値計算上なにかメリットがあるかどうかは定かではない。
