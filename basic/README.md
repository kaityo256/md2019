# 運動方程式とリュービル演算子のエルミート性

ハミルトンの運動方程式が、位相空間にどのような「流れ」を作り、それがどのような意味を持つかを説明する。

## 運動方程式とは何か

簡単のため、一次元系を考えよう。単位質量の物体が位置$x$にいて、力$f$が働いている。すると、運動方程式は

$$
\ddot{x} = f
$$

である。これは時間に関する二階微分方程式になっているが、$\dot{x} = v$を導入し、一階の連立微分方程式にしたほうが後々都合が良いのでそうしよう。

$$
\begin{aligned}
\dot{x} &= v\\
\dot{v} &= f
\end{aligned}
$$

さて、この式の意味を考えてみよう。一階の連立微分方程式にしたことで、この系の状態は$(x,v)$を指定することで一意に定まる。逆に言えば$(x,v)$がこの系を記述する空間を張っている。この空間を位相空間と呼ぶ。

> よく知っている人のための注：ここで張られた空間は、$x$で張られた状態空間$M$の接バンドル(Tangent-Bundle)$TM$になっており、ハミルトニアンを考える余接バンドル$T^*M$とは違う空間である。しかし、以下ではどちらも位相空間として扱う。

運動方程式は時間に関する常微分方程式であるから、初期条件を与えればその後の時間発展は一意に定まる。それは、位相空間で見れば、どこかにトレーサーを置くと、あとはそのトレーサーがどのように動くか決まる、という意味となる。

例えば、原点を中心とする調和振動子系を考えよう。力$f$はポテンシャル力であり、バネ定数を1とすると$f(x) = -x$となる。すると、運動方程式は

$$
\begin{aligned}
\dot{x} &= v\\
\dot{v} &= -x
\end{aligned}
$$

となる。図に書いてみるとわかるが、これは原点を中心とする、反時計周りのベクトル場を作っている。そういう意味において、運動方程式は位相空間に速度場を作る。

さて、講義では詳細に触れないが、ラグランジアン$L$を導入し、その時間積分を作用と呼ぶとき、その作用を最小にせよ、という要請(最小作用の原理)からオイラー・ラグランジュ方程式が導かれるのであった。さらに変数をからルジャンドル変換をすることによりハミルトニアンが導かれ、より対称性の高い運動方程式、ハミルトンの運動方程式が導かれる。以下では、ルジャンドル変換後の一般化座標$q_i$と一般化運動量$p_i$が張る位相空間を考える。

いま、3次元$N$粒子系を考える。すると座標$q_i$,運動量$q_i$がそれぞれ$3N$個あるため、位相空間$\vec{\Gamma} = (q_1, q_2, \cdots, q_{3N}, p_1, p_2, \cdots, p_{3N})$は$6N$次元となる。この空間に
速度場$\vec{\Gamma}$を定めるには、それぞれの成分$\dot{q}_1, \cdots, \dot{p}_{3N}$全てを$\{p_i, q_i\}$の関数として指定しなければならないから、合計$6N$個の関数が必要になる。しかし、変分原理を使うと、その$6N$個の関数がたった一つのスカラー関数から決まるのであった。そのスカラー関数をハミルトニアン$H$と呼び、速度場は

$$
\begin{aligned}
\dot{p}_i &= -\frac{\partial H}{\partial q_i}\\
\dot{q}_i &= \frac{\partial H}{\partial p_i}
\end{aligned}
$$

と決まる。これをハミルトンの運動方程式と呼ぶ。繰り返しになるが、位相空間に一つスカラー関数を定めるだけで、全空間に流れ場を定めている。そう考えると不思議な気がしてくる。

さて、この式を見るとただちにわかることがいくつかある。それは$H$が時間不変量になることだ(ただし$H$は時間に陽に依存しないとする)。

$H$の時間微分を計算してみると、

$$
\begin{aligned}
\dot{H} &= \sum_i 
\left(
\frac{\partial H}{\partial q} \dot{q}
+ \frac{\partial H}{\partial p} \dot{p}
\right) \\
&= \sum_i 
\left(
\frac{\partial H}{\partial q} 
\frac{\partial H}{\partial p} 
-
\frac{\partial H}{\partial p} 
\frac{\partial H}{\partial q} 
\right)
= 0
\end{aligned}
$$

こうして$H$が時間不変量になることがわかる。これはエネルギー保存則に対応している。これはよく知られているであろう。

実は、もう一つ重要な性質がある。それは

$$
\sum_i 
\left(
    \frac{\partial \dot{p}_i}{\partial p_i}
    + \frac{\partial \dot{q}_i}{\partial q_i}
    \right) = 0
$$

が成り立つことだ。以下では、この意味を考えてみたい。

いま、運動量と座標をまとめて$z_i$という$6N$個の座標で表現しよう。先程の式は

$$
\sum_i 
    \frac{\partial \dot{z}_i}{\partial z_i}= 0
$$

となる。運動方程式とは位相空間に速度場を定めるものであった。そこで$\dot{z}_i \equiv v_i$として、これを速度場と呼ぼう。そうして先程の式をもう一度眺めると、速度場の各成分を、座標の各成分で偏微分したものの和であるから、これは速度場の発散(divergence)に他ならない。つまり、

$$
\begin{aligned}
\sum_i  \frac{\partial \dot{z}_i}{\partial z_i} &= 
\sum_i  \frac{\partial v_i}{\partial z_i} \\
&= \nabla \cdot \vec{v}\\
&= 0
\end{aligned}
$$

これは、ハミルトニアンが作る流れ場は非圧縮流である、ということを意味する。

この事実を詳しく見るために、この位相空間に分布関数$f$を定義しよう。これは流体力学の言葉で言えば流体の密度場にあたる。いま、速度場$\vec{v}$を持つ流れがある時、それに密度場をかけたものが流れ場$\vec{J} = \vec{v}f$である。確率の保存則から、流れ場と密度場は以下の連続の式を満たさなければならない。

$$
\frac{\partial f}{\partial t} + \nabla \cdot \vec{J} = 0
$$

これを計算してみよう。

$$
\begin{aligned}
\frac{\partial f}{\partial t} + \nabla \cdot \vec{J} &= 
\frac{\partial f}{\partial t} + \nabla \cdot (\vec{v}f)\\
&=  \underbrace{\frac{\partial f}{\partial t} + \vec{v} \nabla f}_{\frac{D f}{Dt}} + \underbrace{\nabla \cdot \vec{v}}_{=0} f \\
&= 0
\end{aligned}
$$

つまり、

$$
\frac{Df}{Dt} = 0
$$

となった。ここで$D/Dt$はラグランジュ微分であり、流れに沿って見た物理量の変化である。流れに沿って見ると密度が変化しない、と言っているのだから、これは非圧縮流であることを意味する。

## リュービル演算子と時間発展演算子

ハミルトンの運動方程式は位相空間に非圧縮流を作ることがわかった。それが、演算子の言葉ではどう見えるかを見てみよう。そのために、時間発展演算子とリュービル演算子を定義する。

以下簡単のため一自由度系を考える。ハミルトニアン$H(p,q)$で記述される系の運動を考えよう。ここで、ハミルトニアンは時間に陽に依存しないものとする。運動方程式は以下のように記述される。

$$
\begin{aligned}
\dot{p} &= -\frac{\partial H}{\partial q} \\
\dot{q} &= \frac{\partial H}{\partial p}
\end{aligned}
$$

この運動方程式に従い、座標$(p,q)$が変化していく。さて、この系に物理量$A(p(t),q(t))$が定義されているとしよう。$A$は$p,q$のみの関数であり、以下$A(p(t),q(t))$を$A(t)$と略記する。系が時間発展するにつれて、物理量$A$も時間変化する。$A$を用いて、この系の時間発展演算子$U(h)$を、以下のように定義する。

$$
A(t+h) = U(h) A(t)
$$

ここで、$A(t+h)$を$t$のまわりでテイラー展開しよう。

$$
\begin{aligned}
A(t+h) &= A(t) + h \frac{dA}{dt} + \frac{h^2}{2} \frac{d^2A}{dt^2} + \cdots \\
&= \sum_{k=0} \frac{h^k}{k!} \frac{d^k}{dt^k} A\\
&= \underbrace{\exp\left(h \frac{d}{dt}\right)}_{U(h)} A \\
&= U(h) A(t)
\end{aligned}
$$

ここから、時間発展演算子$U(h)$は以下のように書けることがわかった。

$$
U(h) = \exp\left(h \frac{d}{dt}\right)
$$

さて、$A$は$p,q$にのみ依存する関数であったから、その時間微分は$p$や$q$の偏微分として表現できる。

$$
\begin{aligned}
\frac{dA}{dt} &= 
\frac{\partial A}{\partial q}\dot{q}
+\frac{\partial A}{\partial q}\dot{q} \\
&= 
\underbrace{
\left(
\dot{q} \frac{\partial}{\partial q} +
\dot{p} \frac{\partial}{\partial p}
\right)
}_{-iL}
A \\
&= -iLA
\end{aligned}
$$

ここに表れた$L$をリュービル演算子と呼ぶ。$-i$をつけるのは、リュービル演算子をエルミートにするためである。

先程の式は任意の物理量$A$で成り立つから、形式的に

$$
\frac{d}{dt} = -iL
$$

となっている。従って、リュービル演算子は時間微分を与える演算子であることがわかる。先程のテイラー展開の式に代入すると、

$$
A(t+h) = \exp(-i hL) A(t)
$$

以上から、時間発展演算子$U$は、時間微分演算子$iL$を指数の肩に乗せたものであることがわかった。

リュービル演算子についてもう少し見てみよう。リュービル演算子は

$$
-i L = \dot{q} \frac{\partial}{\partial q} +
\dot{p} \frac{\partial}{\partial p}
$$

と書けた。$p,q$はハミルトンの運動方程式に従うため、$\dot{p}, \dot{q}$をハミルトニアンを用いて書き直すと、

$$
-i L = \frac{\partial H}{\partial p} \frac{\partial}{\partial q} - \frac{\partial H}{\partial q} \frac{\partial}{\partial p}
$$

と書ける。ここで、ハミルトニアンが時間に陽に依存しないため、リュービル演算子も時間に陽に依存しないことに注意。リュービル演算子を用いると、$p,q$で記述される任意の量の時間微分を表現できる。もちろん$p,q$自身の時間微分も表現できるので、運動方程式は以下のように記述できる。

$$
\frac{d}{dt}
\begin{pmatrix}
p \\
q
\end{pmatrix}
=
-i L
\begin{pmatrix}
p \\
q
\end{pmatrix}
$$

両辺を形式的に積分してみよう。$P=p(t+h)$,$Q=q(t+h)$と表記すると、

$$
\begin{pmatrix}
P \\
Q
\end{pmatrix}
=
\exp\left(-i h L\right)
\begin{pmatrix}
p \\
q
\end{pmatrix}
$$

これは、先程得られた「時間発展演算子はリュービル演算子を指数関数の肩に乗せたもの」

$$
U(h) = \exp\left(-i h L\right)
$$

という結果と同じである。

## リュービル演算子のエルミート性

さて、やや唐突だが、演算子エルミート性について考えてみよう。何か線形空間があり、その要素$f$と$g$の組に対してスカラー量を定める写像$(f,g)$を、$f$と$g$の内積と呼ぶのであった。この要素に作用する演算子$X$がエルミートであるとは、$(Xf, g) = (f, Xg)$が成り立つことである。以下では、リュービル演算子がエルミートであることを示す。

簡単のため、1自由度系を考える。位相空間は$\vec{\Gamma}=(p,q)$で張られている。この空間に住む関数$f(p,q), g(p,q)$を考え、その内積を以下のように定める。

$$
(f,g) \equiv \int d \Gamma f^* g 
$$

ただし$f^*$は複素共役、$d \Gamma = dp dq$である。

さて、リュービル演算子がエルミートであることを示すには、$L$に対して、$(Lf,g) = (f, Lg)$を示せば良い。

リュービル演算子は

$$
- i L = \dot{q} \frac{\partial }{\partial q}
+\dot{p} \frac{\partial }{\partial p}
$$

であったから、両辺に$i$をかけて、

$$
L = i \left(
\dot{p} \frac{\partial }{\partial p}
+\dot{q} \frac{\partial }{\partial q}
\right)
$$

となる。これを$(f, Lg)$に代入すると、

$$
\begin{aligned}
(f, Lg) &= \int d\Gamma f^* Lg \\
&= \int d\Gamma f^*
i \left(
\dot{p} \frac{\partial }{\partial p}
+\dot{q} \frac{\partial }{\partial q}
\right)g
\end{aligned}
$$

ここで、部分積分をすることで、$g$にかかっている微分を$f$に移す。

$$
\begin{aligned}
(f, Lg) &= \int d\Gamma f^* Lg \\
&= -\int d\Gamma 
i
\left(
\frac{\partial}{\partial p} \left(f^*  \dot{p} \right)g
+  \frac{\partial}{\partial q} \left(f^* \dot{q} \right)g
\right) \\
&= -\int d\Gamma 
i 
\left(
\dot{p} \frac{\partial f^*}{\partial p} + \dot{q}\frac{\partial f^*}{\partial q}
+
\underbrace{
f^*\frac{\partial \dot{q}}{\partial q}
+f^*\frac{\partial \dot{p}}{\partial p}
}_{=0}
\right)g \\
&=  \int d\Gamma 
\left(
\mathcal{L^*}f^*
\right) g  \\
&= (Lf, g)
\end{aligned}
$$

ただし、途中で部分積分で出てくるゴミが境界条件でゼロになることを仮定した。以上から、

$$
(Lf, g) = (f,Lg)
$$

であるから、リュービル演算子$L$がエルミートであることが示された。

逆に、リュービル演算子がエルミートであるためには、部分積分で出てくる余計な項がゼロ、つまり

$$
\frac{\partial \dot{q}}{\partial q} +
\frac{\partial \dot{p}}{\partial p} =0
$$

でなければならないことがわかる。つまり、リュービル演算子がエルミートであると、そのリュービル演算子が作る位相空間での流れは非圧縮になり、位相空間の流れが非圧縮であれば、対応するリュービル演算子がエルミートであることがわかる。

流れが非圧縮であれば、密度場は定数である。密度場に対応するのは分布関数であり、これは流れにともなって確率密度が変化しない、すなわち「ミクロカノニカル分布」に対応していることがわかる。

以上をまとめると、以下のようになる。

* ハミルトンの運動方程式に対応するリュービル演算子はエルミート演算子になる
* エルミートであるリュービル演算子が作る流れは非圧縮流となる
* 流れが非圧縮であるから密度場は変化しない
* これはミクロカノニカル分布$f=\mathrm{const.}$を意味する

ただし、この結果からわかることは、運動方程式に伴う確率流れの密度が一定であるということであり、これがミクロカノニカル分布になるかどうかは、別途エルゴード性の議論が必要になる。
