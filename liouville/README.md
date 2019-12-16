# 8. Generalized Liouville's Theorem of non-Hamiltonian systems

ハミルトンの運動方程式に付随するLiouville演算子は必ずエルミートとなり、エルミートなリュービル演算子が位相空間に作る「流れ」は非圧縮となることを見た。分布関数は位相空間における流れの「密度場」に対応するため、非圧縮な流れ場では密度は一定、すなわち、分布関数が不変に保たれることがわかる(Liouvilleの定理)。

Liouvilleの定理は、時間発展に伴って位相空間体積が変化しないことを主張する。例えば一自由度系$(p,q)$において、時間$h$後に$(P,Q)$に系が変化した時、その前後で位相空間体積が変化しないとは、

$$
dPdQ = dpdq
$$

が成立することを指す。さて、$P$も$Q$も$p, q$の関数であるから、

$$
dP dQ = 
\underbrace{
\begin{pmatrix}
\frac{\partial P}{\partial p} & \frac{\partial P}{\partial q} \\
\frac{\partial Q}{\partial p} & \frac{\partial Q}{\partial q} 
\end{pmatrix}
}_J
dpdq
$$

である。ただし、$J$は$(p,q)$から$(P,Q)$への変換のヤコビアンである。$dPdQ = dpdq$が成り立つのであるから、Liouvilleの定理は力学系に付随するヤコビアンが常に1となることを主張する。

では、熱浴をつけ、定常状態がカノニカル分布となる場合はどうだろうか？当然、分布関数が時間変化するため、位相空間体積も一定ではなく、ヤコビアンも時間変化する。この時、Liouvilleの定理と同様なことが言えないだろうか？

以下では、一般化された形でのLiouvilleの定理を導き、温度制御された系においてLiouville演算子がどのような性質を持つか議論する。

## 8.1 Jacobi's Formula

まず、位相空間$(p_i, q_i)$をまとめて$\vec{z} = (p_1, p_2, \cdots, p_{3N}, q_1, q_2, \cdots, q_{3N})=(z_1, \cdots, z_{6N})$と記述しよう。この系が時刻$t$から$t+h$まで時間発展発展し、$\vec{z}(t)$から$\vec{z}(t+h)$になったとしよう。

$$
\vec{z}(t) \xrightarrow{h} \vec{z}(t+h)
$$

時間発展は決定論的であるから、$\vec{z}(t+h)$は$\vec{z}(t)$を決めれば一意に決まる。従って、$\vec{z}(t+h)$は$\vec{z}(t)$の関数となるから、両者の間にヤコビアンを考えることができる。

$$
J(t+h|t) \equiv \frac{\partial \vec{z}(t+h)}{\partial \vec{z}(t)}
$$

我々は一般の力学系におけるヤコビアンの時間発展を知りたいので、ヤコビアンが従う方程式を書き下したい。すなわち、$dJ/dt$が知りたい。位相空間の流れ$\dot{\vec{z}}$が指定されれば力学系が定まるのであるから、$dJ/dt$と$\dot{\vec{z}}$の間になんらかの関係があるはずである。実際、以下の等式が成り立つ。

$$
\frac{dJ}{dt} = J \nabla \frac{d \vec{z}}{dt}
$$

上記を証明するために、まずJacobiの公式と呼ばれる式を証明する。

$A$を時間依存する正方行列とする時、以下が成り立つ。

$$
\frac{d}{dt}\det A = \det A \mathrm{Tr} \left(A^{-1} \frac{d A}{dt} \right)
$$

これを証明するため、まずは$\det A$の時間微分を定義どおりに計算し、$A(t+h)$を一次までテイラー展開しよう。

$$
\begin{aligned}
\frac{d}{dt}\det A &= \lim_{h\rightarrow 0} \frac{\det A(t+h) - \det A(t)}{h} \\
&= \lim_{h\rightarrow 0} \frac{\det \left(A + \frac{dA}{d t} h\right) - \det A(t)}{h}
\end{aligned}
$$

さて、行列の積の行列式は、行列式の積でかける。

$$
\det (XY) = \det X \det Y
$$

したがって、以下の恒等式が成り立つ。

$$
\begin{aligned}
\det X &= \det A A^{-1} X \\
&= \det A \det A^{-1}X
\end{aligned}
$$

この恒等式を使って、先程の四季から$\det A$をくくりだすことを考える。

$$
\begin{aligned}
\frac{d}{dt}\det A &= \lim_{h\rightarrow 0} \frac{\det \left(A + \frac{dA}{d t} h\right) - \det A(t)}{h} \\
&= \lim_{h\rightarrow 0} \frac{\det \left[A \left(I + A^{-1}\frac{dA}{d t} h\right)\right] - \det (A(t)I)}{h} \\
&= \lim_{h\rightarrow 0} \frac{\det A \det\left(I + A^{-1}\frac{dA}{d t} h\right) - \det A \det I}{h} \\
&= \det A \lim_{h\rightarrow 0}\frac{ \det\left(I + A^{-1}\frac{dA}{d t} h\right) - \det I}{h}
\end{aligned}
$$

ここで、行列の性質として以下が成り立つ。

$$
\det(I + h X) = 1 + h \mathrm{Tr} X + O(h^2)
$$

これを$X = A^{-1} dA/dt$として先の式に代入すると、

$$
\frac{d}{dt}\det A = \det A \mathrm{Tr} \left(A^{-1} \frac{d A}{dt} \right)
$$

これが求めたい式であった。

## 8.2 Dynamics of Jacobian

さて、位相空間$\vec{z}$上に、ダイナミクス$\vec{z}$が与えられているとしよう。時刻$0$から時刻$t$への時間発展

$$
\vec{z}(0) \xrightarrow{h} \vec{z}(t)
$$

を考える(先程は$t$から$t+h$を考えたが、今後時間微分が煩雑にならないように時刻$0$からのダイナミクスを考えている)。

この時間発展のヤコビ行列$M_{ij}$を以下のように定義する。

$$
M_{ij} = \frac{\partial z_t^j}{\partial z_0^i}
$$

ただし、$\vec{z}(0)$の$i$成分を$z_0^i$、$\vec{z}(t)$の$i$成分を$z_t^i$と表記している。

ヤコビアンは$J=\det M$であるから、先程のヤコビの公式により、

$$
\begin{aligned}
\frac{dJ}{dt} &= \frac{d}{dt} \det M \\
&= \underbrace{\det M}_J \mathrm{Tr} \left(M^{-1} \frac{dM}{dt}\right) \\
&= J \underbrace{\mathrm{Tr} \left(M^{-1} \frac{dM}{dt}\right)}_{(*)}
\end{aligned}
$$

以下、(*)の部分の計算を実行しよう。まず、ヤコビ行列の性質として、逆行列は、各要素の偏微分をひっくり返したものになる。

$$
\left(M\right)^{-1}_{ij} = \frac{\partial z_0^j}{\partial z_t^i}
$$

また、ヤコビ行列の時間依存性は、$z_t^j$にしかないので、

$$
\left(\frac{d M}{dt}\right)_{ij} = \frac{\partial \dot{z_t^j}}{\partial z_0^i}
$$

以上から、

$$
\left(M^{-1} \frac{dM}{dt}\right)_{ij} = \sum_k \frac{\partial z_0^k}{\partial z_t^i} \frac{\partial \dot{z_t^j}}{\partial z_0^k}
$$

この対角和を取ると、

$$
\mathrm{Tr} \left(M^{-1} \frac{dM}{dt}\right) = \sum_i \sum_k \frac{\partial z_0^k}{\partial z_t^i} \frac{\partial \dot{z_t^i}}{\partial z_0^k}
$$

さて、ライプニッツ則により、

$$
\frac{\partial \dot{z_t^i}}{\partial z_0^k} = \sum_j \frac{\partial \dot{z_t^i}}{\partial z_t^j}
\frac{\partial z_t^j}{\partial z_0^i}
$$

これを先程の式に代入して、

$$
\mathrm{Tr} \left(M^{-1} \frac{dM}{dt}\right) = \sum_i \sum_j \sum_k 
\underbrace{\frac{\partial z_0^k}{\partial z_t^i}}_{M^{-1}}
\underbrace{\frac{\partial \dot{z_t^i}}{\partial z_t^j}}_{X}
\underbrace{\frac{\partial z_t^j}{\partial z_0^i}}_{M}
$$

この右辺は$\mathrm{Tr}(M^{-1}XM)$の形になっている。対角和の性質から、

$$
\mathrm{Tr} (M^{-1}XM) = \mathrm{Tr} (MM^{-1}X) = \mathrm{Tr}X
$$

ここで行列$X$の要素は

$$
X_{ij} = \frac{\partial \dot{z_t^i}}{\partial z_t^j}
$$

この行列の対角和を取ると、

$$
\begin{aligned}
\mathrm{Tr} X &= \sum_i X_{ii} \\
& = \sum_i \frac{\partial \dot{z_t^i}}{\partial z_t^i}\\
& = \nabla \cdot \dot
{\vec{z}}
\end{aligned}
$$

これは、ダイナミクス$\dot{\vec{z}}$の発散に他ならない。

以上から、時間発展のヤコビアンの時間微分は

$$
\frac{dJ}{dt} = J  \nabla \cdot  \dot{\vec{z}}
$$

となる。

## 8.3 Generalized Liouvllie's Theorem

今、位相空間が$\vec{z}$で張られており、ダイナミクス$\dot{\vec{z}}$が与えられているとする。そして、時刻$t=0$で$\vec{z}(0)$であった状態が、時間発展により時刻$t$で$\vec{z}(t)$になったとしよう。

この系がハミルトンダイナミクスであった場合は、時間発展の前後で位相空間体積が保存する、すなわち

$$
d\vec{z}(t) = d\vec{z}(0)
$$

が成り立ち、これがLiouvilleの定理であった。

さて、一般のダイナミクスでは位相空間体積が変化する。その変化をヤコビアンで表現しよう。

位相空間の速度場の発散(divergence)を$\kappa$と表現しよう。

$$
\kappa \equiv \nabla \cdot  \dot{\vec{z}}
$$

先程示したヤコビアンの時間発展の式から

$$
\frac{dJ}{dt} = \kappa J 
$$

形式的に積分すると

$$
J(t|0) = \int_0^t \kappa dt
$$

$\dot{\omega} = \kappa$という変数を導入すると、

$$
J(t|0) = \exp\left(\omega(t) - \omega(0)  \right)
$$

さて、ヤコビアンの定義から

$$
d \vec{z}(t) = J d \vec{z}(0)
$$

先程求めたヤコビアンの式を使うと、

$$
d \vec{z}(t) =  \exp\left(\omega(t) - \omega(0)  \right) d \vec{z}(0)
$$

整理して、

$$
\mathrm{e}^{-w(t)}d \vec{z}(t) = \mathrm{e}^{-w(0)} d \vec{z}(0)
$$

上式が時刻に寄らず成立する。これを一般化Liouvilleの定理(Generalized Liouville's Theorem)と呼ぶ。

これは、$\mathrm{e}^{-w}d\vec{z}$という量が不変計量(invariant measure)になっていることを意味する。

ハミルトンダイナミクスでは、速度場が非圧縮流、すなわち

$$
\kappa \equiv \nabla \cdot \dot{\vec{z}} = 0
$$

になるのであった。$\dot{\omega} = \kappa$であったから、$\dot{\omega} = 0$、すなわち$\omega$は時間に依存しない定数となる。ここから、$w(t) = w(0)$であるから、

$$
d \vec{z}(t) =  d \vec{z}(0)
$$

通常のLiouvilleの定理が導かれた。先程の一般化Liouvilleの定理が、通常のLiouvilleの定理の自然な一般化になっていることがわかるであろう。

* M. E. Tuckerman et al. "On the classical statistical mechanics of non-Hamiltonian systems", Europhys. Lett. 45, 149 (1999).
