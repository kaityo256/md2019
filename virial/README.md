# Generalized Virial Theorem and Microcanonical Temperature

これまで、Nose-Hoover熱浴やLangevin熱浴で、定常状態が指定の温度のカノニカル分布となることを見てきた。分子動力学法においては、運動エネルギー$K$の定数倍を温度と定義していた。3次元$N$粒子系なら

$$
T = \frac{2K}{3N k_B}
$$

を温度と呼ぶのであった。念のためこれをもう一度導いておこう。

まず、ビリアルと呼ばれる以下の量を定義する。

$$
I = \sum_{i=1}^{3N} p_i \frac{\partial H}{\partial p_i}
$$

なお、以下では次元の自由度は$i$に含める。すなわち、和は$1$から$3N$まで取ることにする。ハミルトニアンが

$$
H = \sum_i \frac{p_i^2}{2m} + V(\{q_i\})
$$

という形をしている時

$$
I = \sum_i 
p_i \frac{\partial H}{\partial p_i} = \sum_i \frac{p_i^2}{m} = 2 K
$$

が成り立つ。これでビリアルが運動エネルギー$K$と結びついた。

また、カノニカル分布が

$$
f = Z^{-1} \exp(-\beta H)
$$

という形をしているとき、物理量$A$のアンサンブル平均は

$$
\left< A \right> \equiv \int A f d\Gamma 
$$

で定義されるのであった。ここで、部分積分から容易に

$$
\left<p_i \frac{\partial H}{\partial p_i} \right> = k_B T
$$

であることから、

$$
\left< I \right> = \sum_i \left<p_i \frac{\partial H}{\partial p_i} \right>  = 3 N k_B T
$$

以上の二式から$I$を消去すると、

$$
2K = 3 N k_B T
$$

こうして運動エネルギーと温度が結びついた。

本質的にはカノニカル分布と部分積分しか使っていないので、$p_i$ではなく$q_i$で微分して

$$
\left< q_i \frac{\partial H}{\partial q_i}\right> = k_B T
$$

も成立し、こちらを状態温度と呼ぶのであった。

さて、上記をさらに一般化しよう。位相空間$\Gamma = (p_i, q_i)$をまとめて、$\vec{z} = (p_1, p_2, \cdots, p_{3N}, q_1, q_2, \cdots, q_{3N})$で表現しよう。

今、位相空間上に、ベクトル場$\vec{B}$を定義する。$\vec{B}$は$3N$次元のベクトルで、その成分$B_i$それぞれについて以下が成り立つことを示すのは容易であろう。

$$
\left<B_i \frac{\partial H}{\partial z_i}\right> = k_B T \left<\frac{\partial B_i}{\partial z_i} \right>
$$

これをベクトル場の言葉で書けば、

$$
\left<\vec{B} \cdot \nabla H \right>
= k_B T \left< \nabla \cdot \vec{B}\right>
$$

この左辺を一般化ビリアル、この式を一般化ビリアル定理(Generalized Virial Theorem)と呼ぶ。
