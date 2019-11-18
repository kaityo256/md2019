# 圧力

## 圧力の定義

しつこいが、分子動力学法において規定されているのは全エネルギーのみである。それ以外の物理量は我々が定義しなければならない。ここでは、圧力について考えてみよう。

## 古典的なビリアル定理

まずは教科書によく書いてある古典的なビリアル定理から圧力を定義してみよう。天下りだが、以下のような量を考えよう。

$$
G = \sum_i \vec{p}_i \cdot \vec{q}_i
$$

両辺を時間微分してみる。

$$
\begin{aligned}
\dot{G} &= \sum_i \underbrace{\vec{p}_i \cdot \dot{\vec{q}}_i}_{(1)} + \sum_i \underbrace{\dot{\vec{p}}_i \cdot \vec{q}_i}_{(2)}
\end{aligned}
$$

まず(1)の項だが、$\dot{\vec{q}}_i = \vec{p}_i/m$であるから、

$$
\vec{p}_i \cdot \dot{\vec{q}}_i  = \frac{\vec{p}_i^2}{m}
$$

である。エネルギー等分配則から

$$
\frac{\vec{p}_i^2}{m} = 3 k_B T
$$

和を取ると、

$$
\sum_i \frac{\vec{p}_i^2}{m} = 3 Nk_B T
$$

次に、(2)だが、$\dot{\vec{p}}_i = \vec{f}_i$であるから、

$$
\dot{\vec{p}}_i \cdot \vec{q}_i = \vec{q}_i \cdot \vec{f}_i 
$$

ここで$\vec{f}_i$は、粒子$i$にかかる力である。ここで、粒子間に働く力と、外力$\vec{f}^\mathrm{ext}_i$にわけよう。粒子間力は二体力とし、粒子$j$から$i$に働く力を$\vec{f}_{ij}$とすると、

$$
\dot{\vec{p}}_i \vec{q}_i = \vec{q}_i \sum_{i\neq j} \vec{f}_{ij} + \vec{q}_i  \cdot \vec{f}^\mathrm{ext}_i
$$

と書ける。$i$に関する和を取ると、

$$
\begin{aligned}
\sum_i \dot{\vec{p}}_i \vec{q}_i &= \underbrace{\sum_i \vec{q}_i \sum_{i\neq j} \vec{f}_{ij}}_{(3)} + \underbrace{\sum_i \vec{q}_i \cdot \vec{f}^\mathrm{ext}_i}_{(4)}
\end{aligned}
$$

まずは(3)について計算しよう。

$$
\begin{aligned}
(3) &=  \sum_i \vec{q}_i \sum_{i\neq j} \vec{f}_{ij} \\
&= \sum_{i< j} \left( \vec{q}_i  \cdot \vec{f}_{ij} + \vec{q}_j  \cdot \vec{f}_{ji} \right) \\
&=  \sum_{i< j}\left( \vec{q}_i  \cdot \vec{f}_{ij} - \vec{q}_j  \cdot \vec{f}_{ij} \right) \\
&=  \sum_{i< j} \vec{q}_{ij} \cdot \vec{f}_{ij}
\end{aligned}
$$

ただし、$\vec{q}_{ij} = \vec{q}_i - \vec{q}_j$である。

次に、外力からの寄与(4)を計算する。$\vec{f}^\mathrm{ext}_i$は、外力から粒子$i$に働く力であるから、その反作用が壁に働く。従って、

$$
\begin{aligned}
(4) &= \sum_i \vec{q}_i \cdot \vec{f}^\mathrm{ext}_i \\
&= -  \int_{\partial V} \vec{r} \cdot (P\vec{n}) d A \\
&= -P \int_V \underbrace{\nabla \cdot \vec{r}}_3 dV \\
&= -3PV
\end{aligned}
$$

ただし、途中でガウスの定理を使った。また、

$$
\nabla \cdot \vec{r} = \frac{\partial x}{\partial x} +
\frac{\partial y}{\partial y} +
\frac{\partial z}{\partial z}  = 3
$$

を用いた。

以上から、

$$
\dot{G} = 3Nk_B T + \sum_{i< j} \vec{q}_{ij} \cdot \vec{f}_{ij} -3PV
$$

となった。

次に、両辺の時間平均を取ろう。時間変化する物理量$A(t)$の時間平均$\bar{A}$は、

$$
\bar{A} \equiv \lim_{\tau \rightarrow \infty} \frac{1}{\tau}
\int_0^\tau dt A(t)
$$

で定義される。$\dot{G}$の時間平均は

$$
\begin{aligned}
\bar{\dot{G}} &= \lim_{\tau \rightarrow \infty} \frac{1}{\tau}
\int_0^\tau dt \dot{G}(t)\\
&= \lim_{\tau \rightarrow \infty} \frac{G(\tau) - G(0)}{\tau}
\end{aligned}
$$

ここで、$G$という量に上限があるなら、十分に長い$\tau$を取ればゼロになる。以上から、

$$
\bar{\dot{G}} = 3Nk_B T + \overline{\sum_{i< j} \vec{q}_{ij} \cdot \vec{f}_{ij}} -3PV = 0
$$

時間平均$\bar{A}$とアンサンブル平均$\left< A\right>$を同一視すると、

$$
PV = N k_BT + \frac{1}{3} \left< \sum_{i< j} \vec{q}_{ij} \cdot \vec{f}_{ij}\right>
$$

これが古典的なビリアル定理による、粒子系の圧力の導出である。最初の$N k_B T$が理想気体からの寄与、つまり粒子の運動から圧力への寄与であり、右辺第二項が相互作用からの寄与である。

## 分配関数からの導出

さて、先ほどのビリアル定理からの圧力の導出では、

* いきなり$G$という量を定義して時間微分する意味が不明瞭
* 粒子にかかる力を粒子間力と外力に分けたが、周期境界条件ではどうなるかわかりにくい

という問題がある。そこで、「分子動力学法の世界にはハミルトニアンしか存在しない」という立場から、分配関数経由で導出してみよう。

まずはこの系のハミルトニアンを定義しよう。

$$
H = \sum_i \frac{\vec{p}_i^2}{2m} + \sum_{i < j} \Phi(q_{ij})
$$

これがこの系の全エネルギーを定義する。

さて、熱力学関係式

$$
P = - \left(\frac{\partial F}{\partial V} \right)_T
$$

から出発しよう。$P$は圧力、$T$が温度、$F$はヘルムホルツ自由エネルギーである。

ヘルムホルツ自由エネルギーは分配関数$Z$を用いて

$$
F = - k_B T \ln Z
$$

と書けるので、

$$
\begin{aligned}
P &= - \left(\frac{\partial F}{\partial V} \right)_T\\
&= k_B T \frac{\partial\ln Z }{\partial V}
= k_B T \frac{1}{Z} \frac{\partial Z}{\partial V}
\end{aligned}
$$

となる。さて、分配関数はハミルトニアン$H$を用いて

$$
Z = \int d \Gamma \mathrm{e}^{-\beta H}
$$

と書けるから、

$$
\frac{\partial Z}{\partial V} = \int d \Gamma \left(-\beta \frac{\partial H}{\partial V} \right) \mathrm{e}^{-\beta H}
$$

である。これを先程の式に代入して整理すると、

$$
\begin{aligned}
P &= - Z^{-1} \int d \Gamma \left(\frac{\partial H}{\partial V} \right) \mathrm{e}^{-\beta H} \\
&= -\left< \frac{\partial H}{\partial V} \right>
\end{aligned}
$$

すなわち、圧力を求めるにはハミルトニアンの体積微分のアンサンブル平均を取れば良い。
もちろんこの式は、熱力学関係式

$$
P = - \left(\frac{\partial U}{\partial V} \right)_S
$$

に対応している。

さて、ハミルトニアンの体積微分を取るために、系のサイズを変化させることを考えよう。一辺$L$の立方体の系を考えると、$V=L^3$である。この系の長さを一様に$L \rightarrow \alpha L$と拡大することを考える。

拡大された世界のハミルトニアンは

$$
H(\alpha) = \sum_i \frac{\alpha^2 \vec{p}_i^2}{2m} + \sum_{i < j} \Phi(\alpha q_{ij})
$$

となる。

やや紛らわしいが、$\alpha = 1$の時の体積を$V$とすると、$V \rightarrow \alpha^3 V$となるから

$$
dV = 3 \alpha^2 V d \alpha
$$

となる(微分する変数としての$V$と基準体積$V$で同じ記号を使っていることに注意)。これにより、

$$
\frac{\partial H}{\partial V} = \lim_{\alpha \rightarrow 1} \frac{\partial H}{\partial \alpha} \frac{d \alpha}{dV} = \lim_{\alpha \rightarrow 1} \frac{\partial H}{\partial \alpha} \frac{1}{3 \alpha^2 V}
$$

と、$V$による微分を$\alpha$による微分に置き換えることができる。

以下、$\partial H/\partial \alpha$を計算しよう。

まず、運動エネルギー部分を考える。空間を一辺$\alpha$倍にすると、$\vec{q}_i \rightarrow \alpha \vec{q}_i$となる。ここで、

$$
\vec{p}_i = \frac{\partial L}{\partial \dot{\vec{q}}_i}
$$

であったから、$\vec{p}_i \rightarrow \vec{p}_i / \alpha$となることに注意すると、

$$
\begin{aligned}
\lim_{\alpha \rightarrow 1}\frac{\partial K(\alpha)}{\partial \alpha} &=
\lim_{\alpha \rightarrow 1} \sum_i \frac{\partial}{\partial \alpha } \frac{\vec{p}_i^2}{2m \alpha^2} \\
&= -2K \\
&= -3 N k_B T
\end{aligned}
$$

ただし、途中でエネルギー等分配則を用いた。

次にポテンシャル部分の計算であるが、まずポテンシャル項の$\alpha$微分は

$$
\lim_{\alpha \rightarrow 1}  \frac{\partial\Phi(\alpha q_{ij})}{\partial \alpha} 
= \Phi'(q_{ij}) q_{ij} $$

ここで、

$$
q_{ij} = |\vec{q}_i - \vec{q}_j|
$$

である。粒子$j$から$i$に働く力$\vec{f}_{ij}$は

$$
\vec{f}_{ij} = - \Phi'(q_{ij}) \frac{\vec{q}_i - \vec{q}_j}{|\vec{q}_i - \vec{q}_j|}$$

である(ベクトルの向きに注意)。以上をまとめると、

$$
\lim_{\alpha \rightarrow 1}  \frac{\partial \Phi(\alpha q_{ij})}{\partial \alpha}  = -\vec{q}_{ij} \cdot \vec{f}_{ij}
$$

となるので、両辺和を取れば、

$$
\lim_{\alpha \rightarrow 1} \sum_{i< j} \frac{\partial \Phi(\alpha q_{ij})}{\partial \alpha} = -\sum_{i < j} \vec{q}_{ij} \cdot \vec{f}_{ij}\\
$$


以上から、

$$
\begin{aligned}
-P &= \left<\frac{\partial H}{\partial V} \right>_S\\
&= \lim_{\alpha \rightarrow 1} \left<\frac{\partial H}{\partial \alpha}\right> \underbrace{\frac{d \alpha}{d V}}_{1/3V} \\
&= \frac{1}{3V} \left(- 3 Nk_B T -  \right<\sum_{i < j} \vec{q}_{ij} \cdot \vec{f}_{ij} \left>\right)
\end{aligned}
$$

整理すると、

$$
PV = N k_BT + \frac{1}{3} \left< \sum_{i< j} \vec{q}_{ij} \cdot \vec{f}_{ij}\right>
$$

先程ビリアル定理で導かれた圧力が導出された。導出を見れば、分母の$3V$は$d\alpha /d V$から、$Nk_B T$の項は運動エネルギー由来、ビリアル項は相互作用由来であることがわかるであろう。また、境界条件に依存しない導出であることもわかるであろう。

分子動力学法では$\vec{q}_{ij} \cdot \vec{f}_{ij}$は容易に計算できるため、これで圧力が計算できることになる。

## 局所圧力

