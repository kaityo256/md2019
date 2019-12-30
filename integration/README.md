# 4 Numerical Integration

## 4.1 Integration of ODE

何か物理系が、変数の組$\vec{x} = (x_1, x_2, \cdots)$で記述されているとしよう。この変数が時間に依存しており、その時間微分$\dot{\vec{x}}$が、$\vec{x}$の関数として書かれている時、この系を *力学系(Dynamical System)* と呼ぶ。特に、外力がなく、系を記述する変数だけで方程式が閉じている時、この系を *自励的(autonomous)* と言う。以下、自励的な系のみを扱う。

我々は、時刻$t=0$における系の状態$\vec{x}(0)$を指定した時、任意の時刻$t$の系の状態$\vec{x}(t)$を知りたい。このような問題設定を初期値問題と呼び、数値計算における基本課題となっている。以下、簡単のために一次元系で考えよう。方程式は以下のように表現されている。

$$
\frac{dx }{dt} = f(x)
$$

これを形式的に積分すると、任意の時刻$t$における$x$の値は以下のように求めることができる。

$$
x(t) = x(0) + \int_0^t f(x) dt
$$

さて、右辺の積分は一般には厳密に求積することはできない。そこでなんらかの近似を行うことにしよう。最も簡単な近似は、短い時間$h$の間であれば$f(t)$はあまり変化しないと思って、定数として扱ってしまうことであろう。

$$
x(t+h) = x(t) + \int_t^{t+h} f(x) dt \sim x(t) + f(x)h
$$

こうして、$x(t)$から$x(t+h)$が求まるので、$t=0$から$h$刻みで次々と状態を更新していけば、任意の時刻の$x(t)$を求めることができる。このような方法を数値積分と呼び、特に今回の積分方法はオイラー法と呼ばれる。

このオイラー法の精度を調べてみよう。$x(t+h)$を$t$のまわりでテイラー展開してみよう。

$$
x(t + h) = x(t) + \dot{x}(t) h + O(h^2)
$$

ここで、$\dot{x} = f(x)$であったから、

$$
x(t + h) = x(t) + f(x) h + O(h^2)
$$

オイラー法は、$h$に関して一次まで正しい。これを一次精度であると呼ぶ。さて、オイラー法は非常に精度が悪く、数値解が指数関数的に厳密解から離れていくことが知られている。

さて、オイラー法は、本来は時刻$t$から$t+h$まで時々刻々と変化する微分係数$f(x)$を、時刻$t$での値で代表させたのだが、さすがにこれは乱暴に過ぎた。そこで、時間刻み$h$ではなく、その中点$h/2$での時間微分を使うことにしよう。

まず、系を一次のオイラー法で$h/2$だけ時刻を進める。すると、その場所$x_m$は

$$
x_m = x(t) + \frac{f(x)h}{2}
$$

となる。この地点での微分係数$f(x_m)$を使って、あらためて現在地点$x(t)$から$x(t+h)$の位置を推定すると、

$$
x(t+h) = x(t) + f(x_m) h
$$

となる。これは中点法と呼ばれ、二次精度となる。念の為に確認しておこう。$x(t+h)$を二次までテイラー展開しよう。

$$
\begin{aligned}
\ddot{x} &= \frac{df}{dt} \\
&=\frac{df}{dx} \frac{dx}{dt}  \\
&= f'(x) f(x)
\end{aligned}
$$

であることに注意すると、

$$
\begin{aligned}
x(t+h) &= x(t) + \dot{x}(t) h + \frac{\ddot{x}h^2}{2} +O(h^3)\\
&= x(t) + f(x) h + \frac{f'(x)f(x)h^2}{2} +O(h^3)
\end{aligned}
$$

となる。また、中点法は、

$$
\begin{aligned}
x(t+h) &= x(t) + f(x_m) h \\
&= x(t) + f\left( x(t) + \frac{f(x)h}{2} \right)h \\
&= x(t) + f(x)h + \frac{f'(x) f(x)h^2}{2} + O(h^3)
\end{aligned}
$$
となり、$x(t+h)$のテイラー展開と$h$の二次まで一致する。ここでは中点を一つだけとったが、これを4点とるのが古典的なRunge-Kutta法であり、比較的実装が容易で4次精度と精度が高いために広く使われている。

その他、これまでの軌跡を覚えておいて、そこから線形補完して場所を予測し、予測点を使って改めて将来の点を修正する、予測子-修正子法も広く使われている。

## 4.2 Integration of Equations of Motion

### 4.2.1 Euler method

さて、分子動力学法の運動方程式を考えよう。簡単のために一次元系を考える。ハミルトニアン$H(p,q)$に支配されている系の運動方程式は以下のように書けるのだった。

$$
\begin{aligned}
\dot{p} &= - \frac{\partial H}{\partial q} \\
\dot{q} &= \frac{\partial H}{\partial p} \\
\end{aligned}
$$

この系は$(p,q)$の変数の組で記述されており、ハミルトニアン$H$は$p,q$の関数であるから、変数の時間微分が自身の関数として表現されている。つまり、ハミルトンの運動方程式も力学系である。力学系であるから、通常の常微分方程式の数値積分法を使うことができる。

いま、系として調和振動子を考えよう。ハミルトニアンは$H=p^2/2 + q^2/2$であり、運動方程式は以下のように書ける。

$$
\begin{aligned}
\dot{p} &= - q \\
\dot{q} &= p \\
\end{aligned}
$$

さて、運動方程式における初期値問題とは、ある時刻$t=0$における初期値$(p(0),q(0))$から、任意の時刻$t$における値$(p(t),q(t))$を求めることである。まず、オイラー法を適用してみよう。時刻$t$の時の微分係数を使って$t+h$の座標を予測するのであるから、

$$
\begin{aligned}
p(t+h) &= p - q h \\
q(t+h) &= q + p h \\
\end{aligned}
$$

と求めることができる。なお、煩雑なので$p(t)$や$q(t)$は$p$、$q$と記述した。さて、時間非依存なハミルトンの運動方程式は、

$$
\dot{H} = \frac{\partial H}{\partial p} \dot{p} + \frac{\partial H}{\partial q} \dot{q} = 0
$$

であるから、ハミルトニアンが保存量となる。調和振動子系なら、$p^2/2+q^2/2$が保存されなければならない。計算してみよう。

$$
p(t+h)^2 + q(t+h)^2 = (p-qh)^2 + (q + ph)^2 = (1+h^2)(p^2+q^2)
$$

つまり、時刻$t$においてエネルギーが$p^2+q^2$であった系は、1ステップ後に$(1+h^2)$倍に、$n$ステップ進むと$(1+h^2)^n$倍になる。「オイラー法が厳密解から指数関数的にずれる」という意味がわかるであろう。

### 4.2.2 Velocity Verlet method

さて、調和振動子の運動方程式をオイラー法で数値積分するとエネルギーがあっという間に発散することがわかった。そこで、Runge-Kuttaや予測子-修正子法のような高次の数値積分法を使っても良いのだが、より簡便で、かつ非常に安定な数値積分法が発見された。velocity Verlet (VV)法である。以下の運動方程式を考えよう。

$$
\begin{aligned}
\dot{p} &= f \\
\dot{q} &= p \\
\end{aligned}
$$

簡単のために一次元系で、質量を1としている。$f(q)$は力であり、ポテンシャル$V(q)$によるものなら$f(q) = -V'(q)$である。この運動方程式に対して、VV法は、以下のように構成される。

まず、位置については、二次までテイラー展開する。

$$
\begin{aligned}
q(t+h) &= q(t) + \dot{q}(t) h + \frac{\ddot{q}(t)h^2}{2} \\
&= q + p h + \frac{f h^2}{2}
\end{aligned}
$$

ここで、運動方程式では位置の時間微分は速度、速度(運動量)の時間微分は力であるから、それぞれ既知なのがミソである。

速度も二次まで展開したいが、速度の微分は力であり、力の微分まで計算するのは面倒だ。そこで、差分を工夫する。先程、すでに時刻$t+h$における位置$q(t+h)$がわかっているので、その場所における力$f(t+h)$を使うことができる。すると、

$$
p(t+h) = p + \frac{f(t+h)+ f(t)}{2}h
$$

と表すことができる。これが二次まで正しいテイラー展開になっていることは容易に確認できるであろう。

さて、このようにして構築されたVV法は、位置に関しても運動量に関しても二次まで正しい展開になっているため、二次精度の数値積分法になっていることが予想される。事実、VV法は二次精度なのだが、この時間発展を行うと **エネルギーが厳密に保存する**。

実際に調和振動子系で確認してみよう。VV法を適用すると

$$
\begin{aligned}
q(t+h) &= q + p h - \frac{q h^2}{2} \\
&= h p + \left( 1 - \frac{h^2}{2}\right)q \\
p(t+h) &= p - \frac{q(t+h) + q(t)}{2} h\\
&= \left(1 - \frac{h^2}{2}\right)p + \left(-h + \frac{h^3}{4} \right)q \end{aligned}
$$

エネルギーの保存を確認すると、

$$
p(t+h)^2 + \left(1-\frac{h^2}{4}\right)q(t+h)^2 = p^2 + \left(1-\frac{h^2}{4}\right)q^2
$$

となっている。つまり、$H = (p^2+q^2)/2$の代わりに、

$$
\tilde{H} = \frac{p^2}{2} + \left(1-\frac{h^2}{4}\right)\frac{q^2}{2}
$$

が厳密に保存される。$h\rightarrow 0$で、$\tilde{H} \rightarrow H$になることに注意しよう。VV法はステップが進んでも、もとのエネルギーからややずれた量が厳密に保存される。実はVV法は、シンプレクティック積分と呼ばれる手法の一種になっている。シンプレクティック積分は、軌道は厳密解からずれるものの、ずれたエネルギーが厳密に保存するため、長時間計算してもエネルギーが一方的にずれたりせず、安定した計算を可能となる。この、シンプレクティック積分が厳密に保存する「少しずれたハミルトニアン」を、「影のハミルトニアン」と呼んだりする。「影のハミルトニアン」の厳密な形が知られているのは線形の時のみである。これについては後述する。

VV法でエネルギーが保存する理由をもう少し詳しく見てみよう。シミュレーションによる時間発展とは、次の時刻$t+h$における$p(t+h), q(t+h)$を、時刻$t$の座標$p(t),q(t)$で表現することである。簡単のため、以下では$P = p(t+h), Q = q(t+h)$と表記することにすると、時間発展は$(p,q)$から$(P,Q)$への座標変換とみなすことができる。一般のこの変換は非線形だが、調和振動子の場合にはこの変形が線形となり、行列で表現することができる。

調和振動子にオイラー法を適用した場合を見てみよう。

$$
\begin{aligned}
P &= p - q h \\
Q &= q + p h \\
\end{aligned}
$$

これが変数変換であることが見やすいように行列表示してみよう。

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

ただし、$U(h)$は以下のような行列である。

$$
U(h)
= 
\begin{pmatrix}
1 &- h \\
-h & 1
\end{pmatrix}
$$

$U(h)$は、$(p,q)$に作用して時間を$h$だけ進める行列(演算子)であるから、時間発展行列(演算子)になっている。さて、これは式だけ見れば一次変換であるから、その変換の性質は行列$U(h)$で決まる。特に重要なのが行列式$|U(h)|$である。一次変換において変換行列の行列式は、変換の前後で面積要素の変化率を表すのであった。オイラー法の場合、

$$
|U(h)| = 1 + h^2 > 1
$$

であるから、この変換が、面積要素を拡大することがわかる。

同様に、VV法の時間発展行列を行列表現すると、

$$
U(h)
= 
\begin{pmatrix}
1 - h^2/2 & h \\
-h + h^3/4 & 1 - h^2/2
\end{pmatrix}
$$

となり、明らかに$|U(h)| = 1$となる。つまり、この行列による変換は、画像を歪めても、面積要素は保存する。実は、この時間発展行列の行列式が1であるということが、エネルギー、すなわちハミルトニアンが厳密に保存することと対応している。調和振動子におけるエネルギー(の2倍)$p^2+q^2$は、円の面積に対応する。オイラー法に対応する変換では、面積要素が毎ステップ拡大してしまうため、エネルギーが増えていってしまうが、VV法に対応する変換では、空間は歪むものの面積要素は保存されるため、エネルギーも保存する、という仕組みになっている。

## 4.3 Symplectic Integrator

### 4.3.1 Matrix Form

VV法が調和振動子の場合にエネルギーを保存するのは、時間発展を記述する行列の行列式が1であることに対応していた。実は、VV法は「Symplectic Integrator」と呼ばれる手法の一種になっており、Symplectic Integratorは指数分解公式から作られる。

ここでは、まずは調和振動子において指数分解公式から時間発展演算子を構築する様子を見てみよう。

運動方程式は行列の形で

$$
\begin{aligned}
\frac{d}{dt}
\begin{pmatrix}
p \\
q
\end{pmatrix}
&=
\underbrace{
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix}
}_{L}
\begin{pmatrix}
p \\
q
\end{pmatrix} \\
&=
L
\begin{pmatrix}
p \\
q
\end{pmatrix}
\end{aligned}
$$

と書ける。ここで$L$は時間微分を表す行列だ。これを形式的に積分すると、

$$
\begin{pmatrix}
P \\
Q
\end{pmatrix}
= 
\underbrace{\mathrm{e}^{hL}}_{U(h)}
\begin{pmatrix}
p \\
q
\end{pmatrix}
$$

つまり、時間微分行列を指数の肩の上に乗せると時間発展行列が得られる。

$$
U(h) = \exp(hL)
$$

線形代数の講義でやったように、$L$を対角化するなどすれば厳密に計算できて、

$$
U(h) =
\begin{pmatrix}
\cos h & \sin h \\
-\sin h & \cos h
\end{pmatrix}
$$

となる。つまり、調和振動子の時間発展は回転で表すことができる。しかし、一般に時間発展行列(演算子)は厳密に求めることができないので、なんらかの近似をすることになる。

まず、オイラー法は以下のように近似している。

$$
U_E(h)=
\begin{pmatrix}
1 & h \\
- h & 1
\end{pmatrix}
$$

厳密解と見比べてみて、1次まで正しい近似になっていることがわかるであろう。

さて、VV法はこうなっていた。

$$
U(h)
= 
\begin{pmatrix}
1 - h^2/2 & h \\
-h + h^3/4 & 1 - h^2/2
\end{pmatrix}
$$

VV法による近似は、テイラー展開の二次まで正しい。ただし、三次の項をうまく付け加えることで、この行列の行列式を1にしているのがポイントである。

さて、今回はある行列$L$を指数の肩に乗せた$\exp(hL)$の表式が厳密に求められたが、そのためには$L$を対角化する必要があり、もし対角化できたら問題は解けたと同義である。そこで、$\exp(ihL)$を近似することを考えよう。

いま、行列$L$が、iL = A + B$と二つの行列の和で書けるとしよう。

一般に$A$と$B$は非可換であるので、

$$
\mathrm{e}^{A+B} \neq \mathrm{e}^{A} \mathrm{e}^{B}
$$

である。しかし、以下の等式が成り立つことが知られている。

$$
\mathrm{e}^{A+B} = \lim_{n \rightarrow}
\left(
\mathrm{e}^{A/n} \mathrm{e}^{B/n}
\right)^n
$$

これをLie-Trotter公式と言う。この式は$n$を無限に飛ばすと厳密だが、それを有限で止めることで、以下のように近似できる。

$$
\begin{aligned}
\exp(h L) &= \exp(h A) \exp(h B) + O(h^2) \\
\exp(h L) &= \exp(h/2 B) \exp(h A) \exp(h/2 B)+ O(h^3) \\
\end{aligned}
$$

これを指数分解公式と呼ぶ。最初の分解が一次、次が二次の公式である。

さて、ここで

$$
A^2 = 0, B^2=0
$$

という性質があったとしよう。すると、これらを指数の肩に乗せても二次以降が消えてしまうので、

$$
\begin{aligned}
\exp(h A) &= I+ hA \\
\exp(h B) &= I+ hB \\
\end{aligned}
$$

と、厳密に値を求めることができる。ただし、$I$は単位行列である。このような分解を利用して数値積分を構成するのがシンプレクティック積分である。

調和振動子の時間微分行列(リュービル演算子)は

$$
L =
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix}
$$

であった。これを

$$
\begin{aligned}
L &=
\begin{pmatrix}
0 & -1 \\
1 & 0
\end{pmatrix} \\
&=
\underbrace{
\begin{pmatrix}
0 & -1 \\
0 & 0
\end{pmatrix}}_A
+
\underbrace{
\begin{pmatrix}
0 & 0 \\
1 & 0
\end{pmatrix}}_B \\
& =
A+B
\end{aligned}
$$

と分解しよう。明らかに$A^2 = B^2 = 0$であるから、

$$
\begin{aligned}
\exp(h A) &= I + h A \\
\exp(h B) &= I + h B \\
\end{aligned}
$$

が成り立つ。さて、まずは一次の指数分解公式

$$
\mathrm{e}^{h(A+B)} \sim \mathrm{e}^{hA}\mathrm{e}^{hB}
$$

を考えてみよう。これを時刻$t$の座標$(p,q)$にかけると時刻$t+h$の座標$(P,Q)$が得られる。つまり、

$$
\begin{pmatrix}
P \\
Q
\end{pmatrix}
=
\mathrm{e}^{hA}\mathrm{e}^{hB}
\begin{pmatrix}
p \\
q
\end{pmatrix}
$$

である。まずは$\mathrm{e}^{h(B)}$を考えよう。

$$
\begin{aligned}
\exp(h B) &= I + h B \\
&= 
\begin{pmatrix}
1 & 0 \\
h & 1
\end{pmatrix}
\end{aligned}
$$
であるから、

$$
\begin{aligned}
\mathrm{e}^{hB}
\begin{pmatrix}
p \\
q
\end{pmatrix}
&=
\begin{pmatrix}
1 & 0 \\
h & 1
\end{pmatrix}
\begin{pmatrix}
p \\
q
\end{pmatrix}\\
&= 
\begin{pmatrix}
p \\
q + ph
\end{pmatrix}
\end{aligned}
$$

つまり、現在の速度$p$で位置が等速直線運動をさせたのと同じである。

次に、$\mathrm{e}^{h(A)}$をかけよう。

$$
\begin{aligned}
\mathrm{e}^{hA}
\begin{pmatrix}
p \\
q+ph
\end{pmatrix}
&=
\begin{pmatrix}
1 & -h \\
0 & 1
\end{pmatrix}
\begin{pmatrix}
p \\
q + ph
\end{pmatrix}\\
&= 
\begin{pmatrix}
(1-h^2)p - hq \\
q + ph
\end{pmatrix}
\end{aligned}
$$

行列の形からわかるように、これは運動量しか更新しない。以上の時間発展をまとめると、

$$
\begin{pmatrix}
P \\
Q
\end{pmatrix}
=
\underbrace{
\begin{pmatrix}
1 - h^2 & -h \\
h & 1
\end{pmatrix}
}_{\tilde{U}_1(h)}
\begin{pmatrix}
p \\
q
\end{pmatrix}
$$

となり、近似された時間発展行列$\tilde{U}_1(h)$は

$$
\tilde{U}_1(h) =
\begin{pmatrix}
1 - h^2 & -h \\
h & 1
\end{pmatrix}
$$

となる。$|\tilde{U}_1(h)|= 1$となっているのがわかるであろう。

このように

* 最初に$q$だけオイラー法で更新する
* 次に **更新した位置を使って** $p$をオイラー法で更新する

として時間積分を構築すると、オイラー法を適用したつもりが、一次のシンプレクティック積分になる(個人的に「なんちゃってオイラー法」と呼んでいる)。
正しいオイラー法は

* 最初に$q$だけオイラー法で更新する
* 次に **更新する前の位置を使って** $p$をオイラー法で更新する

と、更新前の$q$を覚えておかなければならない。

さて、指数分解公式がシンプレクティック積分を作る様子を見てみよう。もともと、時間発展行列は、時間微分行列$L$を指数の肩に乗せたものであり、シンプレクティック性とはその行列式が$1$となること、つまり

$$
|\mathrm{e}^{hL}| = 1
$$

を満たすことであった。さて、指数分解公式は時間微分行列 $L$ を$A^2 = B^2=0$を満たすように$L=A+B$と分解し、それを使って$\exp(hA)$と$\exp(hB)$を組み合わせて時間発展行列を作る。ここで、$A^2=0$であるから、

$$
\exp(h A) = I + h A
$$

さて、行列$A$は

$$
A = 
\begin{pmatrix}
0 & -1 \\
0 & 0
\end{pmatrix}
$$

という形であったから、

$$
1 + hA = 
\begin{pmatrix}
1 & -h \\
0 & 1
\end{pmatrix}
$$

明らかに

$$
|\exp(h A)| = |I + h A| = 1
$$

と、行列式1になることがわかるであろう。

念の為、一般的に$A^2=0$なら$|\exp(h A)|=1$であることを証明しておこう。

まず、$X^n=0$と、べき乗してゼロになる行列を冪零行列と言う。冪零行列の固有値は全て0である。なぜなら固有値$\lambda$と固有ベクトル$v$には$Xv = \lambda x$の関係があるが、$X^nv = \lambda^n x$であり、$X^n = 0$であるから$\lambda^n=0$、したがって$\lambda=0$である。さて、ある行列$X$の行列式$|X|$は、固有値の積である。つまり、$X$の固有値を$\lambda_i$とすると、

$$
|X| = \lambda_1 \lambda_2 \cdots \lambda_N
$$

また、行列指数関数$\mathrm{e}^X$の固有値は、固有値を指数の肩に乗せたものだ。したがって

$$
|\mathrm{e}^X| = \mathrm{e}^{\lambda_1}\mathrm{e}^{\lambda_2} \cdots \mathrm{e}^{\lambda_N}
$$

$A^2=0$であるから$A$は冪零行列であり、冪零行列の固有値は全てゼロであるから、

$$
|\exp(h A)| = |\exp(A)|^h = 1^h = 1
$$

以上で$A^2=0$なら$|\exp(h A)|=1$が証明された。$B$も同様である。

行列の積の行列式は、行列式の積になるから、$\mathrm{e}^A$と$\mathrm{e}^B$の積で作られる行列はかならず行列式が1となる。つまり、変換の面積要素が保存される。これが指数分解公式がシンプレクティック積分を作る理由となる。

以上を行列の言葉でまとめておこう。

* ハミルトンの運動方程式に対応する時間微分行列$L$は歪エルミート行列となる(これを嫌って、通常は時間微分演算子を$iL$として、$L$をエルミートに取る)
* 時間発展行列$U=\mathrm{e}^{L}$は、歪エルミート行列$L$を指数関数の肩に乗せたものなので、ユニタリ行列になる
* ユニタリ行列の行列式は1となる。
* 行列式が1となる行列による変換は、面積要素を保存する。これによりエネルギーが保存する。
* 指数分解公式は、時間微分行列$L$を冪零行列の和$A+B$で表し、時間発展演算子を冪零行列を指数の肩に乗せたもので表現する方法である。
* 冪零行列を指数の肩に乗せると行列式が1となるので、$|\exp(h A)|=|\exp(h B)|=1$である。以上からこの積で作られた行列の行列式が1となり、時間発展がシンプレクティックとなる

### 4.3.2 Liouville Operator

調和振動子の場合は、時間微分演算子、時間発展演算子が行列で書けた。しかし、一般に時刻$t$の座標$(p,q)$から時刻$t+h$の座標$(P,Q)$への写像は非線形となり、それぞれの演算子が行列では書けない。この場合の指数分解公式と、シンプレクティック性について見てみよう。

一般の時間微分演算子を考える。

$$
\begin{pmatrix}
\dot{p} \\
\dot{q}
\end{pmatrix}
= iL 
\begin{pmatrix}
p\\
q
\end{pmatrix}
$$

ここに現れる演算子$iL$を、Liouville Operatorと呼ぶ。ハミルトンの運動方程式におけるリュービル演算子は以下のように書ける。

$$
i \mathcal{L} = \underbrace{\frac{\partial H}{\partial p} \frac{\partial}{\partial q}}_{i\mathcal{L}_K} + \underbrace{ \left(-
\frac{\partial H}{\partial q} \frac{\partial}{\partial p}\right)}_{i\mathcal{L}_V}
$$

これを見て

$$
i\mathcal{L} = i \mathcal{L}_K + i \mathcal{L}_V
$$

という分解を自然に思いつくであろう。さて、いまハミルトニアンが以下の様に、運動量のみに依存する項$K$と座標のみに依存する項$V$の和で書けていたとする。

$$
H(p,q) = K(p) + V(q)
$$

ただし、$K(p) = p^2/2m$である。これを自然ハミルトニアンと呼ぶ。この場合、先程分解した二つの演算子が以下のようになる。

$$
\begin{aligned}
i\mathcal{L}_K = \frac{\partial H}{\partial p} \frac{\partial}{\partial q} = \frac{\partial K}{\partial p} \frac{\partial}{\partial q} \\
i\mathcal{L}_V = - \frac{\partial H}{\partial q} \frac{\partial}{\partial p} = - \frac{\partial V}{\partial p} \frac{\partial}{\partial p} \\
\end{aligned}
$$

ここで、$\partial_q$の係数が$p$のみに依存し、$\partial_p$の係数が$q$のみに依存することに注意。さて、この演算子を$p$や$q$に演算してみよう。

まず$i\mathrm{L}_K$を$p$にかけると、$\partial_q$で消えるのでゼロである。$q$にかけると、

$$
i\mathcal{L}_K q = \frac{\partial K}{\partial p} \frac{\partial q}{\partial q} = \frac{\partial K}{\partial p} 
$$

$K$は$p$のみの関数であるから、さらに$q$で偏微分するとゼロになる。従って

$$
(i\mathcal{L}_K)^2 q = 0 
$$

全く同様に、

$$
(i\mathcal{L}_V)^2 p = 0 
$$

ここから、この演算子を指数関数の肩に乗せたものを$p$や$q$に演算した結果を厳密に計算することができる。

$$
\begin{aligned}
\mathrm{e}^{i h \mathcal{L}_K} p &= 0 \\
\mathrm{e}^{i h \mathcal{L}_K} q &= q + h \underbrace{\frac{\partial K}{\partial p}}_v \\
\mathrm{e}^{i h \mathcal{L}_V} q &= 0 \\
\mathrm{e}^{i h \mathcal{L}_V} q &= p - h \underbrace{\frac{\partial V}{\partial q}}_f \\
\end{aligned}
$$

ここで、$\partial_p K = p/m = v$は速度、$\partial_q V = V'(q) = f$は力であるから、それぞれ「時間$h$の間等速直線運動をした時の座標の変化」「時間$h$の間、力$f$を受け続けた運動量の変化」を表している。あとは全く同様に指数分解公式を用いることで、数値積分法を構築できる。

一次のシンプレクティック積分法であれば、

$$
\mathrm{e}^{i h \mathcal{L}} \sim  \mathrm{e}^{i h \mathcal{L}_V}\mathrm{e}^{i h \mathcal{L}_K}
$$

と分解できるので、

$$
\begin{pmatrix}
P \\
Q
\end{pmatrix}
=
\mathrm{e}^{i h \mathcal{L}_V}\mathrm{e}^{i h \mathcal{L}_K}
\begin{pmatrix}
p \\
q
\end{pmatrix}
$$

を計算すれば良い。これは、

* 最初に現在の速度で時間$h$だけ等速直線運動をさせて
* 次に、更新された座標を使って計算される力が$h$だけ持続した時の力積により運動量を変化させる

というアルゴリズムになっている。

次に、二次の分解公式を考えよう。これは

$$
\mathrm{e}^{i h \mathcal{L}} \sim  \mathrm{e}^{i h \mathcal{L}_V/2}\mathrm{e}^{i h \mathcal{L}_K} \mathrm{e}^{i h \mathcal{L}_V/2}
$$

と分解する方法だ。対応する数値積分アルゴリズムは

1. 最初に現在の力が$h$だけ持続した時の力積により運動量を変化させ
2. 次に更新された運動量で時間$h$だけ等速直線運動をさせ、
3. 最後に更新された座標における力が$h/2$だけ持続した場合の力積変化を計算する

というアルゴリズムになっている。指数分解公式を使うというと難しく感じるが、要するに座標と運動量を更新する際、どちらかが止まっている(定数である)と思って、片方を更新するのを繰り返しているだけである。

### 4.3.3. Symplecity of Velocity Verlet Algorithm

では、最後にVV法が二次のシンプレクティック積分になっていることを示そう。自然ハミルトニアン

$$
H = p^2 + V(q)
$$

を考える。簡単のため、質量を$1$としている。

二次のシンプレクティック積分のアルゴリズムは

1. 最初に現在の力が$h$だけ持続した時の力積により運動量を変化させ
2. 次に更新された運動量で時間$h$だけ等速直線運動をさせ、
3. 最後に更新された座標における力が$h/2$だけ持続した場合の力積変化を計算する

となっていた。時刻$t$において、座標が$(p, q)$であったとしよう。このアルゴリズムにより更新された時刻$t+h$における座標を$(P,Q)$とする。

まず、現在かかっている力が$h/2$だけ持続したとして運動量を変化させる

$$
p(t+h/2) = p + \frac{f(t)h}{2}
$$

次に、更新された運動量で時間$h$だけ等速直線運動をさせる。

$$
Q = q + p(t+h/2) h
$$

最後に、更新された座標$Q$における力が$h/2$だけ持続した場合の運動量変化を考える。

$$
P = p(t+h/2) + \frac{f(t+h)h}{2}
$$

以上で、$(p,q)$から$(P,Q)$への写像、すなわち時間積分が完成した。$p(t+h/2)$を消去すると、

$$
\begin{aligned}
Q &= q + p h + \frac{f h^2}{2} \\
P &= p + \frac{f(t) + f(t+h)}{2} h
\end{aligned}
$$

これは、VV法に他ならない。

先程、二次の指数分解公式として

$$
\mathrm{e}^{i h \mathcal{L}} \sim  \mathrm{e}^{i h \mathcal{L}_V/2}\mathrm{e}^{i h \mathcal{L}_K} \mathrm{e}^{i h \mathcal{L}_V/2}
$$

を考えた。この分解を逆にして、

$$
\mathrm{e}^{i h \mathcal{L}} \sim  \mathrm{e}^{i h \mathcal{L}_K/2}\mathrm{e}^{i h \mathcal{L}_V} \mathrm{e}^{i h \mathcal{L}_K/2}
$$

とすると、数値積分法として

1. まず$h/2$だけ等速直線運動をさせる
2. 現在の力が$h$だけ持続したとして運動量を変化させる
3. 最後に更新された運動量で$h/2$だけ等速直線運動をさせる

というステップを繰り返すアルゴリズムが構築できる。このステップを繰り返すと、同じ速度で座標を$h/2$の時間だけ二回更新するのが無駄である。そこで、

1. $h$だけ等速直線運動をさせる
2. 現在の力が$h$だけ持続したとして運動量を変化させる

というステップを繰り返しつつ、もし座標の情報が欲しい場合は時刻を$h/2$だけずらす、という方法が考えられた。これは座標と運動量が時間$h/2$だけずれて交互に更新されるように見えることからLeap-frog法と呼ばれる。

やっている計算は一次のシンプレクティック積分と変わらないのだが、観測のタイミングが異なると二次になるのが面白い点である。

### 4.3.4 Shadow Hamiltonian

シンプレクティック積分は、元のハミルトニアンから少しだけずれた「影のハミルトニアン」を厳密に保存する。一般に影のハミルトニアンが存在するか、存在するとして時間刻みに対する収束半径はどれくらいかは知られていないが、系が線形の場合は影のハミルトニアンを厳密に求めることができる。

以下、一次元調和振動子系において影のハミルトニアンを求めてみよう。

運動方程式は以下の通り。

$$
\begin{aligned}
\dot{p} &= -q\\
\dot{q} &= p
\end{aligned}
$$

リュービル演算子に対応する時間微分行列は

$$
L =
\begin{pmatrix}
0 & -1\\
1 & 0
\end{pmatrix}
$$

となる。これを

$$
\begin{aligned}
L &= L_A + A_B \\
&=
\begin{pmatrix}
0 & -1\\
0 & 0
\end{pmatrix}
+
\begin{pmatrix}
0 & 0\\
1 & 0
\end{pmatrix}
\end{aligned}
$$

と分離しよう。時間発展行列$U$は

$$
U(h) = \mathrm{e}^{h(L_A+L_B)}
$$

であるが、これを

$$
\begin{aligned}
A &\equiv \mathrm{e}^{hL_A} =  
\begin{pmatrix}
0 & -h\\
0 & 0
\end{pmatrix} \\
B &\equiv \mathrm{e}^{hL_B} =  
\begin{pmatrix}
0 & 0\\
h & 0
\end{pmatrix} 
\end{aligned}
$$

で近似することを考えよう。

もともと、$p^2+q^2$が保存量であったことから、影のハミルトニアンが、二次形式

$$
\tilde{H} = 
\begin{pmatrix}
p & q
\end{pmatrix}
X
\begin{pmatrix}
p\\
q
\end{pmatrix}
$$

という形を仮定しよう(煩雑になるので1/2のファクターは除いてある)。厳密な保存量は$X=I$である。

さて、時間発展演算子として、一次の近似、

$$
\tilde{U}_1(h) = BA
$$

を考えよう。これにより$(p,q)$が$(P,Q)$になったとすると、

$$
\begin{pmatrix}
P\\
Q
\end{pmatrix}
=
\tilde{U}_1(h)
\begin{pmatrix}
p\\
q
\end{pmatrix}
$$


影のハミルトニアンが保存することを要請すると、

$$
\begin{pmatrix}
P & Q
\end{pmatrix}
X
\begin{pmatrix}
P\\
Q
\end{pmatrix}
=
\begin{pmatrix}
p & q
\end{pmatrix}
X
\begin{pmatrix}
p\\
q
\end{pmatrix}
$$

ここから、

$$
\tilde{U}_1^t X \tilde{U}_1 = X
$$

であるから、

$$
A^t B^t X BA = X
$$

となる。ここで、$BA^t = AB^t = I$であることから、両辺に左から$AB$をかけると、

$$
XBA = ABX
$$

この$X$は一意には決まらないが、例えば$X=A$と取れば良いことがわかる。ここから、

$$
\begin{aligned}
\tilde{H}_1 &= 
\begin{pmatrix}
p & q
\end{pmatrix}
X
\begin{pmatrix}
p\\
q
\end{pmatrix} \\
&=
p^2 - h pq + q^2
\end{aligned}
$$

これが一次の指数分解公式で作られたシンプレクティック積分の「影のハミルトニアン」である。もとのハミルトニアンから$h$の一次でずれていることがわかる。

二次の場合も同様に計算できる。計算の便利のために、以下の行列を定義する。

$$
\begin{aligned}
A_\mathrm{h} &= 
\begin{pmatrix}
0 & -h/2\\
0 & 0
\end{pmatrix} \\
B_\mathrm{h} &= 
\begin{pmatrix}
0 & 0\\
h/2 & 0
\end{pmatrix} \\
\end{aligned}
$$

ちょっとややこしいが、$A_\mathrm{h}$の添え字はHalfを表す。ここで、

$$
\begin{aligned}
A_\mathrm{h}^2 &= A\\
B_\mathrm{h}^2 &= B\\
\end{aligned}
$$

であることに注意。

VV法の時間発展演算子は

$$
\tilde{U}_2(h) = A_\mathrm{h} B A_\mathrm{h}
$$

であるから、

$$
A_\mathrm{h}^t B^t A_\mathrm{h}^t X A_\mathrm{h} B A_\mathrm{h} = X
$$

を満たす行列$X$を見つければ良い。

ここで、

$$
\begin{aligned}
B_\mathrm{h} A_\mathrm{h}^t &= I\\
A_\mathrm{h} B_\mathrm{h}^t &= I\\
A_\mathrm{h}^2 &= A\\
B_\mathrm{h}^2 &= B
\end{aligned}
$$

であることを使って、両辺に左から$B_\mathrm{h} A_\mathrm{h}^2 B_\mathrm{h}$をかけると、

$$
X A_\mathrm{h} B_\mathrm{h} B_\mathrm{h} A_\mathrm{h}
= B_\mathrm{h} A_\mathrm{h} A_\mathrm{h} B_\mathrm{h} X
$$

両辺が等しくなるような$X$は、例えば

$$
X = B_\mathrm{h} A_\mathrm{h}
$$

とすれば良い。この時、影のハミルトニアンは

$$
\begin{aligned}
\tilde{H}_2 &= 
\begin{pmatrix}
p & q
\end{pmatrix}
B_\mathrm{h} A_\mathrm{h}
\begin{pmatrix}
p\\
q
\end{pmatrix} \\
&=
p^2 + \left(1 - \frac{h^2}{4} \right)q^2
\end{aligned}
$$

と求まる。

影のハミルトニアンの形

$$
\begin{aligned}
\tilde{H}_1 &= p^2 - hpq + q^2\\
\tilde{H}_2 &= p^2 + \left(1-\frac{h^2}{4}\right)q^2
\end{aligned}
$$

の形を見ると、$h$が$2$より大きくなると、形が楕円型から双曲型に変化し、数値計算が破綻することが予想される。

実際、調和振動子の場合には収束半径まで含めて影のハミルトニアンが厳密に計算できることが知られている(Kobayasih 2007)。

ここでは、影のハミルトニアンが$p,q$の二次形式であることを仮定して発見法的に求めたが、一般の場合において影のハミルトニアンが厳密に求められた例や、時間刻みに対する収束半径が求められた例を筆者は知らない。

[1] H. Kobayashi, Phys. Lett. A, Vol. 371, Issues 5–6, 26 November 2007, Pages 360-362, [doi:10.1016/j.physleta.2007.06.037](https://doi.org/10.1016/j.physleta.2007.06.037}
