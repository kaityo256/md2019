# 分子動力学法の理論と実装

<a href="https://github.com/kaityo256/md2019"> <div class="btn-square"><i class="fab fa-github"></i> View on GitHub</div></a>

## この文書について

これは金沢大学で行われる集中講義の講義ノートにする予定。

## 内容

* 本講義の概要と目的
* 分子動力学法の理論的背景
  * 圧力と界面張力
  * 温度制御とエルゴード性
  * 分子動力学法における数値積分法
* 分子動力学法の実装
  * 分子動力学法の実装と基本的アルゴリズム
  * 分子動力学法の高速化手法について
  * 分子動力学法の並列化とプログラム設計

## 分子動力学法とは

分子動力学法(Moleclar Dynamics method, MD)とは、ニュートンもしくはハミルトンの運動方程式を数値的に解くことで粒子系を時間発展させる手法である。原子や分子、あるいはそれを粗視化した「つぶつぶ」の間にかかる力を計算し、その力によって速度を更新し、その速度によって位置を更新する、というプロセスを繰り返す。MDの基礎方程式が運動方程式という身近なものであることから、なんとなく「とっつきやすい」「わかりやすい」手法のように見える。しかし、MDという手法を突き詰めて考えてみるとだんだんよくわからなくなってくる。ここで温度と言っているのはどういう量なのか？圧力とは何か？そもそもMDにおける時間発展とはなんなのか？

### 変数と観測量

数値計算は、支配方程式を数値的に解く方法論だ。ここで重要となるのは、変数と観測量の区別である。変数(Variable)とは「我々がa prioriに認める量」であり、観測量(Observable)は、「変数から導かれる量」だ。何が変数で何が観測量かは支配方程式によって異なる。

例えば、一次元熱伝導方程式を考えよう。

$$
\frac{\partial T}{\partial t} = \kappa \frac{\partial^2 T}{\partial x^2}
$$

この方程式では、温度$T$は変数である。つまり、「この世界には温度という量がある」ということを無条件に認めている。

また、ナビエ・ストークス方程式を考えてみよう。粘性率一定、非圧縮を仮定し、外場を$\vec{F}$とすると以下のようになるだろう。

$$
\frac{\partial \vec{v}}{\partial t}
+ \left(\vec{v}\cdot \nabla \right)
= -\frac{1}{\rho} \nabla p + \nu \Delta \vec{v} + \vec{F}
$$

この方程式では、速度場$v$に加えて、圧力場$p$も変数である。

さて、MDではどうだろうか？支配方程式としてハミルトンの運動方程式を採用するとこんな方程式となる。

$$
\begin{aligned}
\dot{p_i} &= -\frac{\partial H}{\partial q_i} \\
\dot{q_i} &= \frac{\partial H}{\partial p_i} 
\end{aligned}
$$

変数は座標$q_i$と運動量$p_i$であり、ここには温度や圧力は表れない。したがって、温度や圧力は$q_i$や$p_i$から定義される観測量となる。熱伝導方程式においては、温度を変数にとっているから「温度とは何か」という問いは意味をもたない。しかし、分子動力学法においては温度は定義するものであるので、「私はこの量を温度と呼ぶ」と宣言することになる。そこには様々な温度の定義があり得る。もちろん平衡状態では熱力学と整合するように選ばれるべきだが、後で見るように、異なる温度の定義は非平衡状態では一致しない。MDにおいては「温度とは何か」という問いは本質的な問いとなる。同様にMDにおいては圧力も定義するものであり、その定義には複数の立場があり得る。特に界面がある系における圧力定義は微妙な問題となるので、後でそれにも触れる。

## 温度とは

### マクスウェル分布と温度

温度とは何であろうか？我々は日常的に「熱い/冷たい」「暑い/寒い」と、温度を感じており、温度という量が自然に存在することに疑問を抱かない。しかし、以前に述べたように分子動力学法における変数は運動量$p_i$と座標$q_i$のみであるから、この変数の組$\{p_i, q_i\}$で温度を定義しなければならない。

広く受け入れられているのは「温度は運動エネルギーに比例する」という事実を利用した定義だ。三次元$N$粒子系を考えると、分子動力学法における温度$T$は運動エネルギーを$K$として

$$
T = \frac{2K}{3N k_\mathrm{B}}
$$

という量の時間平均で定義される。ただし、$k_\mathrm{B}$はボルツマン定数である。簡単のため、全ての粒子の質量$m$が等しいとすると、運動エネルギー$K$の定義は

$$
K = \sum_i \frac{p_i^2}{2m}
$$

で与えられる。この量が温度と結びつくのは、運動量がマクスウェル分布に従うことから説明される。すなわち、温度$T$の粒子系における運動量分布は、

$$
f(\{p_i \}) \propto \exp\left( \frac{\sum_i p_i^2}{2m k_\mathrm{B} T} \right)
$$

に従うことが知られており、ここからいわゆるエネルギー当分配則

$$
\left< \frac{p_i^2}{2m}\right> = \frac{1}{2} k_\mathrm{B} T
$$

が導かれる。これを認めれば運動エネルギーと温度が結びつく、というストーリーである。

### カノニカル分布と温度

さて、マクスウェル分布は観測事実であるが、これを別の基本原理から「導く」ことを考えてみよう。

まず、ハミルトニアンを議論の出発点としよう。三次元$N$粒子系で、全ての粒子の質量が等しい場合は、以下のようなハミルトニアンになるだろう。

$$
H = \sum_i \frac{p_i^2}{2m} + V(\{q_i\})
$$

ここで、$V(\{q_i\})$は座標にのみ依存するポテンシャル関数とする。さて、$p_i$と、$q_i$が張る位相空間を$\Gamma$としよう。この$\Gamma$上に、分布関数$f(\{p_i, q_i\})$を定める。これは、位相空間上の点$\{p_i, q_i\}$における系の存在確率をあらわす。したがって、分布関数$f$は規格化されていなければならない。

$$
\int f d \Gamma = 1
$$

次に、天下りだが、この分布関数の対数をとったものの平均に$- k_\mathrm{B}$をかけた量を考え、$S$と名前をつけよう。

$$
S = -k_\mathrm{B} \int f \ln f  d \Gamma
$$

もちろんこの量はエントロピーであるが、それがどういう量かはよくわからない。さて、なぜだかわからないが、あなたははこの$S$という量を最大化したくなった。いま、系のエネルギーが$E$であるとしよう。この条件で$S$を最大化するには、規格化条件をラグランジュの未定乗数法で扱って、

$$
I = \alpha \int f d \Gamma + \int f \ln f d \Gamma
$$

の変分をとればよい。計算の便利のために、ラグランジュの未定乗数$\alpha$にボルツマン定数を含めている。変分条件$\delta I = 0$よりただちに、

$$
f = \frac{1}{\Omega}
$$

が導かれる。ただし、$\Omega$は規格化条件から決まる以下の量である。

$$
\Omega = \int \delta(H-E) d \Gamma
$$

これは等エネルギー面の面積に他ならない。エネルギー一定条件のもとで$S$分布関数が一定になることが導かれた

ミクロカノニカル分布であり、ミクロカノニカルであれば分布関数は一定である、ということに対応している。

ミクロカノニカルでは、エネルギー一定条件を考えたが、これを少し緩めて、エネルギーの期待値一定の条件にしてみよう。

分布関数$f$により、この系における物理量の期待値$\left<A\right>$を以下のように定義する。

$$
\left< A \right> \equiv \int A f d \Gamma
$$

ハミルトニアン$H$の期待値が内部エネルギー$U$である。

$$
U \equiv \left< H \right> = \int Hf d \Gamma
$$

分布関数$f$が規格化条件を満たし、内部エネルギー$U$が一定である、という条件で$S$を最大化することを考える。すると、ラグランジュの未定定数$\alpha, \beta$を用いて、

$$
I = \alpha \int f d \Gamma
+ \beta \int H f d \Gamma
+ \int f \ln f d \Gamma
$$

の変分を取ればよい。先ほどと同様、後の便利のために$\alpha, \beta$にボルツマン定数を吸収させ、かつ$S$に関する条件式に負符号をつけていることに注意(これが自由エネルギー$F = U - TS$になるように、$S$の定義の負符号と変分の負符号がキャンセルするように選んでいる)。変分条件$\delta I = 0$より、

$$
\alpha + \beta H + 1 + \ln f = 0
$$

ここからただちに分布関数が

$$
f = \mathrm{e}^{-(\alpha + 1)} \mathrm{e}^{-\beta H} d \Gamma
$$

となる。規格化条件から

$$
\mathrm{e}^{\alpha + 1} = \int  \mathrm{e}^{-\beta H} d \Gamma \equiv Z
$$

と$\alpha$が求まるため、最終的に分布関数は

$$
f = Z^{-1} \mathrm{e}^{-\beta H}
$$

となる。これはボルツマン重みに従う分布、すなわちカノニカル分布に他ならない。

さて、ここで$S$をエントロピーだと思うことにしよう。$S$と$U$の間に熱力学関係式

$$
dS = \frac{dU}{T}
$$

が成り立つことを要請すると、残るラグランジュの未定定数$\beta$は

$$
\beta = \frac{1}{k_\mathrm{B} T}
$$

と、温度と結びつく。このストーリーでは、「まず世の中にはエネルギー$U$とエントロピー$S$というものが存在し、$U$一定の条件で$S$を最大化しようとするとカノニカル分布が導かれ、エネルギーとエントロピーの熱力学関係式を要請すると、ラグランジュの未定定数$\beta$が逆温度であることがわかる」という流れになっている。

### 一般化ビリアル

さて、先程のような温度の導入をすると、運動エネルギー以外にも温度の定義ができることがわかる。

カノニカル分布を持つ系において、ビリアルと呼ばれる量$p_i \partial_{p_i} H$を考えてみよう。このビリアルについて、以下の関係が成り立つ。

$$
\left<p_i \frac{\partial H}{\partial p_i} \right> = \frac{1}{\beta}
$$

これは部分積分を使うだけで容易に証明できる。

$$
\begin{aligned}
\left<p_i \frac{\partial H}{\partial p_i} \right> &= Z^{-1} \int p_i \frac{\partial H}{\partial p_i} \exp(-\beta H) d \Gamma \\
&=  Z^{-1} \int \frac{p_i}{(-\beta)} \frac{\partial}{\partial  p_i} \left(\exp(-\beta H)\right)  d\Gamma\\
&= \frac{Z^{-1}}{\beta} \int \frac{\partial p_i}{\partial p_i}  \exp(-\beta H)  d \Gamma\\
&= \frac{1}{\beta}
\end{aligned}
$$

ここで、ハミルトニアンの定義から

$$
p_i \frac{\partial H}{\partial p_i} = \frac{p_i}{m}
$$

であるから、等分配則

$$
\left<  \frac{p_i}{2m} \right> = \frac{1}{2 \beta}
$$

が導かれた。

さて、導出を見ると、$p_i$を$q_i$に変えても全く同様なことができる。

$$
\left<q_i \frac{\partial H}{\partial q_i}  \right> = \frac{1}{\beta}
$$

運動量のビリアルから定まる温度を **運動温度 (Kinetic Temperature)** 、**座標から決まる温度を状態温度(Configuration Temperature)** と呼んで区別することがある。この二つの温度は、系がカノニカル分布に従う場合、すなわち平衡状態にある場合は一致するが、一般に非平衡状態においては一致しない。

ここまでの議論において、ビリアルの期待値から温度が出てくるのは、要するに分布関数が$\exp(-\beta H)$という形をしていることと部分積分を使っているだけなので、さらに一般化することができる。位相空間$\{p_i, q_i\}$をまとめて$\{z_i\}$と表記しよう。これは運動量と座標をまとめて連番を振ったもの、つまり、

$$
\{z_i\} = \{p_1, p_2, \cdots, q_1, q_2, \cdots \}
$$

である。何か適当な量$B_i$と、$\partial_{z_i} H$との積$B_i \partial_{z_i} H$を考える。部分積分により、

$$
\left<B_i \frac{\partial H}{\partial z_i} \right>
= \frac{1}{\beta}
\left<\frac{\partial B_i}{\partial z_i} \right>
$$

さて、いま位相空間$\Gamma = \{z_i\}$中に、ベクトル場$\vec{B}$が定義されているとしよう。先ほどの式はそれぞれの座標成分$z_i$ごとに成り立つので、全ての成分について和をとると、以下のようにベクトルの形で書くことができる。

$$
\left<\vec{B} \cdot \nabla H \right>
= \frac{1}{\beta}
\left<\nabla \cdot \vec{B} \right>
$$

この式の左辺を一般化ビリアルと呼ぶ。$\vec{B}$が運動量成分しか含まない場合、すなわち

$$
\vec{B} = (p_1, p_2, \cdots, p_N, 0, \cdots, 0)
$$

の場合に導かれる温度は運動温度となり、いわゆるエネルギー等分配則を表す。逆に、$\vec{B}$が座標成分しか含まない場合は、状態温度を導く。

### 温度の幾何学的定義

先ほど、エントロピー最大化条件からカノニカル分布を導き、そこからビリアルによって温度を定義した。しかし、分子動力学法において、ハミルトンの運動方程式に従う時間発展を考えるとエネルギーが一定に保たれるのであるから、結果として得られる分布は(系がエルゴード的であれば)ミクロカノニカル分布となる。もちろん、系の自由度が大きくなれば分布がデルタ関数的になり、ミクロカノニカルとカノニカルの区別はほとんどなくなるのだが、ハミルトンの運動方程式が厳密にエネルギーを保存する以上、まずミクロカノニカルで温度を定義し、それが$N \rightarrow \infty$の極限でカノニカル分布における温度と一致することを確認したくなる。

以下では、「温度の幾何学的定義」により、ミクロカノニカル分布における温度を定義してみよう。

TODO: 温度の幾何学的定義の説明

## [分子動力学法と圧力](pressure/README.md)


## 温度制御

熱力学における自然な変数には、内部エネルギーや圧力、体積、温度やエントロピーなどがある。これらの変数は示量性の量と示強性の量にわけることができる。さて、「お互いにかけてエネルギーになる示量性の量と示強性の量の組」は共役な量と呼ばれる。

例えば体積$V$と圧力$P$、エントロピー$S$と温度$T$、粒子数$N$と化学ポテンシャル$\mu$が互いに共役であり、それぞれ前者が示量性、後者が示強性の量である。

さて、世の中の量は「a priori」に認める量と、その量から導かれる量の二種類があるのであった。我々は一般に示量性の量をa prioriに認めることが多い。例えば「長さ」は知っているものとするから、体積$V$はa prioriに認め、共役な量である圧力$P$はそこから導かれる量とする。個数$N$も基本的な量として、相方である化学ポテンシャル$\mu$はそこから定義される量である。

分子動力学法では、示強性の量を制御したい場合がある。先ほど、圧力$P$を制御するのに、共役な量の相方である体積$V$をコントロールした。これはわかりやすい。

では、温度$T$を制御するにはどうすればよいだろうか？圧力制御からの類推では、エントロピーをコントロールすることになるが、エントロピーとはどうやってコントロールすればよいのだろうか。そもそもエントロピーと温度、どちらが基本的な量であり、どちらを従属的な量であろうか？

この問いへの回答は筆者の能力を超える。ここでは「両方の立場があり得る」とだけコメントしておく。例えば田崎さんの熱力学の教科書は温度を基本的な量に取ってエントロピーを導く形式であり、清水さんの教科書はエントロピーを基本的な量に取って温度を導く形式である。

TODO: 能勢のハミルトニアン

TODO: 時間の制御？

TODO: Nose-Hoover法

TODO: Nose-Hooverの一般化、Nose-Hoover-Family

TODO: Nose-Hooverの保存量とは

TODO: Nose-Hooverと調和振動子とエルゴード性
