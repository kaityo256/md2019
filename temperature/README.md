# 分子動力学法における温度制御

熱力学における自然な変数には、内部エネルギーや圧力、体積、温度やエントロピーなどがある。これらの変数は示量性の量と示強性の量にわけることができる。さて、「お互いにかけてエネルギーになる示量性の量と示強性の量の組」は共役な量と呼ばれる。

例えば体積$V$と圧力$P$、エントロピー$S$と温度$T$、粒子数$N$と化学ポテンシャル$\mu$が互いに共役であり、それぞれ前者が示量性、後者が示強性の量である。

さて、世の中の量は「a priori」に認める量と、その量から導かれる量の二種類があるのであった。我々は一般に示量性の量をa prioriに認めることが多い。例えば「長さ」は知っているものとするから、体積$V$はa prioriに認め、共役な量である圧力$P$はそこから導かれる量とする。個数$N$も基本的な量として、相方である化学ポテンシャル$\mu$はそこから定義される量である。

分子動力学法では、示強性の量を制御したい場合がある。先ほど、圧力$P$を制御するのに、共役な量の相方である体積$V$をコントロールした。これはわかりやすい。

では、温度$T$を制御するにはどうすればよいだろうか？圧力制御からの類推では、エントロピーをコントロールすることになるが、エントロピーとはどうやってコントロールすればよいのだろうか。そもそもエントロピーと温度、どちらが基本的な量であり、どちらを従属的な量であろうか？

この問いへの回答は筆者の能力を超える。ここでは「両方の立場があり得る」とだけコメントしておく。例えば田崎さんの熱力学の教科書は温度を基本的な量に取ってエントロピーを導く形式であり、清水さんの教科書はエントロピーを基本的な量に取って温度を導く形式である。

## Nose-Hoover法以外の方法

* Velocity Scaling
* Gaussian Thermostat
* Berendsen Thermostat

TODO: 書く

## 能勢の方法

TODO: 書く？さらっと流すだけ？

## 能勢のハミルトニアンからNoseHoover法へ

能勢のハミルトニアンは以下で与えられる。

$$
H_\mathrm{Nose} = H(p/s, q) + \frac{p_s^2}{2Q} + \frac{\ln s}{\beta}
$$

運動方程式は以下の通り。後で時間のスケールをする都合上、時間微分を陽に書いている。

$$
\begin{aligned}
\frac{dq}{dt} &= \frac{1}{s} \partial_1 H(p/s, q) \\
\frac{dp}{dt} &= -\frac{\partial H(p/s,q)}{\partial q} \\
\frac{d p_s}{dt} &= \frac{p}{s^2}\partial_1 H(p/s, q) - \frac{1}{s\beta}\\
\frac{ds}{dt} &= \frac{p_s}{Q}
\end{aligned}
$$

ただし、$\partial_1$は多変数関数の一つ目の変数に関する微分という意味である。

ここで、$p/s = p'$、$dt' = dt/s$という変換を考える。

まず、$q$について、$p/s = p'$を考えると、

$$
\frac{dq}{dt} = \frac{1}{s} \frac{\partial H}{\partial p'}
$$

時間を$t$から$t'$に変えると、

$$
\frac{dq}{dt'} = \frac{\partial H}{\partial p'}
$$

次に、$p'=p/s$について考えてみよう。$p'$を時間微分すると、

$$
\begin{aligned}
\frac{d p'}{dt} &= \frac{1}{s}\frac{d p}{dt} -\frac{p}{s^2} \frac{ds}{dt} \\
&= -\frac{1}{s} \frac{\partial H}{\partial q} - \frac{p'}{s} \frac{p_s}{Q}
\end{aligned}
$$

時間を$t$から$t'$に変えると、

$$
\frac{dp'}{dt'} = - \frac{\partial H}{\partial q} -p' \zeta
$$

ただし$\zeta = p_s/Q$である。

次に、$\zeta$の時間微分を考えよう。

$$
\begin{aligned}
\frac{d \zeta}{dt} &= \frac{1}{Q} \frac{d p_s}{d t} \\
&= \frac{1}{Q} \left(
    \frac{p}{s^2}\partial_1 H(p/s, q) - \frac{1}{s\beta}
\right) \\
&= \frac{1}{Qs}
\left(
p' \frac{\partial H}{\partial p} - \frac{1}{\beta}
\right)
\end{aligned}
$$

時間を$t$から$t'$に変えると、

$$
\frac{d \zeta}{dt'} = \frac{1}{Q}
\left(
p' \frac{\partial H}{\partial p} - \frac{1}{\beta}
\right)
$$

さて、あらためて$t', p'$を$t,p$と表記すると、運動方程式は

$$
\begin{aligned}
\dot{q} &= \frac{\partial H}{\partial q} \\
\dot{p} &= - \frac{\partial H}{\partial p} - p\zeta \\
\dot{\zeta} &= \frac{1}{Q}\left(
p \frac{\partial H}{\partial p} - \frac{1}{\beta}
\right)
\end{aligned}
$$
となり、これをNose-Hoover法と呼ぶ。この方程式が定常状態として指定の温度のカノニカル分布を持つことは後で証明する。

## 能勢の保存量

先程の運動方程式は$(p,q,\zeta)$で閉じてしまい、$s$に関する式が含まれていなかった。これについて見てみよう。

もともと能勢のハミルトニアンはこのような形であった。

$$
H_\mathrm{Nose} = H(p/s, q) + \frac{p_s^2}{2Q} + \frac{\ln s}{\beta}
$$

この、最後の$\ln s/ \beta$を改めて$\eta$と定義しよう。$\eta$の時間微分は

$$
\begin{aligned}
\frac{d\eta}{dt} &= \frac{1}{\beta s}\frac{ds}{dt} \\
&= \frac{1}{\beta s} \frac{p_s}{Q} \\
&= \frac{\zeta}{\beta s}
\end{aligned}
$$

さらに$t$から$t'$に移ると

$$
\frac{d \eta}{d t'} = \frac{\zeta}{\beta}
$$

これも含めれば運動方程式は、

$$
\begin{aligned}
\dot{q} &= \frac{\partial H}{\partial q} \\
\dot{p} &= - \frac{\partial H}{\partial p} - p\zeta \\
\dot{\zeta} &= \frac{1}{Q}\left(
p \frac{\partial H}{\partial p} - \frac{1}{\beta}
\right) \\
\dot{\eta} &= \frac{\zeta}{\beta}
\end{aligned}
$$

となる。$(p,q,\zeta)$で運動方程式が閉じているので、$\eta$の時間発展を計算する必要は無いが、$\eta$まで考えると、もともとの能勢のハミルトニアン

$$
H_\mathrm{Nose} = H(p/s, q) + \frac{p_s^2}{2Q} + \eta
$$

が時間不変量になっていることがわかる。もともと能勢の方法では、ハミルトンの運動方程式に従うために、能勢のハミルトニアンが保存量となっていたのだが、変数変換を行ってハミルトンの運動方程式でなくなった今でも、これは時間不変量のままとなっている。この量を能勢の保存量、もしくはNose-Hoover保存量と呼ぶことがある。

$\eta$そのものは時間発展には不必要だが、Nose-Hoover保存量を見ることで時間発展の精度を確認するために計算される場合がある。

## Nose-Hoover法の別導出

先程は能勢のハミルトニアンから導出された運動方程式を、変数変換することでNose-Hoover法が導出された。以下では逆に、「定常状態として指定の温度のカノニカル分布が実現するとしたら、運動方程式はどのような形でなければならないか」を考えてみよう。以下、簡単のために一自由度系を考える。

今、カノニカル分布を実現したいハミルトニアン$H_0(p,q)$があるとする。実現したい分布は

$$
f(p,q) \sim \exp(-\beta H)
$$

である。ただし、$\beta = 1/kT$は逆温度である。この位相空間は$(p,q)$で張られている。

さて、この分布を直接実現するのは難しそうなので、自由度$\zeta$を追加し、拡大された位相空間$(p,q,\zeta)$を考える。この空間で、$\zeta$も含めたカノニカル分布

$$
f_\mathrm{ex}(p,q,\zeta) \sim \exp(-\beta H)\exp\left(-\beta \frac{Q \zeta^2}{2}\right)
$$

を考えよう。$Q$の意味は後述する。もしこの分布が実現されたなら、$\zeta$に関して積分してしまうことで、所望の分布$f$を得ることができる。

$$
f_0 = \int_{-\infty}^{\infty} f_\mathrm{ex} d \zeta \sim \exp(-\beta H) 
$$

さて、拡大された位相空間$(p,q,\zeta)$に、先程の分布関数$f$を定常状態に持つような運動方程式を導入したい。ハミルトンの運動方程式の場合には「作用を最小化する」という変分原理から運動方程式が導けたが、温度制御された系にはそのような変分原理は存在しないので、適当に決めることになる。

とりあえずハミルトニアンの運動方程式をなるべく修正しない方向で検討しよう。温度制御のため、運動量$p$と追加自由度$\zeta$の相互作用は必要であろう。そこで、以下のような運動方程式を考えてみる。

$$
\begin{aligned}
\dot{p} &= -\frac{\partial H}{\partial q} - \phi_p(p,\zeta) \\
\dot{q} &= \frac{\partial H}{\partial p} \\
\dot{\zeta} &= \phi_\zeta(p,q,\zeta)
\end{aligned}
$$

我々の目標は、このダイナミクスが拡張された空間でのカノニカル分布$f_\mathrm{ex}$を定常状態に保つように$\phi_p$や$\phi_\zeta$を決めることである。

いま、位相空間が$\Gamma = (p,q,\zeta)$で張られており、そこに速度場$\dot{\Gamma} = (\dot{p},\dot{q},\dot{\zeta})$が定義されているとしよう。この空間の分布関数$f_\mathrm{ex}$を考えると、分布関数と速度場の積$\dot{\Gamma} f_\mathrm{ex}$が流れ場$J$となる。確率の保存則から、分布関数は以下の連続の式を満たす。

$$
\frac{\partial f_\mathrm{ex}}{\partial t} = 
- \mathrm{div} \underbrace{J}_{\dot{\Gamma} f_\mathrm{ex}}
$$

もし$f_\mathrm{ex}$が定常状態なら時間微分がゼロとなるので、

$$
\mathrm{div} \left( \dot{\Gamma} f_\mathrm{ex}\right) = 
  \frac{\partial}{\partial p} (\dot{p} f_\mathrm{ex})
+ \frac{\partial}{\partial q} (\dot{q} f_\mathrm{ex})
+ \frac{\partial}{\partial p} (\dot{\zeta} f_\mathrm{ex}) 
=0
$$

ここで、

$$
\begin{aligned}
\frac{\partial f_\mathrm{ex}}{\partial p} &= - \beta \frac{\partial H}{\partial p} f_\mathrm{ex}\\
\frac{\partial f_\mathrm{ex}}{\partial q} &= -\beta \frac{\partial H}{\partial q} f_\mathrm{ex}\\
\frac{\partial f_\mathrm{ex}}{\partial \zeta} &= - \beta Q \zeta f_\mathrm{ex}\\
\end{aligned}
$$

であることに注意して、一つ一つ愚直に計算していくと、

$$
\begin{aligned}
\frac{\partial}{\partial p} (\dot{p}f_\mathrm{ex}) &=
\left(
-\frac{\partial^2 H}{\partial p \partial q} - \frac{\partial \phi_p}{\partial p} +\beta\frac{\partial H}{\partial p} \frac{\partial H}{\partial q} + \beta \frac{\partial H}{\partial p} \phi_p   
\right) f_\mathrm{ex} \\
\frac{\partial}{\partial q} (\dot{q}f_\mathrm{ex}) &= \left(
\frac{\partial^2 H}{\partial p \partial q}
-\beta\frac{\partial H}{\partial p} \frac{\partial H}{\partial q}
\right)f_\mathrm{ex} \\
\frac{\partial}{\partial \zeta} (\dot{\zeta}f_\mathrm{ex}) &=
\left( \frac{\partial \phi_\zeta}{\partial \zeta} -\beta Q \zeta \phi_\zeta
\right)f_\mathrm{ex}
\end{aligned}
$$

となる。整理すると、

$$
-\frac{\partial \phi_p}{\partial p} + \beta \frac{\partial H}{\partial p} \phi_p
+ \frac{\partial \phi_\zeta}{\partial \zeta} - \beta Q \zeta \phi_\zeta = 0
$$

が満たされなければならない。ここで、ハミルトンの運動方程式由来の項が消えていることに注意しよう。ハミルトンの運動方程式が作る流れは非圧縮であるので、圧縮性流れに寄与しない。

さて、逆に上式が満たされれば、どのような$\phi_p, \phi_\zeta$を与えようとも、カノニカル分布が定常状態に持つような運動方程式を作ることができる。

まず、簡単のために$\phi_\zeta$が$\zeta$に依存しないとしよう。すると

$$
\frac{\partial \phi_\zeta}{\partial \zeta} = 0
$$

となる。次に、$p$と$\zeta$の相互作用を決める$\phi_p$について、$p$と$\zeta$を含む最も簡単な非線形関数である$p \zeta$としてしまおう。すると、満たすべき式は、

$$
-\zeta + \zeta \beta p \frac{\partial H}{\partial p} - \beta Q \zeta \phi_\zeta = 0
$$

となる。$\phi_\zeta$について解くと、

$$
\phi_\zeta = \frac{1}{Q}\left(
p\frac{\partial H}{\partial p} - \frac{1}{\beta}
\right)
$$

以上から、運動方程式は

$$
\begin{aligned}
\dot{p} &= - \frac{\partial H}{\partial q} - p \zeta \\
\dot{q} &= \frac{\partial H}{\partial p} \\
\dot{\zeta} &= \frac{1}{Q}\left(
p\frac{\partial H}{\partial p} - \frac{1}{\beta}
\right)
\end{aligned}
$$

これはNose-Hoover法にほかならない。要するにNose-Hoover法とは、

* 系に一つ自由度$\zeta$を追加し、
* 自由度を追加した世界でのカノニカル分布を実現するように運動方程式を修正したもの

に過ぎない。そこになんらかの物理的な意味を認めるかどうかは、研究者の間で意見が別れている。

## Nose-Hooverの保存量

さて、Nose-Hoover法には、Noseのハミルトニアンに由来する保存量がある。この量を知らないものとして、導出してみよう。

まず、Nose-Hoover法が実現する、拡張された空間におけるカノニカル分布は以下のように書かれる。

$$
f_\mathrm{ex}(p,q,\zeta) \sim \exp(-\beta H)\exp\left(-\beta \frac{Q \zeta^2}{2}\right)
$$

ここで$H$は、もともと我々がカノニカル分布を実現したいハミルトニアンであった。さて、これを見ると、拡張されたハミルトニアン

$$
H_\mathrm{ex} = H + \frac{Q \zeta^2}{2}
$$

に対するカノニカル分布に見える。そこで、この拡張されたハミルトニアンの時間微分を計算してみよう。

$$
\begin{aligned}
\dot{H_\mathrm{ex}} &= \dot{H} + Q \zeta \dot{\zeta} \\
&= \frac{\partial H}{\partial p} \dot{p}
+ \frac{\partial H}{\partial q} \dot{q}
+ \zeta \left( p \frac{\partial H}{\partial p} - \frac{1}{\beta}\right) \\
&= -p\zeta \frac{\partial H}{\partial p} + p \zeta \frac{\partial H}{\partial p} - \frac{\zeta}{\beta}\\
&= - \frac{\zeta}{\beta}
\end{aligned}
$$

ほとんどの項がキャンセルするのだが、最後に少しゴミが残る。そこで、時間微分がこのゴミとキャンセルするような新たな自由度$\eta$を導入しよう。

$$
\dot{\eta} = \frac{\zeta}{\beta}
$$

定義から自明だが、$H_\mathrm{ex} + \eta$は時間保存量となる。これは、Nose-Hoover保存量と一致する。能勢のハミルトニアンからNose-Hoover法を導出した時には、Nose-Hoover保存量は能勢のハミルトニアン由来という意味があったが、「分布関数がカノニカル分布になるべし」という立場からNose-Hoover法を導くと、拡張されたハミルトニアンの時間微分のゴミをキャンセルしただけのように見える。その物理的解釈については読者に委ねる。

TODO: Nose-Hooverの一般化、Nose-Hoover-Family
