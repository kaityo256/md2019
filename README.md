# 分子動力学法の理論と実装

<a href="https://github.com/kaityo256/md2019"> <div class="btn-square"><i class="fab fa-github"></i> View on GitHub</div></a>

## この文書について

これは金沢大学で行われる集中講義の講義ノートにする予定。

[HTML版はこちら](https://kaityo256.github.io/md2019/)

## 内容(予定)

* 本講義の概要と目的
* 分子動力学法の理論的背景
  * 圧力と界面張力
  * 温度制御とエルゴード性
  * 分子動力学法における数値積分法
* 分子動力学法の実装
  * 分子動力学法の実装と基本的アルゴリズム
  * 分子動力学法の高速化手法について
  * 分子動力学法の並列化とプログラム設計

## 本講義の概要と目的

分子動力学法(Moleclar Dynamics method, MD)とは、ニュートンもしくはハミルトンの運動方程式を数値的に解くことで粒子系を時間発展させる手法である。原子や分子、あるいはそれを粗視化した「つぶつぶ」の間にかかる力を計算し、その力によって速度を更新し、その速度によって位置を更新する、というプロセスを繰り返す。MDの基礎方程式が運動方程式という身近なものであることから、なんとなく「とっつきやすい」「わかりやすい」手法のように見える。しかし、MDという手法を突き詰めて考えてみるとだんだんよくわからなくなってくる。ここで温度と言っているのはどういう量なのか？圧力とは何か？そもそもMDにおける時間発展とはなんなのか？

とりあえず分子動力学法の物理的な意味はさておこう。実際に使うにはプログラムを書かなければならない。LAMMPSなど既存のコードを使うにしても、「温度制御」や「圧力制御」が実際に何をやっているのかわからなければ、何か間違った結果を得たとしてもそれに気づくことは難しいであろう。また、完全にブラックボックスにするよりも、裏で使われている高速化手法等について知っておくのは有用であろう。

以上から、本講義は、分子動力学法の理論的側面と実装の両面について扱う。分子動力学法には様々な種類があるが、ここでは短距離相互作用をする古典系のみ考慮する。分子動力学法は長い歴史があり、「枯れた技術」と思われがちだが、中身は意外に奥が深い。本講義を通じて、少しでも分子動力学法は面白いな、と思ってもらえたら幸いである。

## [分子動力学法における観測量](observables/README.md)

* 変数と観測量

## [分子動力学法と圧力](pressure/README.md)

## [分子動力学法における時間積分](integration/README.md)

## [分子動力学法における温度制御](temperature/README.md)
