molecular oriented RMSD with branch and bound

@todo

- 3dsp の失敗する座標を調べる <br>

- Add Usage <br>

- Add benchmark results <br>

@endtodo

# Algorithm

### 分子座標表現

 \( m \) 分子からなる同種分子が \( N \) 個からなる分子集合体の原子座標を考える。
分子座標を

\begin{equation}
  X _ I = \{ x _ {jI} \} _ {j} \in \mathbb{R} ^ {d \times m}
\end{equation}

で表す。ここで \(x _ {jI}\in\mathbb{R}^d\) は
 \(I\) 番目の分子の \(j\) 番目の原子である。
指示ベクトル \(e _ I\) を

\begin{equation}
  e _ I = \{ \delta _ {IK} \} _ K \in \mathbb{R} ^ {1 \times N}
\end{equation}

で定義する。ここで \(\delta _ {IK}\) はクロネッカーのデルタである。
指示ベクトルに対して

\begin{equation}
  e _ Ie _ K^T = \delta _ {IK} \in \mathbb{R} ^ {1 \times 1} \\
  e _ I^Te _ K = \{ \delta _ {jK} \delta _ {Il} \} _ {jl} \in \mathbb{R} ^ {N \times N}
\end{equation}

などが成り立つ。
分子集合体の原子座標は

\begin{equation}
  X = \sum _ {I=1} ^ N e _ I \otimes X _ I \in \mathbb{R} ^ {d \times N m}
\end{equation}

で表わす。
ここで \( \otimes \) はクロネッカー積である。

### 分子対称性補正最小化RMSD

分子集合体座標の組 \( X, Y \in \mathbb{R} ^ {d \times m N} \) に対して、
分子対称性補正RMSDを以下で定義する。

\begin{equation}
  \text{RMSD}(X,Y)=
  \min _ {R, \mu,\nu } \sqrt{\frac{1}{Nm}\sum _ {I=1}^N
  \sum _ {j=1}^m\|x _ {jI}-Ry _ {\mu(j)\nu(I)}\|^2}
\end{equation}

ここで \( R \in \mathbb R ^ {d \times d} \) は \( d \) 次元回転行列、すなわち

\begin{equation}
  RR ^ T = E _ d, \det R = 1
\end{equation}

を満たす。
ここで $E _ n$ は \( n \) 次元単位行列を表す。 $\nu,\mu$ はそれぞれ同一のインデックス集合を始域と終域とする単射 $\nu:\{1,2,\dots,N\}\mapsto\{1,2,\dots,N\}$ , $\mu:\{1,2,\dots,m\}\mapsto\{1,2,\dots,m\}$ である。
考慮する置換の集合を $\mathcal{N}\subset S _ N,\mu\in\mathcal{M}\subset S _ m$ と記載する。
ここで $S _ n$ は \( n \) 次の置換全体の集合である。
RMSDの最小化は二乗和誤差SDの最小化と等価なので、ここからはSDの最小化を考える。

\begin{equation}

  \underset{ R, \nu, \mu } { \min } \text{ \text{SD}}(X, Y, R, \nu, \mu)

  & = & \underset{R, \nu, \mu}{\min} \sum _ {I=  1} ^ N \sum _ {j = 1} ^ m \| x _ {jI} - Ry _ {\mu(j) \nu(I)} \| ^ 2\\

  & = & \underset{R, \nu, \mu}{\min} \sum _ {I = 1} ^ N\sum _ {j = 1} ^ m \left ( \text{tr} \left [ x _ {jI} ^ {\top} x _ {jI} \right ]
    + \text{tr} \left [ y _ {\mu(j) \nu(I)} ^ {\top} R ^ {\top} R y _ {\mu(j) \nu(I)} \right ] \right . \\

  &   & \left . - \text{tr} \left [ x _ {jI} ^ {\top} R y _ {\mu(j) \nu(I)} \right ]
                - \text{tr} \left[y _ {\mu(j) \nu(I)}^{\top} R ^ {\top} x _ {jI} \right ] \right ) \\

  & = & \text{tr} [ X ^ {\top} X ] + \text{tr} [ Y ^ {\top} Y ] - 2 \underset{R, \nu, \mu}{\max}
        \sum _ {I = 1} ^ M \sum _ {j = 1} ^ N \text{tr} [ x _ {jI}^{\top} R y _ {\mu(j) \nu(I)} ] \\

\end{equation}

と変形できるから、問題は右辺第三項の最大化問題となり、

\begin{equation}

  R ^ {\*}, \nu ^ {\*}, \mu ^ {\*}
  = \underset{R, \nu, \mu}{\arg \max} \sum _ {I = 1} ^ N \sum _ {j = 1} ^ m \text{tr} \left [ x _ {jI} ^ {\top} R y _ {\mu(j) \nu(I)} \right ]

\end{equation}

が求まればよいことになる。

### 分子置換

分子集合体構造 $X\in\mathbb{R}^{d\times mN}$ に対して、その分子間置換行列と分子内置換は置換行列 $P\in\mathbb{R}^{N\times N},Q _ I\in\mathbb{R}^{m\times m}$ を用いて

$$
P^{(mN)}=P\otimes E _ M
$$

$$
Q^{(mN)}=\sum _ {I=1}^Ne _ I^Te _ I\otimes Q _ I
$$

で表現できる。さらに、分子内置換 $Q _ I$ の候補は $I$ に依存せず、列挙できるとする。

$$
Q _ I\in \{Q^{(s)};s=1,\dots,S\}
$$

写像 $\sigma:\{1,2,\dots,I\}\mapsto\{1,2,\dots,S\}$ を用いて

$$
Q^{(mN)}=\sum _ {I=1}^N\sum _ {s=1}^S\delta _ {s\sigma(I)}\cdot e _ I^Te _ I\otimes Q^{(s)}
$$

と書き直しておく。任意の分子置換操作は置換の積

$$
\begin{split}
&&P^{(mN)}Q^{(mN)}\\
&=&\left(P\otimes E _ M\right)\left(\sum _ {I=1}^N\sum _ {s=1}^S\delta _ {s\sigma(I)}\cdot e _ I^Te _ I\otimes Q^{(s)}\right)\\
&=&\sum _ {I=1}^N\sum _ {s=1}^S\delta _ {s\sigma(I)}\cdot Pe _ I^Te _ I\otimes Q _ I\\
&=&\sum _ {I=1}^N\sum _ {s=1}^S\delta _ {s\sigma(I)}\cdot P _ I\otimes Q _ I
\end{split}
$$

で表現できる。

一般に $P^{(mN)},Q^{(nN)}$ は非可換で $Q^{(mN)}=E _ N\otimes Q$ のとき、つまり $Q _ 1=Q _ 2=\dots=Q _ N$ のときのみ可換になる。

$$
P^{(MN)}Q^{(MN)}=(P\otimes E _ M)(E _ N\otimes Q)=(E _ N\otimes Q)(P\otimes E _ M)=Q^{(MN)}P^{(MN)}
$$

### トレースの変形

トレースを変形する。

$$
\begin{split}
&\textrm{tr}\left[X^TRYP^{(mN)}Q^{(mN)}\right]\\
=&\textrm{tr}\left[\left(\sum _ {I=1}^Ne_I\otimes X_I\right)^T\left(E _ 1\otimes R\right)\left(\sum _ {J=1}^Ne _ J\otimes Y _ J\right)\left(\sum _ {K=1}^N\sum _ {s=1}^S\delta _ {s\sigma(K)}\cdot P _ K\otimes Q^{(s)}\right)\right]&\leftarrow 定義より\\
=&\sum _ {I=1}^N\sum _ {J=1}^N\sum _ {K=1}^N\textrm{tr}\left[\left(e _ I\otimes X _ I\right)^T\left(E _ 1\otimes R\right)\left(e _ J\otimes Y _ J\right)\left(P _ K\otimes Q^{(s)}\right)\right]&\leftarrowトレースの線形性\\
=&\sum _ {I=1}^N\sum _ {J=1}^N\sum _ {K=1}^N\sum _ {s=1}^S \delta _ {s\sigma(K)} \textrm{tr}\left[e _ I^Te _ JP _ K\otimes X _ I^TRY _ JQ^{(s)}\right]&\leftarrow クロネッカー積の混合積性\\
=&\sum _ {I=1}^N\sum _ {K=1}^N\sum _ {s=1}^S \delta _ {s\sigma(K)}\textrm{tr}\left[\sum _ {J=1}^Np _ {JK}\left(e _ I^Te _ K\otimes X _ I^TRY _ JQ^{(s)}\right)\right]&\leftarrow e _ JP _ K=p _ {JK}e _ K\\
=&\sum _ {I=1}^N\sum _ {K=1}^N\sum _ {s=1}^S \delta _ {s\sigma(K)}\textrm{tr}\left[(e _ I^Te _ K)\otimes \sum _ {J=1}^Np _ {JK}\left(X _ I^TRY _ JQ^{(s)}\right)\right]&\leftarrow クロネッカー積の線形性\\
=&\sum _ {I=1}^N\sum _ {K=1}^N\sum _ {s=1}^S \delta _ {s\sigma(K)}\delta _ {IK}\cdot \textrm{tr}\left[ \sum _ {J=1}^Np _ {JK}X _ I^TRY _ JQ^{(s)}\right]&\leftarrow \textrm{tr}[(e _ I^Te _ K)\otimes A]=\delta _ {IK}\textrm{tr}[A]\\
=&\sum _ {I=1}^N\sum _ {J=1}^N\sum _ {s=1}^S \delta _ {s\sigma(I)}p _ {JI}\cdot\textrm{tr}\left[X _ I^TRY _ JQ^{(s)}\right]&\leftarrow\textrm{tr}[(e _ I^Te _ K)\otimes A]=\delta _ {IK}\textrm{tr}[A]\\
=&\sum _ {I=1}^N\sum _ {J=1}^N\sum _ {s=1}^S \delta _ {s\sigma(I)}p_{JI}\cdot\textrm{tr}\left[RY _ JQ^{(s)}X _ I^T\right]&\leftarrow 行列の線形性とトレースの循環性\\
=&\sum _ {I=1}^N\sum _ {J=1}^N\sum _ {s=1}^S \delta _ {s\sigma(I)}\delta _ {I\nu^{-1}(J)}\cdot\textrm{tr}\left[RY _ JQ^{(s)}X _ I^T\right]&\leftarrow Pは置換行列\\
=&\sum _ {I=1}^N\sum _ {J=1}^N\sum _ {s=1}^S \delta _ {s\sigma(I)}\delta _ {I\nu^{-1}(J)}\cdot\textrm{tr}\left[RC _ {sIJ}\right]&\leftarrow 定義\ C _ {sIJ}:=Y _ JQ^{(s)}X _ I^T\\
=&\sum _ {I=1}^N\textrm{tr}\left[RC _ {\sigma(I)I\nu(I)}\right]\\
\end{split}
$$

ここで登場した共分散行列 $C _ {sIJ}\in\mathbb{R}^{d\times d}$ は分子内置換 $Q^{(s)}$ が既知なら与えられた構造に対して簡単に評価できる。したがって、まず与えられた構造に対して分散共分散テンソルを計算しておいて、その後 $R,\sigma,\nu$ を変分的に最大化すると効率が良い。 

### 変分下限

組み合わせ最適化のため、分枝計画法を用いる。分枝計画法に必要なのは部分集合への分割と部分集合に対する下限（上限）の計算である。そこで、話を再びRMSD値の最小化へと戻し、分割しやすい形へと変形する。評価関数式Xは\( n \)個の定数と行列トレースの和であり、そのインデックス毎に割当 $\sigma,\nu$ を一つずつ重複組合せと順列で選んでいけば、重複なく木構造を作ることができる。

$$
\begin{split}
&&\text{\ RMSD}(X,Y,R,\sigma,\mu)\\
&=&\text{tr}[X^TX]+\text{tr}[Y^TY]-2\sum _ {I=1}^N\textrm{tr}\left[RC _ {\sigma(I)I\nu(I)}\right]\\
&=&\sum _ {I=1}^N\left(\text{tr}[X _ I^TX _ I]+\text{tr}[Y_{\nu(I)}^TY _ {\nu(I)}]-2\text{tr}\left[RC _ {\sigma(I)I\nu(I)}\right]\right)\\
&=&\sum _ {I=1}^N\left(S _ {I\nu(I)}-2\text{tr}[RC _ {\sigma(I)I\nu(I)}]\right)\\
\end{split}
$$

ここで定義した自己相関関数のトレース $S _ {IJ}:=\text{tr}[X  _  I^TX _ I]+\text{tr}[Y _ {J}^TY _ {J}]$ も $C _ {sIJ}$ と同様に事前に計算しておくことができる。

まず、全体の下限を求める。下限の作り方は任意性があるが、ここでは回転行列をインデックス毎に分割することで確認する。

$$
\begin{split}
&&\min _ {R,\sigma,\nu}\sum _ {I=1}^N\left(S _ {I\nu(I)}-2\text{tr}\left[RC _ {\nu}\sum _ {I=1}^N\left(S _ {I\nu(I)}-2\max _ {R _ I,\sigma}\text{tr}\left[R _ IC _ {\sigma(I)I\nu(I)}\right]\right)&\leftarrow 最小値の線形性\\
&=&\min _ {\nu}\sum _ {I=1}^N S _ {I\nu(I)}-2\text{tr}\left[R^b _ IC _ {\sigma^b(I)I\nu(I)}\right]&\leftarrow\arg\maxの代入\\
&=&\min _ {\nu}\sum _ {I=1}^N\left(S _ {I\nu(I)}-2T _ {I\nu(I)}\right)&\leftarrow T _ {IJ}:=\text{tr}\left[R^b _ IC _ {\sigma^b(I)IJ}\right]\\
\end{split}
$$


行列のトレースを最大化する元 $R _ I^b,\sigma^b$ はKabschのアルゴリズム、あるいは四元数を用いたキー行列の最大固有値の推定によって簡単に求まる。また、行列 

$$
  L=\{S _ {IJ}-2T _ {IJ}\} _ {IJ}
$$

を定義してハンガリアンアルゴリズムを用いることによって、これを最小化する組み合わせ $\nu _ 0^{\*}$ は効率的に発見できる。これより

$$
  \text{\ RMSD}(X,Y,R,\sigma,\mu)\ge\sum _ {I=1}^N L _ {I\nu _ 0^*(I)}
$$

が示され、解は必ずこの数値以上であることを保証できる。
これを可能な $J,s$ に対して評価することで、親ノードの下限を引き上げることができる。

# 分枝

$\Sigma$ を $\sigma$ 全体の集合、 $\mathcal{N}$ を $\nu$ 全体の集合（ $=S _ N$ ）として、次の分割を考える。

$$
\Sigma _ s=\{\sigma\in\Sigma;\sigma(1)=s\}
$$

$$
\mathcal{N} _ J=\{\nu\in\mathcal{N};\nu(1)=J\}
$$

このとき $\Sigma=\Sigma _ 1\cup\Sigma _ 2\cup\dots\cup\Sigma _ S, \mathcal{N}=\mathcal{N} _ 1\cup\mathcal{N} _ 2\cup\dots\cup\mathcal{N} _ N$ となっており、それぞれの部分集合に重複はない。 
したがって、これらの組集合 $(\sigma_ _ s,\nu _ J)\in(\Sigma _ s,\mathcal{N} _ J)$ も重複はなく、それらの直積集合の和集合は $\sigma,\nu$ 全体と一致する。
この分割は $SN$ 個 の部分集合を与える。
$\sigma\in\Sigma _ s,\nu\in\mathcal{N} _ J$ に対する評価関数は $\sigma\in\Sigma _ s,\nu\in\mathcal{N} _ J$ に対して

$$
\begin{split}
&\min _ {R,\sigma,\nu} \left(S _ {1J}-2\text{tr}\left[RC _ {s1J}\right]+\sum _ {I=2}^N\left(S _ {I\nu(I)}-2\text{tr}\left[RC _ {\sigma(I)I\nu(I)}\right]\right)\right)\\
\ge& S _ {1J}-2\max _ {R _ 1}\text{tr}\left[R _ 1C _ {s1J}\right]+\min _ {\nu}\sum _ {I=2}^N\left(S _ {I\nu(I)}-2\max _ {R _ I,\sigma}\text{tr}[R _ IC _ {\sigma(I)I\nu(I)}]\right)\\
=&S _ {1J}-2\text{tr}\left[R _ 1^bC _ {s1J}\right]+\sum _ {I=2}^NL _ {I\nu^b _ 1(I)}\\
\end{split}
$$

で与えられる。この右辺第二項もKabschアルゴリズムとハンガリアンアルゴリズムによって簡単に計算できるから、容易に評価可能である。また、この値は

$$
S _ {IJ}-2\text{tr}\left[R _ 1^bC _ {sIJ}\right]+\sum _ {I=2}^NL _ {I\nu^b _ 1(I)}=L _ {IJ}+\sum _ {I=2}^NL _ {I\nu^b _ 1(I)}=\sum _ {I=1}^NL _ {I\nu^b _ 1(I)}\ge\sum _ {I=1}^NL _ {I\nu^b _ 0(I)}
$$

だから、部分集合の下限は集合全体の下限よりも大きな値を持つことが確認できる。

同様に 

$$
\Sigma _ s=\Sigma _ {s1}\cup\Sigma _ {s2}\cup\dots\cup\Sigma _ {sS}
$$

$$
\mathcal{N} _ {J}=\mathcal{N} _ {J1}\cup\mathcal{N} _ {J2}\cup\dots\cup\mathcal{N} _ {J(N-1)}
$$

となる分割を考えて、下限を $\sigma\in\Sigma _ {ss'},\nu\in\mathcal{N} _ {JJ'}$ に対して

$$
\begin{split}
&&\min _ {R,\sigma,\nu}\left(S _ {1J}+S _ {2J'}-2\text{tr}\left[RC _ {s1J}\right]-2\text{tr}\left[RC _ {s'2J'}\right]+\sum _ {I=3}^N\left(S _ {I\nu(I)}-2\text{tr}\left[RC _ {\sigma(I)I\nu(I)}\right]\right)\right)\\
&\ge&S _ {1J}+S _ {2J'}-2\max _ {R _ _ {s1J}+C_{s'2J'}\right)\right\]+\min _ {\nu\in\mathcal{N}\ _ {JJ'}}\sum\ _ {I=3}^N\left(S _ {I\nu(I)}-2\max _ {R _ _ IC\ _ {\sigma(I)I\nu(I)}\right\]\right)\\
&=&S_ _ {1J}+S _ {2J'}-2\text{tr}\left[R _ 2^b\left(C _ {s1J}+C _ {s'2J'}\right)\right]+\sum _ {I=3}^NL _ {I\nu^b _ 2(I)}\\
\end{split}
$$



と定めれば、必ず親集合よりも大きな値を持つことが確認できる。

## 限定法


## RMSD最小化のメモ


$A,B\in\mathbb {R}^{d\times n},\mathbb{R}\in\mathbb{R}^{^{d\times d}},P\in\{0,1\}^{n\times n}$ とする。$R^TR=I,P^TP=I$である。一般化RMSDは $\textrm{tr}(RBPA^T)$の最大化問題に帰着される。

コーシーシュワルツの不等式の一般化により
$$
\textrm{tr}(RBPA^T)^2\le\textrm{tr}(R^TR)\textrm{tr}((BPA^T)^TBPA^T)\\
=d\cdot\textrm{tr}\left(S^TS\right)\\
$$
がなりたつ。ここで$S=BPA^T$とおいた。

※これ以上分解すると
$$
\textrm{tr}(P^T\Sigma _ bP\Sigma _ A)\le\sqrt{\textrm{tr}(P^T\Sigma _ b^TPP^T\Sigma _ bP)\textrm{tr}(\Sigma _ A^T\Sigma _ A)}\\
=\textrm{tr}(\Sigma _ b)\textrm{tr}(\Sigma _ A)\\
=\left(\sum _ {i=1}^n\lambda _ i^{(A)}\right)\left(\sum _ {i=1}^n\lambda _ i^{(B)}\right)
$$
となり自明な情報しか得られない。

ここで正規直交行列$\tilde P\in\mathbb{R}^{n\times n}$が$\tilde P^T\tilde P=I$を満たすとする。このとき、Kabschのアルゴリズムを用いて、特異値分解
$$
A^TB=U\Sigma V^T
$$
に対して
$$
\tilde P=VX U^T
$$
が定まる（はずである、ok）。ここで$X=\mathrm{Diag}(1,1,\dots,1,\det(VU^T))$ である。

$P$は部分空間に分けられるので分枝限定法に利用できるはずである。
$$
P=\begin{pmatrix}{}
P _ {1} &\mathbf 0\\
\mathbf 0 &P _ {2}\\
\end{pmatrix}
$$
とする。このとき...

Rを固定してPの変分下界を計算することが必要？





## 変分上界

最小化RMSDを見つけるのは難しい。そこでまず $\mathrm{tr}[X^TRY\tilde S]$ の変分上界を評価することにする。 $\mathcal O _ n$ を \( n \) 次元直交行列全体の集合として、$P\in\mathcal O _ N,P _ I=e^T _ Ie _ IP,Q _ I\in\mathcal O _ M$ に対して
$$
S=\sum _ {I=1}^NP _ I\otimes Q _ I\\
$$
と定義する。このとき、$\mathcal S _ n\subset\mathcal O _ n$から、
$$
\underset{\tilde S\in\mathcal S _ N\times(\sigma _ M)^N,R\in\mathcal{SO _ d}}{\max\arg}\mathrm{tr}[X^TRY\tilde S]\le\underset{S\in\mathcal O _ N\times(\mathcal O _ M)^N,R\in\mathcal{SO _ d}}{\max\arg}\mathrm{tr}[X^TRYS]
$$
が成り立つ。

$P^{(MN)}=P\otimes E _ M$ は任意の分子置換、$Q^{(MN)}=\sum _ {I=1}^Ne _ I^Te _ IQ _ I$ は任意の分子内回転を表す置換行列を部分集合に含む直交行列である。したがって $P^{(MN)}Q^{(MN)}$は任意の分子対称性を制限された置換行列を部分集合に含む。一般に $P^{(MN)},Q^{(MN)}$ は非可換で$Q^{(MN)}=E _ N\otimes Q$ のとき、つまり $Q _ 1=Q _ 2=\dots=Q _ N$ のときのみ可換になる（$P^{(MN)}Q^{(MN)}=(P\otimes E _ M)(E _ N\otimes Q)=(E _ N\otimes Q)(P\otimes E _ M)=Q^{(MN)}P^{(MN)}$）。

トレースを変形する。
$$
\begin{split}
&\textrm{tr}[X^TRYS]\\
=&\textrm{tr}\left[\left(\sum _ {I=1}^Ne _ I\otimes X _ I\right)^T\left(E _ 1\otimes R\right)\left(\sum _ {J=1}^Ne _ J\otimes Y _ J\right)\left(\sum _ {K=1}^NP _ K\otimes Q _ K\right)\right]\leftarrow 定義より\\
=&\sum _ {I=1}^N\sum _ {J=1}^N\sum _ {K=1}^N\textrm{tr}\left[\left(e _ I\otimes X _ I\right)^T\left(E _ 1\otimes R\right)\left(e _ J\otimes Y _ J\right)\left(P _ K\otimes Q _ K\right)\right]\leftarrowトレースの線形性\\
=&\sum _ {I=1}^N\sum _ {J=1}^N\sum _ {K=1}^N\textrm{tr}\left[e _ I^Te _ JP _ K\otimes X _ I^TRY _ JQ _ K\right]\leftarrow クロネッカー積の混合積性\\
=&\sum _ {I=1}^N\sum _ {K=1}^N\textrm{tr}\left[\sum _ {J=1}^Np _ {JK}\left(e _ I^Te _ K\otimes C _ {IJ}Q _ K\right)\right]\leftarrow e _ JP _ K=p _ {JK}e _ K\\
=&\sum _ {I=1}^N\sum _ {K=1}^N\textrm{tr}\left[(e _ I^Te _ K)\otimes \sum _ {J=1}^Np _ {JK}(C _ {IJ}Q _ K)\right]\leftarrow クロネッカー積の線形性\\
=&\sum _ {I=1}^N\textrm{tr}\left[\sum _ {J=1}^Np _ {JI}(C _ {IJ}Q _ I)\right]\leftarrow\textrm{tr}[(e _ I^Te _ J)\otimes A]=\delta _ {IJ}\textrm{tr}[A]\\
=&\sum _ {I=1}^N\textrm{tr}\left[Q _ I\sum _ {J=1}^Np _ {JI}C _ {IJ}\right]\leftarrow 行列の線形性とトレースの循環性\\
\end{split}
$$
ここで
$$
C _ {IJ}=X^T _ IRY _ J\in\mathbb{R}^{M\times M}
$$
を定義した。$Q _ I,C _ {IJ}\in\mathbb{R}^{M\times M}$ であり、その行列積は\( n \)回評価すればよい（$\mathcal O(NM^3)$）。これは$NM\times NM$の行列積（$\mathcal O(N^3M^3)$）を計算するよりははるかに評価が簡単になる。

## 次元削減

$P,Q$を同時に解くのは難しいが、行列
$$
C _ {Q _ I}(R,P)=\sum _ {J=1}^Np _ {JI}C _ {IJ}
$$
$$
C _ P(R,Q)=\left\{\mathrm{tr}\left[Q _ I\sum _ {J=1}^N C _ {IJ}\right]\right\} _ {IJ}
$$

に対して直交プロクラステス問題を繰り返し解くことで、$P,Q$を自己無撞着に求めることができる。ここで $d\le n,d\le m$ の場合、$\sum _ {J=1}^Np _ {JI}C _ {IJ}$ や $\{\mathrm{tr}[Q _ I\sum _ {J=1}^N C _ {IJ}]\} _ {IJ}$ の持つ次元はたかだか \( d \) に制限されるため、それぞれフル次元で求解することは非効率で、解の不定性と不安定性をもたらす。そこで事前処理として、$X,Y$の有効次元を落としておく。
直交行列 $U\in\mathbb{R}^{d\times d},V\in\mathbb{R}^{n\times n},W\in\mathbb{R}^{m\times m}$ に対して
$$
X'=UX(V\otimes W)
$$
$$
Y'=UY(V\otimes W)
$$
とする。これはトレースの循環性より評価値を変えない。ここで $X,Y\in\mathbb R^{d\times n}$ に対して、$Z=X-Y$ とおくと、$U$ は $Z$ の左特異ベクトル、$V\otimes W$ は右特異ベクトルを取れば良い。
$$
U^TZT
$$

分散が0の軸$w\in\mathbb R^{n\times 1}$は
$$
w^TZ^TZw=0
$$
を満たす。したがって
$$
0=w^T(X^TX+Y^TY+X^TY+Y^TX)w=x^Tx+y^Ty+x^Ty+y^Tx
$$
$$
x^Ty=\frac{1}{2}(x^Tx+y^Ty)
$$
が成り立っている。


そこでまずPCAによる次元削減を行い、\( d \) 次元以下の直交プロクラステス問題を解くことを考える。

まずは $C _ Q$ の次元削減を考える。
$$
Z _ {IJ}=X _ {I}-RY _ {J}\in\mathbb{R}^{d\times m}
$$
に対して、対称行列
$$
S _ {IJ}=Z _ {IJ}^TZ _ {IJ}
$$
を作る。このとき固有ベクトルとなる直交行列 $W _ {IJ}$ を取ることができる。
$$
\left(\sum _ {J=1}^NW^T _ {IJ}Z _ {IJ}^TZ _ {IJ}W _ {IJ}\right)=\sum _ {J=1}^N\Lambda _ {IJ}
$$
実効次元が制限されるため、任意の$IJ$と$k>d$ に対して、$\lambda _ {IJ,1}=\dots=\lambda _ {IJ,M}=0$ が成り立つように取れば、任意の$\lambda _ k=0$ に対応する固有ベクトル $w _ k$ は
$$
X _ Iw=RY _ Jw
$$

$$
w _ k^T(X _ I^TX _ I+Y _ J^TY _ J-X^T _ IRY _ J-Y^T _ JR^TX _ I)w _ k=0
$$

を満たすから、式変形すると
$$
w _ k^TC _ {IJ}w _ k=w _ k^TX _ I^TX _ Iw _ k
$$
となり、直交プロクラステス問題は
$$
\mathrm{tr}[Q _ IC _ {P _ I}]\\
=\mathrm{tr}\left[\sum _ {J=1}^Np _ {JI}C _ {IJ}\right]\\
=\mathrm{tr}[W^TQ _ IWW^TC _ {IJ}W]\\
=\mathrm{tr}\left[Q _ I'W'^T(X _ I^TRY _ J)W'\right]+\mathrm{tr}\left[Q _ I''W''^T(X _ I^TRY _ J)W''\right]\\
=\mathrm{tr}\left[Q' _ I(W'^TX _ I^T)(RY _ JW')\right]+\mathrm{tr}\left[Q _ I''X _ I^TX _ I\right]\\
$$
とすることができる。ここで右辺第二項を最大化する直交行列は $Q _ I''=E _ {m-d}$ であると直ちに定まるから、$Q _ I'\in\mathbb{R}^{d\times d}$ を求める問題に帰着する。


$$
W'^TC _ {IJ}W'
$$


## 目的関数の勾配

ラグランジュの未定乗数 $\lambda=(\lambda _ P,\lambda _ Q,\lambda _ R,\lambda _ {|R|})$ に対して
$$
f(P,Q,R,\lambda)=\sum _ {I=1}^N\mathrm{tr}\left[Q _ I\sum _ {J=1}^Np _ {JI}X^T _ IRY _ J\right]-\lambda _ Pg(P)-\lambda _ Qg(Q)-\lambda _ Rg(R)-\lambda _ {|R|}(\det(R)-1)^2
$$
を目的関数として最大化することを考える。ここで $g$ は$A\in\mathbb{R}^{n\times m}$ に対して
$$
g(A)=\|AA^T-E _ n\|^2=\mathrm{tr}[(AA^T)^2]-2\mathrm{tr}[AA^T]+n
$$
で定義される。

$P,Q$ の拘束条件より
$$
f(P,Q _ 1,Q _ 2,\dots,Q _ N,R,\lambda)\\
=\mathrm{tr}\left[\sum _ {I=1}^NQ _ I\sum _ {J=1}^Np _ {JI}C _ {IJ}\right]-\lambda _ Pg(P)-\lambda _ Q\sum _ {I=1}^Ng(Q _ I)-\lambda _ Rg(R)-\lambda _ {|R|}(\det(R)-1)^2
$$
の最大化問題へと簡略化される。

$P$ について偏微分を行うと
$$
\frac{\partial f}{\partial P}=\frac{\partial}{\partial P}\sum _ {I=1}^N\mathrm{tr}\left[Q _ I\sum _ {J=1}^Np _ {JI}C _ {IJ}\right]-\lambda _ P\frac{\partial g(P)}{\partial P}\\
$$
と計算できる。

第一項は
$$
\frac{\partial}{\partial P}\sum _ {I=1}^N\mathrm{tr}\left[Q _ I\sum _ {J=1}^Np _ {JI}C _ {IJ}\right]=\left\{\frac{\partial}{\partial p _ {mn}}\sum _ {I=1}^N\mathrm{tr}\left[Q _ I\sum _ {J=1}^Np _ {JI}C _ {IJ}\right]\right\} _ {mn}
=\left\{\mathrm{tr}\left[Q _ nC _ {nm}\right]\right\} _ {mn}\\
$$
と求まる。第二項は
$$
\lambda _ P\frac{\partial g(P)}{\partial P}=4\lambda _ P((\mathrm{Tr}[(P^TP)^2])'-P)=4\lambda _ P(PP^T-I)P
$$
と計算でき、ここで、$\mathbb{R}^{n\times n}\ni A=\{a _ {ij}\} _ {ij}$ に対して
$$
(\mathrm{Tr}[(A^TA)^2])'=\left\{a _ {ij}^3+\left(\sum _ {k\ne j}a _ {ik}^2+\sum _ {k\ne i}a _ {kj}^2\right)a _ {ij}+\sum _ {k\ne i,l\ne j}a _ {il}a _ {kj}a _ {kl}\right\} _ {ij}\\
=\left\{\sum _ {k=1}^N\sum _ {l=1}^Na _ {il}a _ {kj}a _ {kl}\right\} _ {ij}\\
=AA^TA
$$
まとめると
$$
\frac{\partial f}{\partial P}=\left\{\mathrm{tr}\left[Q _ nC _ {nm}\right]\right\} _ {mn}-4\lambda _ P(PP^T-I)P\\
$$
が得られる。

$Q _ I$についても同様に
$$
\frac{\partial}{\partial Q _ I}\sum _ {I=1}^N\mathrm{tr}\left[Q _ I\sum _ {J=1}^Np _ {JI}C _ {IJ}\right]=\sum _ {J=1}^Np _ {JI}C _ {IJ}^T\\
$$
であるから
$$
\frac{\partial f}{\partial Q _ {J=1}^Np_ _ {JI}C _ {IJ}^T-4\lambda _ Q(Q _ IQ _ I^T-E)Q _ I\\
$$
と定まる。

$R$については
$$
\sum _ {I=1}^N\textrm{tr}\left[Q _ I\sum _ {J=1}^Np _ {JI}C _ {IJ}\right]=\sum _ {I=1}^N\sum _ {J=1}^Np _ {JI}\textrm{tr}\left[Q _ IX^T _ IRY _ J\right]\\
=\sum _ {I=1}^N\sum _ {J=1}^Np _ {JI}\textrm{tr}\left[Y _ JQ _ IX^T _ IR\right]\\
$$
と変形しておいて
$$
\frac{\partial}{\partial R}\sum _ {I=1}^N\sum _ {J=1}^Np _ {JI}\textrm{tr}\left[Y _ JQ _ IX^T _ IR\right] =\sum _ {I=1}^N\sum _ {J=1}^Np _ {JI}\left(Y _ JQ _ IX^T _ I\right)^T\\
=\sum _ {I=1}^N\left(X _ IQ _ I^T\left(\sum _ {J=1}^Np _ {JI}Y _ J^T\right)\right)
$$
と計算できる。この形式は部分回転を考える際に計算不要な部分をスキップできるためクロネッカー積で表される$\mathbb{R}^{MN\times MN}$ の直行行列を用いるより有利である。まとめて
$$
\frac{\partial f}{\partial R}=\sum _ {I=1}^N\left(X _ IQ _ I^T\left(\sum _ {J=1}^Np _ {JI}Y _ J^T\right)\right)-4\lambda _ R(RR^T-E)R-2\lambda _ {|R|}(\det(R)-1)\mathrm{adj}(R)\\
$$
ここで$\mathrm{adj}(R)$ は $R$ の余因子行列を表す。

未定乗数の勾配についても同様に
$$
\frac{\partial f}{\partial \lambda _ P}=g(P)\\
$$
$$
\frac{\partial f}{\partial \lambda _ Q}=\sum _ {I=1}^Ng(Q _ I)\\
$$
$$
\frac{\partial f}{\partial \lambda _ R}=(R)\\
$$
$$
\frac{\partial f}{\partial \lambda _ {|R|}}=(\det(R)-1)^2\\
$$

と計算される。
