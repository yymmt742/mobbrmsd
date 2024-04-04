# 一般化RMSD-fitを枝刈りする

## 分子座標表現

$m$ 分子からなる同種分子が $N$ 個からなる分子集合体の原子座標を考える。分子座標を

$$
X_I=\\{x_{jI}\\}_{j}\in\mathbb{R}^{d\times m}
$$

で表す。ここで $x_{jI}\in\mathbb{R}^d$は $I$ 番目の分子の $j$ 番目の原子である。指示ベクトル $e_I$ を

$$
e_I=\\{\delta_{IK}\\}_K\in\mathbb{R}^{1\times N}
$$

で定義する。ここで $\delta_{IK}$ はクロネッカーのデルタである。指示ベクトルに対して

$$
e_Ie_K^T=\delta_{IK}\in\mathbb{R}^{1\times1}\\
e_I^Te_K=\{\delta_{jK}\delta_{Il}\}_{jl}\in\mathbb{R}^{N\times N}
$$

などが成り立つ。分子集合体の原子座標は

$$
X=\sum_{I=1}^Ne_I\otimes X_I\in\mathbb{R}^{d\times Nm}
$$

で表わす。ここで $\otimes$ はクロネッカー積である。

## 分子対称性補正最小化RMSD

分子集合体座標の組 $X,Y\in\mathbb{R}^{d\times mN}$ に対して、分子対称性補正RMSDを以下で定義する。

$$
\text{RMSD}(X,Y)=\\
\min_{R, \mu,\nu } \sqrt{\frac{1}{Nm}\sum_{I=1}^N\sum_{j=1}^m\|x_{jI}-Ry_{\mu(j)\nu(I)}\|^2}
$$

ここで$R\in\mathbb R^{d\times d}$は $d$ 次元回転行列、すなわち

$$
RR^T=E_d, \det R=1
$$

を満たす。
ここで $E_n$ は $n$ 次元単位行列を表す。 $\nu,\mu$ はそれぞれ同一のインデックス集合を始域と終域とする単射 $\nu:\{1,2,\dots,N\}\mapsto\{1,2,\dots,N\}$ , $\mu:\{1,2,\dots,m\}\mapsto\{1,2,\dots,m\}$ である。
考慮する置換の集合を $\mathcal{N}\subset S_N,\mu\in\mathcal{M}\subset S_m$ と記載する。
ここで $S_n$ は $n$ 次の置換全体の集合である。
RMSDの最小化は二乗和誤差SDの最小化と等価なので、ここからはSDの最小化を考える。

$$
\begin{split}
\underset{R,\nu,\mu}{\min}\text{\ SD}(X,Y,R,\nu,\mu)
&=&\underset{R,\nu,\mu}{\min}\sum_{I=1}^N\sum_{j=1}^m\|x_{jI}-Ry_{\mu(j)\nu(I)}\|^2\\
&=&\underset{R,\nu,\mu}{\min}\sum_{I=1}^N\sum_{j=1}^m\left(\text{tr}\left[x_{jI}^Tx_{jI}\right]+\text{tr}\left[y_{\mu(j)\nu(I)}^TR^TRy_{\mu(j)\nu(I)}\right]\right.\\
&&\left.-\text{tr}\left[x_{jI}^TRy_{\mu(j)\nu(I)}\right]-\text{tr}\left[y_{\mu(j)\nu(I)}^TR^Tx_{jI}\right]\right)\\
&=&\text{tr}[X^TX]+\text{tr}[Y^TY]-2\underset{R,\nu,\mu}{\max}\sum_{I=1}^M\sum_{j=1}^N\text{tr}[x_{jI}^TRy_{\mu(j)\nu(I)}]\\
\end{split}
$$

と変形できるから、問題は右辺第三項の最大化問題となり、

$$
R^{\*},\nu^{\*},\mu^{\*} = \underset{R,\nu,\mu}{\arg\max} \sum_{I=1}^N \sum_{j=1}^m \text{tr}\left[ x_{jI}^TRy_{\mu(j)\nu(I)} \right]
$$

が求まればよいことになる。

## 分子置換

分子集合体構造 $X\in\mathbb{R}^{d\times mN}$ に対して、その分子間置換行列と分子内置換は置換行列 $P\in\mathbb{R}^{N\times N},Q_I\in\mathbb{R}^{m\times m}$ を用いて

$$
P^{(mN)}=P\otimes E_M
$$

$$
Q^{(mN)}=\sum_{I=1}^Ne_I^Te_I\otimes Q_I
$$

で表現できる。さらに、分子内置換 $Q_I$ の候補は $I$ に依存せず、列挙できるとする。

$$
Q_I\in \{Q^{(s)};s=1,\dots,S\}
$$

写像 $\sigma:\{1,2,\dots,I\}\mapsto\{1,2,\dots,S\}$ を用いて

$$
Q^{(mN)}=\sum_{I=1}^N\sum_{s=1}^S\delta_{s\sigma(I)}\cdot e_I^Te_I\otimes Q^{(s)}
$$

と書き直しておく。任意の分子置換操作は置換の積

$$
\begin{split}
&&P^{(mN)}Q^{(mN)}\\
&=&\left(P\otimes E_M\right)\left(\sum_{I=1}^N\sum_{s=1}^S\delta_{s\sigma(I)}\cdot e_I^Te_I\otimes Q^{(s)}\right)\\
&=&\sum_{I=1}^N\sum_{s=1}^S\delta_{s\sigma(I)}\cdot Pe_I^Te_I\otimes Q_I\\
&=&\sum_{I=1}^N\sum_{s=1}^S\delta_{s\sigma(I)}\cdot P_I\otimes Q_I
\end{split}
$$

で表現できる。

一般に $P^{(mN)},Q^{(nN)}$ は非可換で $Q^{(mN)}=E_N\otimes Q$ のとき、つまり $Q_1=Q_2=\dots=Q_N$ のときのみ可換になる。

$$
P^{(MN)}Q^{(MN)}=(P\otimes E_M)(E_N\otimes Q)=(E_N\otimes Q)(P\otimes E_M)=Q^{(MN)}P^{(MN)}
$$

## トレースの変形

トレースを変形する。

$$
\begin{split}
&\textrm{tr}\left[X^TRYP^{(mN)}Q^{(mN)}\right]\\
=&\textrm{tr}\left[\left(\sum_{I=1}^Ne_I\otimes X_I\right)^T\left(E_1\otimes R\right)\left(\sum_{J=1}^Ne_J\otimes Y_J\right)\left(\sum_{K=1}^N\sum_{s=1}^S\delta_{s\sigma(K)}\cdot P_K\otimes Q^{(s)}\right)\right]&\leftarrow 定義より\\
=&\sum_{I=1}^N\sum_{J=1}^N\sum_{K=1}^N\textrm{tr}\left[\left(e_I\otimes X_I\right)^T\left(E_1\otimes R\right)\left(e_J\otimes Y_J\right)\left(P_K\otimes Q^{(s)}\right)\right]&\leftarrowトレースの線形性\\
=&\sum_{I=1}^N\sum_{J=1}^N\sum_{K=1}^N\sum_{s=1}^S \delta_{s\sigma(K)} \textrm{tr}\left[e_I^Te_JP_K\otimes X_I^TRY_JQ^{(s)}\right]&\leftarrow クロネッカー積の混合積性\\
=&\sum_{I=1}^N\sum_{K=1}^N\sum_{s=1}^S \delta_{s\sigma(K)}\textrm{tr}\left[\sum_{J=1}^Np_{JK}\left(e_I^Te_K\otimes X_I^TRY_JQ^{(s)}\right)\right]&\leftarrow e_JP_K=p_{JK}e_K\\
=&\sum_{I=1}^N\sum_{K=1}^N\sum_{s=1}^S \delta_{s\sigma(K)}\textrm{tr}\left[(e_I^Te_K)\otimes \sum_{J=1}^Np_{JK}\left(X_I^TRY_JQ^{(s)}\right)\right]&\leftarrow クロネッカー積の線形性\\
=&\sum_{I=1}^N\sum_{K=1}^N\sum_{s=1}^S \delta_{s\sigma(K)}\delta_{IK}\cdot \textrm{tr}\left[ \sum_{J=1}^Np_{JK}X_I^TRY_JQ^{(s)}\right]&\leftarrow \textrm{tr}[(e_I^Te_K)\otimes A]=\delta_{IK}\textrm{tr}[A]\\
=&\sum_{I=1}^N\sum_{J=1}^N\sum_{s=1}^S \delta_{s\sigma(I)}p_{JI}\cdot\textrm{tr}\left[X_I^TRY_JQ^{(s)}\right]&\leftarrow\textrm{tr}[(e_I^Te_K)\otimes A]=\delta_{IK}\textrm{tr}[A]\\
=&\sum_{I=1}^N\sum_{J=1}^N\sum_{s=1}^S \delta_{s\sigma(I)}p_{JI}\cdot\textrm{tr}\left[RY_JQ^{(s)}X_I^T\right]&\leftarrow 行列の線形性とトレースの循環性\\
=&\sum_{I=1}^N\sum_{J=1}^N\sum_{s=1}^S \delta_{s\sigma(I)}\delta_{I\nu^{-1}(J)}\cdot\textrm{tr}\left[RY_JQ^{(s)}X_I^T\right]&\leftarrow Pは置換行列\\
=&\sum_{I=1}^N\sum_{J=1}^N\sum_{s=1}^S \delta_{s\sigma(I)}\delta_{I\nu^{-1}(J)}\cdot\textrm{tr}\left[RC_{sIJ}\right]&\leftarrow 定義\ C_{sIJ}:=Y_JQ^{(s)}X_I^T\\
=&\sum_{I=1}^N\textrm{tr}\left[RC_{\sigma(I)I\nu(I)}\right]\\
\end{split}
$$

ここで登場した共分散行列 $C_{sIJ}\in\mathbb{R}^{d\times d}$ は分子内置換 $Q^{(s)}$ が既知なら与えられた構造に対して簡単に評価できる。したがって、まず与えられた構造に対して分散共分散テンソルを計算しておいて、その後 $R,\sigma,\nu$ を変分的に最大化すると効率が良い。 

## 変分下限

組み合わせ最適化のため、分枝計画法を用いる。分枝計画法に必要なのは部分集合への分割と部分集合に対する下限（上限）の計算である。そこで、話を再びRMSD値の最小化へと戻し、分割しやすい形へと変形する。評価関数式Xは$N$個の定数と行列トレースの和であり、そのインデックス毎に割当 $\sigma,\nu$ を一つずつ重複組合せと順列で選んでいけば、重複なく木構造を作ることができる。

$$
\begin{split}
&&\text{\ RMSD}(X,Y,R,\sigma,\mu)\\
&=&\text{tr}[X^TX]+\text{tr}[Y^TY]-2\sum_{I=1}^N\textrm{tr}\left[RC_{\sigma(I)I\nu(I)}\right]\\
&=&\sum_{I=1}^N\left(\text{tr}[X_I^TX_I]+\text{tr}[Y_{\nu(I)}^TY_{\nu(I)}]-2\text{tr}\left[RC_{\sigma(I)I\nu(I)}\right]\right)\\
&=&\sum_{I=1}^N\left(S_{I\nu(I)}-2\text{tr}[RC_{\sigma(I)I\nu(I)}]\right)\\
\end{split}
$$

ここで定義した自己相関関数のトレース $S_{IJ}:=\text{tr}[X_I^TX_I]+\text{tr}[Y_{J}^TY_{J}]$ も $C_{sIJ}$ と同様に事前に計算しておくことができる。

まず、全体の下限を求める。下限の作り方は任意性があるが、ここでは回転行列をインデックス毎に分割することで確認する。

$$
\begin{split}
&&\min_{R,\sigma,\nu}\sum_{I=1}^N\left(S_{I\nu(I)}-2\text{tr}\left[RC_{\sigma(I)I\nu(I)}\right]\right)&\\
&=&\min_{\nu}\sum_{I=1}^N\left(S_{I\nu(I)}-2\max_{R_I,\sigma}\text{tr}\left[R_IC_{\sigma(I)I\nu(I)}\right]\right)&\leftarrow 最小値の線形性\\
&=&\min_{\nu}\sum_{I=1}^N S_{I\nu(I)}-2\text{tr}\left[R^b_IC_{\sigma^b(I)I\nu(I)}\right]&\leftarrow\arg\maxの代入\\
&=&\min_{\nu}\sum_{I=1}^N\left(S_{I\nu(I)}-2T_{I\nu(I)}\right)&\leftarrow T_{IJ}:=\text{tr}\left[R^b_IC_{\sigma^b(I)IJ}\right]\\
\end{split}
$$


行列のトレースを最大化する元 $R_I^b,\sigma^b$ はKabschのアルゴリズム、あるいは四元数を用いたキー行列の最大固有値の推定によって簡単に求まる。また、行列 

$$
L=\{S_{IJ}-2T_{IJ}\}_{IJ}
$$

を定義してハンガリアンアルゴリズムを用いることによって、これを最小化する組み合わせ $\nu_0^{\*}$ は効率的に発見できる。これより

$$
\text{\ RMSD}(X,Y,R,\sigma,\mu)\ge\sum_{I=1}^N L_{I\nu_0^*(I)}
$$

が示され、解は必ずこの数値以上であることを保証できる。
これを可能な $J,s$ に対して評価することで、親ノードの下限を引き上げることができる。

## 分枝

$\Sigma$ を $\sigma$ 全体の集合、 $\mathcal{N}$ を $\nu$ 全体の集合（ $=S_N$ ）として、次の分割を考える。

$$
\Sigma_s=\{\sigma\in\Sigma;\sigma(1)=s\}
$$

$$
\mathcal{N}_J=\{\nu\in\mathcal{N};\nu(1)=J\}
$$

このとき $\Sigma=\Sigma_1\cup\Sigma_2\cup\dots\cup\Sigma_S, \mathcal{N}=\mathcal{N}_1\cup\mathcal{N}_2\cup\dots\cup\mathcal{N}_N$ となっており、それぞれの部分集合に重複はない。 
したがって、これらの組集合 $(\sigma_s,\nu_J)\in(\Sigma_s,\mathcal{N}_J)$ も重複はなく、それらの直積集合の和集合は $\sigma,\nu$ 全体と一致する。
この分割は $SN$ 個 の部分集合を与える。
$\sigma\in\Sigma_s,\nu\in\mathcal{N}_J$ に対する評価関数は $\sigma\in\Sigma_s,\nu\in\mathcal{N}_J$ に対して

$$
\begin{split}
&\min_{R,\sigma,\nu} \left(S_{1J}-2\text{tr}\left[RC_{s1J}\right]+\sum_{I=2}^N\left(S_{I\nu(I)}-2\text{tr}\left[RC_{\sigma(I)I\nu(I)}\right]\right)\right)\\
\ge& S_{1J}-2\max_{R_1}\text{tr}\left[R_1C_{s1J}\right]+\min_{\nu}\sum_{I=2}^N\left(S_{I\nu(I)}-2\max_{R_I,\sigma}\text{tr}[R_IC_{\sigma(I)I\nu(I)}]\right)\\
=&S_{1J}-2\text{tr}\left[R_1^bC_{s1J}\right]+\sum_{I=2}^NL_{I\nu^b_1(I)}\\
\end{split}
$$

で与えられる。この右辺第二項もKabschアルゴリズムとハンガリアンアルゴリズムによって簡単に計算できるから、容易に評価可能である。また、この値は

$$
S_{IJ}-2\text{tr}\left[R_1^bC_{sIJ}\right]+\sum_{I=2}^NL_{I\nu^b_1(I)}=L_{IJ}+\sum_{I=2}^NL_{I\nu^b_1(I)}=\sum_{I=1}^NL_{I\nu^b_1(I)}\ge\sum_{I=1}^NL_{I\nu^b_0(I)}
$$

だから、部分集合の下限は集合全体の下限よりも大きな値を持つことが確認できる。

同様に 

$$
\Sigma_s=\Sigma_{s1}\cup\Sigma_{s2}\cup\dots\cup\Sigma_{sS}
$$

$$
\mathcal{N}\_{J}=\mathcal{N}\_{J1}\cup\mathcal{N}\_{J2}\cup\dots\cup\mathcal{N}\_{J(N-1)}
$$

となる分割を考えて、下限を $\sigma\in\Sigma_{ss'},\nu\in\mathcal{N}_{JJ'}$ に対して

$$
\begin{split}
&&\min_{R,\sigma,\nu}\left(S_{1J}+S_{2J'}-2\text{tr}\left[RC_{s1J}\right]-2\text{tr}\left[RC_{s'2J'}\right]+\sum_{I=3}^N\left(S_{I\nu(I)}-2\text{tr}\left[RC_{\sigma(I)I\nu(I)}\right]\right)\right)\\
&\ge&S_{1J}+S_{2J'}-2\max_{R_2}\text{tr}\left\[R\_2\left(C\_{s1J}+C_{s'2J'}\right)\right\]+\min_{\nu\in\mathcal{N}\_{JJ'}}\sum\_{I=3}^N\left(S_{I\nu(I)}-2\max_{R_I,\sigma}\text{tr}\left\[R\_IC\_{\sigma(I)I\nu(I)}\right\]\right)\\
&=&S_{1J}+S_{2J'}-2\text{tr}\left[R_2^b\left(C_{s1J}+C_{s'2J'}\right)\right]+\sum_{I=3}^NL_{I\nu^b_2(I)}\\
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
\textrm{tr}(P^T\Sigma_bP\Sigma_A)\le\sqrt{\textrm{tr}(P^T\Sigma_b^TPP^T\Sigma_bP)\textrm{tr}(\Sigma_A^T\Sigma_A)}\\
=\textrm{tr}(\Sigma_b)\textrm{tr}(\Sigma_A)\\
=\left(\sum_{i=1}^n\lambda_i^{(A)}\right)\left(\sum_{i=1}^n\lambda_i^{(B)}\right)
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
P_{1} &\mathbf 0\\
\mathbf 0 &P_{2}\\
\end{pmatrix}
$$
とする。このとき...

Rを固定してPの変分下界を計算することが必要？





## 変分上界

最小化RMSDを見つけるのは難しい。そこでまず $\mathrm{tr}[X^TRY\tilde S]$ の変分上界を評価することにする。 $\mathcal O_n$ を $n$ 次元直交行列全体の集合として、$P\in\mathcal O_N,P_I=e^T_Ie_IP,Q_I\in\mathcal O_M$ に対して
$$
S=\sum_{I=1}^NP_I\otimes Q_I\\
$$
と定義する。このとき、$\mathcal S_n\subset\mathcal O_n$から、
$$
\underset{\tilde S\in\mathcal S_N\times(\sigma_M)^N,R\in\mathcal{SO_d}}{\max\arg}\mathrm{tr}[X^TRY\tilde S]\le\underset{S\in\mathcal O_N\times(\mathcal O_M)^N,R\in\mathcal{SO_d}}{\max\arg}\mathrm{tr}[X^TRYS]
$$
が成り立つ。

$P^{(MN)}=P\otimes E_M$ は任意の分子置換、$Q^{(MN)}=\sum_{I=1}^Ne_I^Te_IQ_I$ は任意の分子内回転を表す置換行列を部分集合に含む直交行列である。したがって $P^{(MN)}Q^{(MN)}$は任意の分子対称性を制限された置換行列を部分集合に含む。一般に $P^{(MN)},Q^{(MN)}$ は非可換で$Q^{(MN)}=E_N\otimes Q$ のとき、つまり $Q_1=Q_2=\dots=Q_N$ のときのみ可換になる（$P^{(MN)}Q^{(MN)}=(P\otimes E_M)(E_N\otimes Q)=(E_N\otimes Q)(P\otimes E_M)=Q^{(MN)}P^{(MN)}$）。

トレースを変形する。
$$
\begin{split}
&\textrm{tr}[X^TRYS]\\
=&\textrm{tr}\left[\left(\sum_{I=1}^Ne_I\otimes X_I\right)^T\left(E_1\otimes R\right)\left(\sum_{J=1}^Ne_J\otimes Y_J\right)\left(\sum_{K=1}^NP_K\otimes Q_K\right)\right]\leftarrow 定義より\\
=&\sum_{I=1}^N\sum_{J=1}^N\sum_{K=1}^N\textrm{tr}\left[\left(e_I\otimes X_I\right)^T\left(E_1\otimes R\right)\left(e_J\otimes Y_J\right)\left(P_K\otimes Q_K\right)\right]\leftarrowトレースの線形性\\
=&\sum_{I=1}^N\sum_{J=1}^N\sum_{K=1}^N\textrm{tr}\left[e_I^Te_JP_K\otimes X_I^TRY_JQ_K\right]\leftarrow クロネッカー積の混合積性\\
=&\sum_{I=1}^N\sum_{K=1}^N\textrm{tr}\left[\sum_{J=1}^Np_{JK}\left(e_I^Te_K\otimes C_{IJ}Q_K\right)\right]\leftarrow e_JP_K=p_{JK}e_K\\
=&\sum_{I=1}^N\sum_{K=1}^N\textrm{tr}\left[(e_I^Te_K)\otimes \sum_{J=1}^Np_{JK}(C_{IJ}Q_K)\right]\leftarrow クロネッカー積の線形性\\
=&\sum_{I=1}^N\textrm{tr}\left[\sum_{J=1}^Np_{JI}(C_{IJ}Q_I)\right]\leftarrow\textrm{tr}[(e_I^Te_J)\otimes A]=\delta_{IJ}\textrm{tr}[A]\\
=&\sum_{I=1}^N\textrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]\leftarrow 行列の線形性とトレースの循環性\\
\end{split}
$$
ここで
$$
C_{IJ}=X^T_IRY_J\in\mathbb{R}^{M\times M}
$$
を定義した。$Q_I,C_{IJ}\in\mathbb{R}^{M\times M}$ であり、その行列積は$N$回評価すればよい（$\mathcal O(NM^3)$）。これは$NM\times NM$の行列積（$\mathcal O(N^3M^3)$）を計算するよりははるかに評価が簡単になる。

## 次元削減

$P,Q$を同時に解くのは難しいが、行列
$$
C_{Q_I}(R,P)=\sum_{J=1}^Np_{JI}C_{IJ}
$$
$$
C_P(R,Q)=\left\{\mathrm{tr}\left[Q_I\sum_{J=1}^N C_{IJ}\right]\right\}_{IJ}
$$

に対して直交プロクラステス問題を繰り返し解くことで、$P,Q$を自己無撞着に求めることができる。ここで $d\le n,d\le m$ の場合、$\sum_{J=1}^Np_{JI}C_{IJ}$ や $\{\mathrm{tr}[Q_I\sum_{J=1}^N C_{IJ}]\}_{IJ}$ の持つ次元はたかだか $d$ に制限されるため、それぞれフル次元で求解することは非効率で、解の不定性と不安定性をもたらす。そこで事前処理として、$X,Y$の有効次元を落としておく。
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


そこでまずPCAによる次元削減を行い、$d$ 次元以下の直交プロクラステス問題を解くことを考える。

まずは $C_Q$ の次元削減を考える。
$$
Z_{IJ}=X_{I}-RY_{J}\in\mathbb{R}^{d\times m}
$$
に対して、対称行列
$$
S_{IJ}=Z_{IJ}^TZ_{IJ}
$$
を作る。このとき固有ベクトルとなる直交行列 $W_{IJ}$ を取ることができる。
$$
\left(\sum_{J=1}^NW^T_{IJ}Z_{IJ}^TZ_{IJ}W_{IJ}\right)=\sum_{J=1}^N\Lambda_{IJ}
$$
実効次元が制限されるため、任意の$IJ$と$k>d$ に対して、$\lambda_{IJ,1}=\dots=\lambda_{IJ,M}=0$ が成り立つように取れば、任意の$\lambda_k=0$ に対応する固有ベクトル $w_k$ は
$$
X_Iw=RY_Jw
$$

$$
w_k^T(X_I^TX_I+Y_J^TY_J-X^T_IRY_J-Y^T_JR^TX_I)w_k=0
$$

を満たすから、式変形すると
$$
w_k^TC_{IJ}w_k=w_k^TX_I^TX_Iw_k
$$
となり、直交プロクラステス問題は
$$
\mathrm{tr}[Q_IC_{P_I}]\\
=\mathrm{tr}\left[\sum_{J=1}^Np_{JI}C_{IJ}\right]\\
=\mathrm{tr}[W^TQ_IWW^TC_{IJ}W]\\
=\mathrm{tr}\left[Q_I'W'^T(X_I^TRY_J)W'\right]+\mathrm{tr}\left[Q_I''W''^T(X_I^TRY_J)W''\right]\\
=\mathrm{tr}\left[Q'_I(W'^TX_I^T)(RY_JW')\right]+\mathrm{tr}\left[Q_I''X_I^TX_I\right]\\
$$
とすることができる。ここで右辺第二項を最大化する直交行列は $Q_I''=E_{m-d}$ であると直ちに定まるから、$Q_I'\in\mathbb{R}^{d\times d}$ を求める問題に帰着する。


$$
W'^TC_{IJ}W'
$$


## 目的関数の勾配

ラグランジュの未定乗数 $\lambda=(\lambda_P,\lambda_Q,\lambda_R,\lambda_{|R|})$ に対して
$$
f(P,Q,R,\lambda)=\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}X^T_IRY_J\right]-\lambda_Pg(P)-\lambda_Qg(Q)-\lambda_Rg(R)-\lambda_{|R|}(\det(R)-1)^2
$$
を目的関数として最大化することを考える。ここで $g$ は$A\in\mathbb{R}^{n\times m}$ に対して
$$
g(A)=\|AA^T-E_n\|^2=\mathrm{tr}[(AA^T)^2]-2\mathrm{tr}[AA^T]+n
$$
で定義される。

$P,Q$ の拘束条件より
$$
f(P,Q_1,Q_2,\dots,Q_N,R,\lambda)\\
=\mathrm{tr}\left[\sum_{I=1}^NQ_I\sum_{J=1}^Np_{JI}C_{IJ}\right]-\lambda_Pg(P)-\lambda_Q\sum_{I=1}^Ng(Q_I)-\lambda_Rg(R)-\lambda_{|R|}(\det(R)-1)^2
$$
の最大化問題へと簡略化される。

$P$ について偏微分を行うと
$$
\frac{\partial f}{\partial P}=\frac{\partial}{\partial P}\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]-\lambda_P\frac{\partial g(P)}{\partial P}\\
$$
と計算できる。

第一項は
$$
\frac{\partial}{\partial P}\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]=\left\{\frac{\partial}{\partial p_{mn}}\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]\right\}_{mn}
=\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}\\
$$
と求まる。第二項は
$$
\lambda_P\frac{\partial g(P)}{\partial P}=4\lambda_P((\mathrm{Tr}[(P^TP)^2])'-P)=4\lambda_P(PP^T-I)P
$$
と計算でき、ここで、$\mathbb{R}^{n\times n}\ni A=\{a_{ij}\}_{ij}$ に対して
$$
(\mathrm{Tr}[(A^TA)^2])'=\left\{a_{ij}^3+\left(\sum_{k\ne j}a_{ik}^2+\sum_{k\ne i}a_{kj}^2\right)a_{ij}+\sum_{k\ne i,l\ne j}a_{il}a_{kj}a_{kl}\right\}_{ij}\\
=\left\{\sum_{k=1}^N\sum_{l=1}^Na_{il}a_{kj}a_{kl}\right\}_{ij}\\
=AA^TA
$$
まとめると
$$
\frac{\partial f}{\partial P}=\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}-4\lambda_P(PP^T-I)P\\
$$
が得られる。

$Q_I$についても同様に
$$
\frac{\partial}{\partial Q_I}\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]=\sum_{J=1}^Np_{JI}C_{IJ}^T\\
$$
であるから
$$
\frac{\partial f}{\partial Q_I}=\sum_{J=1}^Np_{JI}C_{IJ}^T-4\lambda_Q(Q_IQ_I^T-E)Q_I\\
$$
と定まる。

$R$については
$$
\sum_{I=1}^N\textrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]=\sum_{I=1}^N\sum_{J=1}^Np_{JI}\textrm{tr}\left[Q_IX^T_IRY_J\right]\\
=\sum_{I=1}^N\sum_{J=1}^Np_{JI}\textrm{tr}\left[Y_JQ_IX^T_IR\right]\\
$$
と変形しておいて
$$
\frac{\partial}{\partial R}\sum_{I=1}^N\sum_{J=1}^Np_{JI}\textrm{tr}\left[Y_JQ_IX^T_IR\right] =\sum_{I=1}^N\sum_{J=1}^Np_{JI}\left(Y_JQ_IX^T_I\right)^T\\
=\sum_{I=1}^N\left(X_IQ_I^T\left(\sum_{J=1}^Np_{JI}Y_J^T\right)\right)
$$
と計算できる。この形式は部分回転を考える際に計算不要な部分をスキップできるためクロネッカー積で表される$\mathbb{R}^{MN\times MN}$ の直行行列を用いるより有利である。まとめて
$$
\frac{\partial f}{\partial R}=\sum_{I=1}^N\left(X_IQ_I^T\left(\sum_{J=1}^Np_{JI}Y_J^T\right)\right)-4\lambda_R(RR^T-E)R-2\lambda_{|R|}(\det(R)-1)\mathrm{adj}(R)\\
$$
ここで$\mathrm{adj}(R)$ は $R$ の余因子行列を表す。

未定乗数の勾配についても同様に
$$
\frac{\partial f}{\partial \lambda_P}=g(P)\\
$$
$$
\frac{\partial f}{\partial \lambda_Q}=\sum_{I=1}^Ng(Q_I)\\
$$
$$
\frac{\partial f}{\partial \lambda_R}=(R)\\
$$
$$
\frac{\partial f}{\partial \lambda_{|R|}}=(\det(R)-1)^2\\
$$

と計算される。

$P,Q$ を固定した状態での $R$ の推定は Kabsch-Umeyama あるいは四元数を用いた陽な解法がある。どっち使ったら良いかはわからない。

