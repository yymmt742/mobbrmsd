# symRMSD
molecular symmetry corrected RMSD by branch and prune

installation
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=[Release, Debug, ...] -DCMAKE_PREFIX_PATH=your_local_cmake_path
make install

RMSD最小化のメモ


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



## 分子座標表現

$M$分子からなる同種分子が $N$ 個からなる分子クラスターの分子構造を、原子座標$x_{ij}\in\mathbb{R}^d$ のベクトル $X\in\mathbb{R}^{d\times MN}$ で表すことにする。$X$ の内部では $M$ が先に回る形で格納される。
$$
X=\{x_{11},\dots,x_{M1},x_{12},\dots,x_{M2},x_{1N},\dots,x_{MN}\},x_{ij}\in \mathbb{R}^{d}\\
$$
指示ベクトル $e_I$ を
$$
e_I=(\underset{I-1}{0,0,\dots},1,\underset{N-I}{\dots,0})\in\mathbb{R}^{1\times N}
$$
で定義する。指示ベクトルに対して
$$
e_Ie_J^T=\delta_{IJ}\in\mathbb{R}^{1\times1}\\
e_I^Te_J=\{\delta_{iJ}\delta_{Ij}\}_{ij}\in\mathbb{R}^{N\times N}
$$
などが成り立つ。ここで $\delta_{ij}$ はクロネッカーのデルタである。

個別の分子座標を
$$
X_I=\{x_{1I},x_{2I},\dots,x_{MI}\}\in\mathbb{R}^{d\times M}
$$
で表すことにすれば
$$
X=\sum_{I=1}^Ne_I\otimes X_I
$$
と表すことができる。ここで $\otimes$ はクロネッカー積である。

## 分子置換

$n$ 次元置換行列全体の集合を $\mathcal S_n\subset\mathbb{Z}^{n\times n}$ と書く。分子集合体構造 $X\in\mathbb{R}_{d\times MN}$ に対して、分子置換は$P'\in\mathcal S_N$ を用いて
$$
\tilde P=P'\otimes E_M
$$
と表現できる。また、分子内の置換（対称性）は任意の $\sigma_M\subset\mathcal S_M$ に対して$Q'_I\in\sigma_M,I=1,\dots,N$ を用いて
$$
\tilde Q=\sum_{I=1}^Ne_I^Te_I\otimes Q'_I
$$
 で表現できる。任意の分子置換操作は
$$
\tilde S=\tilde P\tilde Q\\
=\left(P'\otimes E_M\right)\left(\sum_{I=1}^Ne_I^Te_I\otimes Q'_I\right)\\
=\sum_{I=1}^NPe_I^Te_I\otimes Q'_I\\
=\sum_{I=1}^NP'_I\otimes Q'_I
$$
で表現できる。ここで $P'=\begin{pmatrix}p'_1&p'_2&\cdots&p_i'&\cdots&p'_N\end{pmatrix},p_i\in\mathbb{Z}^{N\times 1}$ と表記することにすれば
$$
P'_I=\begin{pmatrix}0&0&\cdots&p_I'&\cdots&0\end{pmatrix}
$$
である。上記の議論より、$\tilde S\in S_N\times(\sigma_M)^N$ である。

分子置換のみを許した最小化RMSDを求める問題は、
$$
(\tilde S^*,R^*)=\underset{\tilde S\in\mathcal S_N\times(\sigma_M)^N,R\in\mathcal{SO_d}}{\max\arg}\mathrm{tr}[X^TRY\tilde S]
$$
に帰着される。ここで$\mathcal{SO}_d$ は $d$ 次元回転行列全体の集合を表す。

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
