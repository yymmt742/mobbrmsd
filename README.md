# symRMSD
molecular symmetry corrected RMSD by branch and prune


RMSD最小化のメモ


$A,B\in\R^{d\times n},R\in\R^{^{d\times d}},P\in\{0,1\}^{n\times n}$ とする。$R^TR=I,P^TP=I$である。一般化RMSDは $\mathrm{tr}(RBPA^T)$の最大化問題に帰着される。

コーシーシュワルツの不等式の一般化により
$$
\mathrm{tr}(RBPA^T)^2&\le&\mathrm{tr}(R^TR)\mathrm{tr}((BPA^T)^TBPA^T)\\
&=&d\cdot\mathrm{tr}\left(S^TS\right)\\
$$
がなりたつ。ここで$S=BPA^T$とおいた。

※これ以上分解すると
$$
\mathrm{tr}(P^T\Sigma_bP\Sigma_A)&\le&\sqrt{\mathrm{tr}(P^T\Sigma_b^TPP^T\Sigma_bP)\mathrm{tr}(\Sigma_A^T\Sigma_A)}\\
&=&\mathrm{tr}(\Sigma_b)\mathrm{tr}(\Sigma_A)\\
&=&\left(\sum_{i=1}^n\lambda_i^{(A)}\right)\left(\sum_{i=1}^n\lambda_i^{(B)}\right)
$$
となり自明な情報しか得られない。

ここで正規直交行列$\tilde P\in\R^{n\times n}$が$\tilde P^T\tilde P=I$を満たすとする。このとき、Kabschのアルゴリズムを用いて、特異値分解
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

$M$分子からなる同種分子が $N$ 個からなる分子クラスターの分子構造を、原子座標$x_{ij}\in\R^d$ のベクトル $X\in\R^{d\times MN}$ で表すことにする。$X$ の内部では $M$ が先に回る形で格納される。
$$
X=\{x_{11},\dots,x_{M1},x_{12},\dots,x_{M2},x_{1N},\dots,x_{MN}\},&x_{ij}\in \R^{d}\\
$$
指示ベクトル $e_I$ を
$$
e_I=(\underset{I-1}{0,0,\dots},1,\underset{N-I}{\dots,0})\in\R^{1\times N}
$$
で定義する。指示ベクトルに対して
$$
e_Ie_J^T=\delta_{IJ}\in\R^{1\times1}\\
e_I^Te_J=\{\delta_{iJ}\delta_{Ij}\}_{ij}\in\R^{N\times N}
$$
などが成り立つ。ここで $\delta_{ij}$ はクロネッカーのデルタである。

個別の分子座標を
$$
X_I=\{x_{1I},x_{2I},\dots,x_{MI}\}\in\R^{d\times M}
$$
で表すことにすれば
$$
X=\sum_{I=1}^Ne_I\otimes X_I
$$
と表すことができる。ここで $\otimes$ はクロネッカー積である。

## 分子置換

$n$ 次元置換行列全体の集合を $\mathcal S_n\subset\Z^{n\times n}$ と書く。分子集合体構造 $X\in\R_{d\times MN}$ に対して、分子置換は$P'\in\mathcal S_N$ を用いて
$$
\tilde P=P'\otimes E_M
$$
と表現できる。また、分子内の置換（対称性）は任意の $\sigma_M\sub\mathcal S_M$ に対して$Q'_I\in\sigma_M,I=1,\dots,N$ を用いて
$$
\tilde Q=\sum_{I=1}^Ne_I^Te_I\otimes Q'_I
$$
 で表現できる。任意の分子置換操作は
$$
\tilde S&=&\tilde P\tilde Q\\
&=&\left(P'\otimes E_M\right)\left(\sum_{I=1}^Ne_I^Te_I\otimes Q'_I\right)\\
&=&\sum_{I=1}^NPe_I^Te_I\otimes Q'_I\\
&=&\sum_{I=1}^NP'_I\otimes Q'_I
$$
で表現できる。ここで $P'=\begin{pmatrix}p'_1&p'_2&\cdots&p_i'&\cdots&p'_N\end{pmatrix},p_i\in\Z^{N\times 1}$ と表記することにすれば
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
S&=&\sum_{I=1}^NP_I\otimes Q_I\\
$$
と定義する。このとき、$\mathcal S_n\sub\mathcal O_n$から、
$$
\underset{\tilde S\in\mathcal S_N\times(\sigma_M)^N,R\in\mathcal{SO_d}}{\max\arg}\mathrm{tr}[X^TRY\tilde S]\le\underset{S\in\mathcal O_N\times(\mathcal O_M)^N,R\in\mathcal{SO_d}}{\max\arg}\mathrm{tr}[X^TRYS]
$$
が成り立つ。

$P^{(MN)}=P\otimes E_M$ は任意の分子置換、$Q^{(MN)}=\sum_{I=1}^Ne_I^Te_IQ_I$ は任意の分子内回転を表す置換行列を部分集合に含む直交行列である。したがって $P^{(MN)}Q^{(MN)}$は任意の分子対称性を制限された置換行列を部分集合に含む。一般に $P^{(MN)},Q^{(MN)}$ は非可換で$Q^{(MN)}=E_N\otimes Q$ のとき、つまり $Q_1=Q_2=\dots=Q_N$ のときのみ可換になる（$P^{(MN)}Q^{(MN)}=(P\otimes E_M)(E_N\otimes Q)=(E_N\otimes Q)(P\otimes E_M)=Q^{(MN)}P^{(MN)}$）。

トレースを変形する。
$$
&&\mathrm{tr}[X^TRYS]\\
&=&\mathrm{tr}\left[\left(\sum_{I=1}^Ne_I\otimes X_I\right)^T\left(E_1\otimes R\right)\left(\sum_{J=1}^Ne_J\otimes Y_J\right)\left(\sum_{K=1}^NP_K\otimes Q_K\right)\right]&\leftarrow&定義より\\
&=&\sum_{I=1}^N\sum_{J=1}^N\sum_{K=1}^N\mathrm{tr}\left[\left(e_I\otimes X_I\right)^T\left(E_1\otimes R\right)\left(e_J\otimes Y_J\right)\left(P_K\otimes Q_K\right)\right]&\leftarrow&トレースの線形性\\
&=&\sum_{I=1}^N\sum_{J=1}^N\sum_{K=1}^N\mathrm{tr}\left[e_I^Te_JP_K\otimes X_I^TRY_JQ_K\right]&\leftarrow&クロネッカー積の混合積性\\
&=&\sum_{I=1}^N\sum_{K=1}^N\mathrm{tr}\left[\sum_{J=1}^Np_{JK}\left(e_I^Te_K\otimes C_{IJ}Q_K\right)\right]&\leftarrow&e_JP_K=p_{JK}e_K\\
&=&\sum_{I=1}^N\sum_{K=1}^N\mathrm{tr}\left[(e_I^Te_K)\otimes \sum_{J=1}^Np_{JK}(C_{IJ}Q_K)\right]&\leftarrow&クロネッカー積の線形性\\
&=&\sum_{I=1}^N\mathrm{tr}\left[\sum_{J=1}^Np_{JI}(C_{IJ}Q_I)\right]&\leftarrow&\mathrm{tr}[(e_I^Te_J)\otimes A]=\delta_{IJ}\mathrm{tr}[A]\\
&=&\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]&\leftarrow&行列の線形性とトレースの循環性\\
$$
ここで
$$
C_{IJ}=X^T_IRY_J\in\R^{M\times M}
$$
を定義した。$Q_I,C_{IJ}\in\R^{M\times M}$ であり、その行列積は$N$回評価すればよい（$\mathcal O(NM^3)$）。これは$NM\times NM$の行列積（$\mathcal O(N^3M^3)$）を計算するよりははるかに評価が簡単になる。

## 目的関数の勾配

ラグランジュの未定乗数 $\lambda=(\Lambda_P,\Lambda_{Q_1},\dots,\Lambda_{Q_N},\Lambda_R,\lambda_{|R|})$ に対して
$$
f(P,Q_1,Q_2,\dots,Q_N,R,\lambda)=\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}X^T_IRY_J\right]-g(P,\Lambda_P)-\sum_{I=1}^{N}g(Q_I,\Lambda_{Q_I})-g(R,\Lambda_R)-\lambda_{|R|}(\det(R)-1)
$$
を目的関数として最大化することを考える。ここで $g$ は$\R^{m\times n}\ni A=(a_1,\dots,a_n)$ に対して
$$
g(A,\Lambda)=\sum_{i=1}^n\sum_{j=1}^n \lambda_{ij}(a_i^Ta_j-\delta_{ij})
$$
で定義される。

$P$について偏微分を行う。
$$
\frac{\partial f}{\partial P}&=&\frac{\partial}{\partial P}\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]-\frac{\partial g(P,\Lambda_P)}{\partial P}\\
$$
第一項は
$$
\frac{\partial}{\partial P}\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]&=&\left\{\frac{\partial}{\partial p_{mn}}\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]\right\}_{mn}
&=&\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}\\
$$
と求まる。第二項は
$$
\frac{\partial g(P,\Lambda_P)}{\partial P}&=&\left\{\sum_{i=1}^n\sum_{j=1}^n\sum_{k=1}^n \lambda_{ij}\frac{\partial}{\partial p_{mn}}(p_{ki}p_{kj}-\delta_{ij})\right\}_{mn}\\
&=&\left\{\sum_{i=1}^n\lambda_{ni}p_{mi}+\sum_{i=1}^n\lambda_{in}p_{mi}\right\}_{mn}\\
&=&(\Lambda_P+\Lambda_P^T)P
$$
と計算できる。

$Q_I$についても同様に
$$
\frac{\partial}{\partial Q_I}\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]&=&\sum_{J=1}^Np_{JI}C_{IJ}^T\\
$$
である。

$R$については
$$
\sum_{I=1}^N\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]&=&\sum_{I=1}^N\sum_{J=1}^Np_{JI}\mathrm{tr}\left[Q_IX^T_IRY_J\right]\\
&=&\sum_{I=1}^N\sum_{J=1}^Np_{JI}\mathrm{tr}\left[Y_JQ_IX^T_IR\right]\\
$$
と変形しておいて
$$
\frac{\part}{\part R}\sum_{I=1}^N\sum_{J=1}^Np_{JI}\mathrm{tr}\left[Y_JQ_IX^T_IR\right] &=&\sum_{I=1}^N\sum_{J=1}^Np_{JI}\left(Y_JQ_IX^T_I\right)^T\\
&=&\sum_{I=1}^N\left(X_IQ_I^T\left(\sum_{J=1}^Np_{JI}Y_J^T\right)\right)
$$
と計算できる。この形式は部分回転を考える際に計算不要な部分をスキップできるためクロネッカー積で表される$\R^{MN\times MN}$ の直行行列を用いるより有利である。まとめて
$$
\frac{\partial f}{\partial R}&=&\sum_{I=1}^N\left(X_IQ_I^T\left(\sum_{J=1}^Np_{JI}Y_J^T\right)\right)-(\Lambda_{R}+\Lambda_{R}^T)R-\lambda_{|R|}\mathrm{adj}(R)\\
$$
ここで$\mathrm{adj}(R)$ は $R$ の余因子行列を表す。

未定乗数の勾配についても同様に
$$
\frac{\partial f}{\partial \Lambda_P}&=&P^TP-E_N\\
\frac{\partial f}{\partial \Lambda_{Q_I}}&=&Q_I^TQ_I-E_M\\
\frac{\partial f}{\partial \Lambda_R}&=&R^TR-E_d\\
\frac{\partial f}{\partial \lambda_{|R|}}&=&\det(R)-1\\
$$
と計算される。

まとめれば
$$
\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}&=&(\Lambda_{P}+\Lambda_{P}^T)P\\
\sum_{J=1}^Np_{JI}C_{IJ}^T&=&(\Lambda_{Q_I}+\Lambda_{Q_I}^T)Q_I\\
\sum_{I=1}^N\left(X_IQ_I^T\left(\sum_{J=1}^Np_{JI}Y_J^T\right)\right)&=&(\Lambda_R+\Lambda_R^T)R+\lambda_{|R|}\mathrm{adj}(R)\\
P^TP&=&E_N\\
Q_I^TQ_I&=&E_M\\
R^TR&=&E_d\\
\det(R)&=&1\\
$$
連立して
$$
\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}P^T&=&\Lambda_{P}+\Lambda_{P}^T\\
\left(\sum_{J=1}^Np_{JI}C_{IJ}^T\right)Q_I^T&=&\Lambda_{Q_I}+\Lambda_{Q_I}^T\\
\left(\sum_{I=1}^N\left(X_IQ_I^T\left(\sum_{J=1}^Np_{JI}Y_J^T\right)\right)\right)R^T&=&\Lambda_R+\Lambda_R^T+\lambda_{|R|}\det(R)E_d\\
$$
つまり、左辺は対称行列とならなければらない。ここで特異値分解
$$
\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}=U\Sigma V^T
$$
を考えて、ある$SS^T=E_M$ となる対角行列$S$に対して
$$
P&=&U_PSV^T_P\\
$$
とすれば
$$
PP^T=UV^T(UV^T)^T=UV^TVU^T=E_N
$$

$$
\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}P^T&=&U\Sigma V^T(USV^T)^T\\
&=&U\Sigma V^TVS^TU^T\\
&=&U\Sigma SU^T
$$

であるから、$(U\Sigma SU^T)^T=US^T\Sigma^TU^T=U\Sigma SU^T$ となり、これは未定乗数法の解となっている。特に$S=E_M$としたものは直交プロクラステス問題の解と一致するため、トレースの最大化が保証される。同様に
$$
Q_I&=&U_{Q_I}V_{Q_I}^T\\
R&=&U_R\tilde S_RV_R^T\\
$$
はそれぞれの方程式の解となっている。ただし$\tilde S_R$はKabsch-Umeyama アルゴリズムの $\tilde S$ である。

それぞれの方程式は
$$
P'(\{Q_I\},R)&=&U_PV^T_P\\
Q'_I(P,R)&=&U_{Q_I}V^T_{Q_I}\\
R(P,\{Q_I\})&=&U_R\tilde SV^T_R\\
$$
となり解けない。そこでSCFで最適化することになる。特に $P',Q'$ は分散共分散行列行列を一度計算しておけば繰り返し使い回せるため、$R\rightarrow \{P,\{Q_I\}\}\rightarrow R\rightarrow\cdots$ という順番で計算するのが効率的である。また、$d=3$ の場合、四元数を用いた$R$の高速な推定法があるため、それで代用することができる。

## 一般化





ここで、$\R^{n\times n}\ni A=\{a_{ij}\}_{ij}$ に対して
$$
(\mathrm{Tr}[(A^TA)^2])'&=&\left\{a_{ij}^3+\left(\sum_{k\ne j}a_{ik}^2+\sum_{k\ne i}a_{kj}^2\right)a_{ij}+\sum_{k\ne i,l\ne j}a_{il}a_{kj}a_{kl}\right\}_{ij}\\
&=&\left\{\sum_{k=1}^N\sum_{l=1}^Na_{il}a_{kj}a_{kl}\right\}_{ij}\\
&=&AA^TA
$$
