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



$M$分子からなる同種分子が $N$ 個からなる分子クラスターの最小化RMSDの変分下界を求める問題は、分子構造を $X,Y\in\R^{d\times MN}$ で表すことにすれば、
$$
(P^*,Q^*,R^*)=\underset{P^{(MN)},Q^{(MN)},R}\max\arg\mathrm{tr}[X^TRYP'Q']
$$
を求める問題に帰着される。$X,Y$ は $M$ が先に回る形で格納されるとする。
$$
X=\{x_{11},\dots,x_{M1},x_{12},\dots,x_{M2},x_{1N},\dots,x_{MN}\},&x_{IJ}\in \R^{d}\\
$$
ここで $P',Q'\in \R^{MN\times MN},R\in\R^{d\times d}$ で、$P'P'=Q'Q'^T=E_{MN},RR^T=E_d,\det(R)=1$ を満たすとする。ただし、$E_X$ は$\R^{X\times X}$ の単位行列。ここで、問題に合わせて不要な自由度を削除すると、
$$
P'&=&P\otimes E_M,&P\in\R^{N\times N}\\
Q'&=&\sum_{I=1}^N e_I\otimes Q_I,&e_I\in\R^{N\times N},Q_I\in\R^{M\times M},1\le I\le N
$$
ここでクロネッカーのデルタ $\delta_{ij}$ に対して $e_I=\{\delta_{iI}\delta_{Ij}\}_{ij}$ である。定義から$PP^T=E_N,Q_IQ_I^T=E_M$である（多分）。したがって
$$
P'Q'=\sum_{I=1}^N P_Ne_I\otimes Q_I=\sum_{I=1}^NP_{I}\otimes Q_I
$$
と表される。定義より $P=(p_1,p_2,\dots,p_I,\dots,p_N)$ に対して $P_I=(0,0,\dots,p_I,\dots,0)$ である。

$P$は任意の分子置換、$Q$は任意の分子内回転を表す置換行列を部分集合に含み、直交行列$PQ$は任意の分子対称性を制限された置換行列を部分集合に含む。一般に $P,Q$ は非可換で$Q=Q_M\otimes E_N$ つまり $Q_1=Q_2=\dots=Q_N$ のときのみ可換になる（$P'=P\otimes E,Q'=E\otimes Q$）。

目的のトレースを変形する。
$$
\mathrm{tr}[X^TRYP'Q']&=&\mathrm{tr}\left[X^TRY\left(\sum_{I=1}^NP_I\otimes Q_I\right)\right]\\
&=&\sum_{I=1}^N\mathrm{tr}\left[X^TRY\left(P_I\otimes Q_I\right)\right]\\
$$
ここで、部分構造（分子座標） $X_I,Y_I\in\R^{d\times M}$ を$X_I=\{x_{1I},\dots,x_{MI}\},Y_I=\{y_{1I},\dots,y_{MI}\}$ で定義し、
$$
C_{IJ}=X^T_IRY_J\in\R^{M\times M}
$$
を定義すれば、行列の行列 $C=\{C_{kl}\}_{kl} \in\left(\R^{M\times M}\right)^{N\times N}$ は $X^TRY$ を区分行列表現したものになる。

$P_I=(0,0,\dots,p_I,\dots,0)$であるから、
$$
P_I\otimes Q_I=\{\delta_{Il}p_{kI}Q_I\}_{kl} \in\left(\R^{M\times M}\right)^{N\times N}
$$
で表記できるから
$$
\left(Q_I\otimes P_I\right)C=\left\{\sum_{j=1}^N\delta_{Ij}p_{kI}Q_IC_{jl}\right\}_{kl}=\left\{p_{kI}Q_IC_{Il}\right\}_{kl}
$$
と変形できる。トレースを取れば
$$
\mathrm{tr}\left[\left(Q_I\otimes P_I\right)C\right]&=&\mathrm{tr}\left[Q_I\sum_{J=1}^Np_{JI}C_{IJ}\right]\\
$$
最終的に求めるトレースは
$$
\mathrm{tr}\left[P'Q'C\right]&=&\sum_{I=1}^N\mathrm{tr}\left[\left(Q_I\otimes P_I\right)C\right]&=&\mathrm{tr}\left[\sum_{I=1}^NQ_I\sum_{J=1}^Np_{JI}C_{IJ}\right]\\
$$
に帰着される。$Q_I,C_{IJ}$ は $M\times M$ 行列であり、その行列積は$N$回評価すればよい（$\mathcal O(NM^3)$）。これは$NM\times NM$の行列積（$\mathcal O(N^3M^3)$）を計算するよりははるかに評価が簡単になる。

## 目的関数の勾配

ラグランジュの未定乗数 $\lambda=(\lambda_P,\lambda_Q,\lambda_R,\lambda_{|R|})$ に対して
$$
f(P,Q,R,\lambda)=\mathrm{tr}[XPQY^TR]-\lambda_Pg(P)-\lambda_Qg(Q)-\lambda_Rg(R)-\lambda_{|R|}(\det(R)-1)^2
$$
を目的関数として最大化することを考える。ここで $g$ は$A\in\R^{n\times n}$ に対して
$$
g(A)=\|AA^T-I_N\|^2=\mathrm{tr}[(AA^T)^2]-2\mathrm{tr}[AA^T]+n
$$
で定義される。

$P,Q$ の拘束条件よりl
$$
f(P,Q_1,Q_2,\dots,Q_N,R,\lambda)\\
=\mathrm{tr}\left[\sum_{I=1}^NQ_I\sum_{J=1}^Np_{JI}C_{IJ}\right]-\lambda_Pg(P)-\lambda_Q\sum_{I=1}^Ng(Q_I)-\lambda_Rg(R)-\lambda_{|R|}(\det(R)-1)^2
$$
の最大化問題へと簡略化される。

$P$ について偏微分を行うと
$$
\frac{\partial f}{\partial P}&=&\sum_{I=1}^N\frac{\partial\mathrm{tr}\left[\left(P_I\otimes Q_I\right)C\right]}{\partial P}-\lambda_P\frac{\partial g(P)}{\partial P}\\
$$
と計算できる。第一項は
$$
\sum_{I=1}^N\frac{\partial\mathrm{tr}\left[\left(P_I\otimes Q_I\right)Y^TRX\right]}{\partial P}&=&\left\{\frac{\partial}{\partial p_{mn}}\sum_{I=1}^N\sum_{j=1}^Np_{jI}\ \mathrm{tr}\left[Q_IC_{Ij}\right]\right\}_{mn}\\
&=&\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}
$$
と求まる。第二項の偏微分は
$$
\frac{\partial g(P)}{\partial P}=4(h(P)-P)=4(PP^T-I)P
$$
と計算でき、ここで、$\R^{n\times n}\ni A=\{a_{ij}\}_{ij}$ に対して
$$
h(A)&=&\left\{a_{ij}^3+\left(\sum_{k\ne j}a_{ik}^2+\sum_{k\ne i}a_{kj}^2\right)a_{ij}+\sum_{k\ne i,l\ne j}a_{il}a_{kj}a_{kl}\right\}_{ij}\\
&=&\left\{\sum_{k=1}^N\sum_{l=1}^Na_{il}a_{kj}a_{kl}\right\}_{ij}\\
&=&AA^TA
$$
と定義した（これは$\mathrm{Tr}[(A^TA)^2]$ の偏微分/4）。まとめると
$$
\left\{\frac{\partial f}{\partial P_N}\right\}_{mn}&=&\mathrm{tr}\left[Q_nC_{nm}\right]-4\lambda_P\left\{(PP^T-I)P\right\}_{mn}\\
$$
が得られる。

同様に
$$
\frac{\partial}{\partial Q_I}\sum_{J=1}^N\mathrm{tr}\left[\left(P_J\otimes Q_J\right)Y^TRX\right]&=&\frac{\partial}{\partial Q_I}\sum_{J=1}^N\sum_{K=1}^Np_{JK}\ \mathrm{tr}\left[Q_KC_{KJ}\right]
&=&\sum_{J=1}^Np_{JI}C_{IJ}^T\\
$$
であるから
$$
\frac{\partial f}{\partial Q_I}&=&\sum_{J=1}^Np_{JI}C_{IJ}^T-4\lambda_Q(Q_IQ_I^T-E)Q_I\\
$$
と定まる。

$R$については
$$
\frac{\partial f}{\partial R}&=&YP'Q'X^T-4\lambda_R(RR^T-E)R-2\lambda_{|R|}(\det(R)-1)\mathrm{adj}(R)\\
$$
ここで$\mathrm{adj}(R)$ は $R$ の余因子行列を表す。第一項は$F_I=(\underset{1}0,\underset{2}0,\dots,\underset{I}1,\dots,\underset{N}0)\in\R^{1\times N}$ を用いて
$$
YP'Q'X^T&=&Y\left(\sum_{I=1}^NP_I\otimes Q_I\right)X^T\\
&=&\sum_{I=1}^N\left(\sum_{J=1}^NF_J\otimes Y_J\right)(P_I\otimes Q_I)\left(\sum_{K=1}^NF_K^T\otimes X_K^T\right)\\
&=&\sum_{I=1}^N\sum_{J=1}^N\sum_{K=1}^N\left(F_JP_IF_K^T\otimes Y_JQ_IX_K^T\right)\\
&=&\sum_{I=1}^N\sum_{J=1}^N\left(F_JP_IF_J^T\otimes Y_JQ_IX_J^T\right)\\
&=&\sum_{I=1}^N\sum_{J=1}^N\left(p_{JI} Y_JQ_IX_J^T\right)\\
$$
と評価できる（が、直接評価するのとどっちがいいかはわからない）。

未定乗数についても同様に
$$
\frac{\partial f}{\partial \lambda_P}&=&g(P)\\
\frac{\partial f}{\partial \lambda_Q}&=&\sum_{I=1}^Ng(Q_I)\\
\frac{\partial f}{\partial \lambda_R}&=&g(R)\\
\frac{\partial f}{\partial \lambda_{|R|}}&=&(\det(R)-1)^2\\
$$
と計算される。

ただし、$P,Q$ を固定した状態での $R$ の推定は Kabsch-Umeyama あるいは四元数を用いた陽な解法がある。どっち使ったら良いかはわからない。
