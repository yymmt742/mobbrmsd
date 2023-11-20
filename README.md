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
(P^*,Q^*,R^*)=\underset{P^{(MN)},Q^{(MN)},R}\max\arg\mathrm{tr}[XP^{(NM)}Q^{(NM)}Y^TR]
$$
を求める問題に帰着される。ただし $X,Y$ は $M$ が先に回る形で格納される。
$$
X=\{x_{11},\dots,x_{M1},x_{12},\dots,x_{M2},x_{1N},\dots,x_{MN}\},&x_{IJ}\in \R^{d}\\
$$
ここで $P^{(MN)},Q^{(MN)}\in \R^{MN\times MN},R\in\R^{d\times d}$ で、$P^{(MN)}P^{(MN)T}=Q^{(MN)}Q^{(MN)T}=I_N,RR^T=I_d,\det(R)=1$ を満たすとする。ただし
$$
P^{(MN)}&=&P\otimes I_M,&I_M\in\R^{M\times M},P\in\R^{N\times N}\\
Q^{(MN)}&=&\sum_{I=1}^N e_I\otimes Q_I,&Q_I\in\R^{M\times M},I=1,\dots,N,e_I\in\R^{N\times N}
$$
であり、クロネッカーのデルタ $\delta_{ij}$ に対して $e_I=\{\delta_{iI}\delta_{Ij}\}_{ij}$ である。定義から$PP^T=I_N,Q_IQ_I^T=I_M$である（多分）。したがって
$$
PQ=\sum_{I=1}^N P_Ne_I\otimes Q_I=\sum_{I=1}^NP_{I}\otimes Q_I
$$
と表される。ただし$P=(p_1,p_2,\dots,p_I,\dots,p_N)$ に対して $P_I=(0,0,\dots,p_I,\dots,0)$ である。

$P$は任意の分子置換、$Q$は任意の分子内回転を表す置換行列を部分集合に含み、直交行列$PQ$は任意の分子対称性を制限された置換行列を部分集合に含む。一般に $P,Q$ は非可換で$Q=Q_M\otimes I_N$ つまり $Q_1=Q_2=\dots=Q_N$ のときのみ可換。

目的のトレースを変形する。
$$
\mathrm{tr}[XPQY^TR]&=&\mathrm{tr}\left[X\left(\sum_{I=1}^NP_I\otimes Q_I\right)Y^TR\right]\\
&=&\sum_{I=1}^N\mathrm{tr}\left[\left(P_I\otimes Q_I\right)Y^TRX\right]\\
$$
ここで $\R^{n\times n}\ni C=Y^TRX$ とおき、部分行列 $C_{IJ}\in\R^{M\times M}$ によって行列の行列 $C=\{C_{kl}\}_{kl} \in\left(\R^{M\times M}\right)^{N\times N}$ と表記することにすれば、$P_I=(0,0,\dots,p_I,\dots,0)$であるから、
$$
P_I\otimes Q_I=\{\delta_{kl}p_{kI}Q_I\}_{kl} \in\left(\R^{M\times M}\right)^{N\times N}
$$
で表記できるから
$$
\left(Q_I\otimes P_I\right)C=\left\{p_{kI}Q_IC_{Il}\right\}_{kl}
$$
と変形できる。これより
$$
\mathrm{tr}\left[\left(Q_I\otimes P_I\right)C\right]&=&\mathrm{tr}\left[p_{kI}Q_IC_{Ik}\right]\\
&=&\mathrm{tr}\left[Q_I\sum_{k=1}^Np_{kI}C_{Ik}\right]
$$
求めるトレース全体は
$$
\sum_{I=1}^N\mathrm{tr}\left[\left(P_I\otimes Q_I\right)C\right]&=&\sum_{I=1}^N\sum_{j=1}^Np_{jI}\cdot \mathrm{tr}\left[Q_IC_{Ij}\right]\\
$$
に帰着される。$Q_I,C_{Ij}$ は $M\times M$ 行列であるから、$NM\times NM$の行列積を計算するよりははるかに評価が簡単である。

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

$P,Q$ の拘束条件より
$$
f(P,Q_1,Q_2,\dots,Q_N,R,\lambda)\\
=\sum_{I=1}^N\mathrm{tr}\left[\left(P_I\otimes Q_I\right)Y^TRX\right]-\lambda_Pg(P)-\lambda_Q\sum_{I=1}^Ng(Q_I)-\lambda_Rg(R)-\lambda_{|R|}(\det(R)-1)^2
$$
の最大化問題へと簡略化される。

$P$ について偏微分を行うと
$$
\frac{\partial f}{\partial P}&=&\sum_{I=1}^N\frac{\partial\mathrm{tr}\left[\left(P_I\otimes Q_I\right)Y^TRX\right]}{\partial P}-\lambda_P\frac{\partial g(P)}{\partial P}\\
$$
と計算できる。第一項は
$$
\sum_{I=1}^N\frac{\partial\mathrm{tr}\left[\left(P_I\otimes Q_I\right)Y^TRX\right]}{\partial P}&=&\left\{\frac{\partial}{\partial p_{mn}}\sum_{I=1}^N\sum_{j=1}^Np_{jI}\ \mathrm{tr}\left[Q_IC_{Ij}\right]\right\}_{mn}\\
&=&\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}
$$
と求まる。第二項の偏微分は
$$
\frac{\partial g(P)}{\partial P}=4(h(P)-P)
$$
と計算でき、こで、$\R^{n\times n}\ni A=\{a_{ij}\}_{ij}$ に対して
$$
h(A)=\left\{a_{ij}^3+\left(\sum_{k\ne i}a_{ik}^2+\sum_{k\ne j}a_{kj}^2\right)a_{ij}+\sum_{k\ne i}\sum_{l\ne j}a_{il}a_{kj}a_{kl}\right\}_{ij}
$$
と定義した（これは$\mathrm{Tr}[(A^TA)^2]$ の偏微分）。まとめると
$$
\left\{\frac{\partial f}{\partial P_N}\right\}_{mn}&=&\left\{\mathrm{tr}\left[Q_nC_{nm}\right]\right\}_{mn}-4\lambda_P(\{h(P)\}_{mn}-P_{mn})\\
$$
が得られる。



同様に
$$
\frac{\partial}{\partial Q_I}\sum_{J=1}^N\mathrm{tr}\left[\left(P_J\otimes Q_J\right)Y^TRX\right]&=&\frac{\partial}{\partial Q_I}\sum_{J=1}^N\sum_{i=1}^Np_{iJ}\ \mathrm{tr}\left[Q_JC_{Ji}\right]\\
&=&\sum_{j=1}^Np_{jI}C_{Ij}^T
$$
であるから
$$
\frac{\partial f}{\partial Q_I}&=&\sum_{j=1}^Np_{jI}C_{Ij}^T-4\lambda_Q(h(Q)-Q)\\
\frac{\partial f}{\partial R}&=&XPQY^T-4\lambda_R(h(R)-R)-2\lambda_{|R|}(\det(R)-1)\mathrm{adj}(R)\\
$$
と定まる。ここで$\mathrm{adj}(R)$ は $R$ の余因子行列を表す。未定乗数についても同様に
$$
\frac{\partial f}{\partial \lambda_P}&=&g(P)\\
\frac{\partial f}{\partial \lambda_Q}&=&\sum_{I=1}^Ng(Q_I)\\
\frac{\partial f}{\partial \lambda_R}&=&g(R)\\
\frac{\partial f}{\partial \lambda_{|R|}}&=&(\det(R)-1)^2\\
$$
と計算される。
