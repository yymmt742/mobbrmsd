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
