# Theory

## Molecular-oriented RMSD

The molecular-oriented RMSD is defined for molecular assemblies
containing multiple molecular species.

For two molecular-oriented coordinate sets \(X\) and \(X'\),

\[
\operatorname{moRMSD}(X,X')
=
\min_{\mathbf R,c,\mu,\nu}
\sqrt{
\frac{1}{N}
\sum_i
\left\|
x_i-\mathbf R x'_{\nu(i),\mu(i)}-c
\right\|^2
}.
\]

The minimization includes:

- global rotation,
- translation,
- permutation of equivalent molecules,
- intramolecular symmetry permutations.

## Branch-and-bound

A direct enumeration of molecular permutations becomes prohibitively
expensive as the number of molecules increases.

`mobbrmsd` uses a branch-and-bound algorithm to eliminate portions of
the permutation search space that cannot contain the optimum.

See [mobbRMSD](mobbrmsd.md) for more details.
