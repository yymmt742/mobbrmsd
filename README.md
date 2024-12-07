<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>


<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/yymmt742/mobbrmsd">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

<h3 align="center">mobbrmsd</h3>

  <p align="center">
    molecular-oriented RMSD with the branch-and-bound.
    <br />
<!--
    <a href="https://github.com/yymmt742/mobbrmsd"><strong>Explore the docs »</strong></a>
    <br />
-->
    <br />
    <a href="https://github.com/yymmt742/mobbrmsd/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/yymmt742/mobbrmsd/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#background">Background</a></li>
    <li><a href="#benchmark">Benchmark</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#reference">Reference</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
# About The Project

[![CI](https://github.com/yymmt742/mobbrmsd/actions/workflows/ci.yml/badge.svg)](https://github.com/yymmt742/mobbrmsd/actions/workflows/ci.yml)
[![Create Release Branch](https://github.com/yymmt742/mobbrmsd/actions/workflows/create_release.yml/badge.svg)](https://github.com/yymmt742/mobbrmsd/actions/workflows/create_release.yml)

Root mean squared deviation (RMSD) is one of the most common metrics
for comparing the similarity of three-dimensional chemical structures.
The chemical structuresimilarity plays an important role in data chemistry
because it is closely related tochemical reactivity, physical property, and bioactivity.

The molecular-oriented RMSD with branch-and-bound (mobbRMSD) is an RMSD-based metric for 3D chemical structure similarity.
mobbRMSD is formulated in molecular-oriented coordinates
and uses the branch-and-bound method to obtain an exact solution.
It can handle large and complex chemical systems such as molecular liquids, solvationsof solute, and self-assembly of large molecules,
which are difficult to handle using conventional methods.

For a molecular coordinate pair $\mathbf{X}$ and $\mathbf{X}'$ that consisting of $M$ molecules of $n$ atoms,
the moRMSD is deﬁned as follows:
```math
  \text{moRMSD}\left(\mathbf{X},\mathbf{X}'\right)
  =
  \min_{\mathbf{R},c,\mu,\nu_{1},\nu_{2},\ldots,\nu_{M}}
  \sqrt{\frac{1}{Mn} \sum_{I=1}^{M} \sum_{j=1}^{n} \left\|{x}_{j,I}-\mathbf{R}{x}'_{\nu_{I}(j),\mu(I)}-{c}\right\|^2},
```
where, $x_{j,I}$ and ${x'}_{j,I}$ are the Cartesian coordinates corresponding to $j$-th atom of $I$-th molecule, respectively,
$c$ is a translation vecotor, and $\mathbf{R}$ is a rotation matrix.
$\nu\_{I}$ and $\mu$ are permutations on $\\{1,\ldots,n\\}$ and $\\{1,\ldots, M\\}$, respectively.
$\nu\_{I}$ takes the appropriate domain of definition corresponding to the molecular topology.
Since $\nu\_{I}$ and $\mu$ expand the solution space by factorial and exponential costs with respect to $M$, respectively,
it becomes more difficult to find a solution by brute force when $M$ is large.

mobbRMSD practically eliminates this difficulty by using the branch-and-bound method.
See Back Ground and Benchmark for details.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

* gfortran >= 9.4.0
* OpenBLAS (optional)
* OpenMP (optional)

To use the Python interface, you additionally need the following:
* python >= 3.8
* pip

### Installation

   You can use package build via
   ```sh
   pip install git+ssh://git@github.com/yymmt742/mobbrmsd.git
   ```

   Running demonstrations via
   ```sh
   python -m mobbrmsd demo
   ```

   Running via
   ```sh
   python -m mobbrmsd run input
   ```
   Input format is ...

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Usage

<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Background


<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Benchmark

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Roadmap

- [ ] Add Usage
- [x] Enable autovariance sorting
- [ ] Enable skip tree
- [ ] Compatible with compilers (intel)
- [x] Compatible with compilers (nv)
- [ ] Add detail documentation
- [ ] Add detail documentation (Python interface)
- [ ] Add benchmarks
- [x] Internalize lapack

See the [open issues](https://github.com/yymmt742/mobbrmsd/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

This project is open source and we invite contributions.
If you have a suggestion that would make this better,
please fork the repo and create a pull request.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the MIT License.
See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

YYMMT742 - yymmt@kuchem.kyoto-u.ac.jp

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- Reference -->
## Reference

Molecular superposition is based on the following algorithm

* [A solution for the best rotation to relate two sets of vectors](https://scripts.iucr.org/cgi-bin/paper?S0567739476001873)
* [Rapid calculation of RMSDs using a quaternion-based characteristic polynomial](https://scripts.iucr.org/cgi-bin/paper?S0108767305015266)
* [RMSD and Symmetry](https://onlinelibrary.wiley.com/doi/10.1002/jcc.25802)

The solution to the linear assignment problem for estimating the variational lower bound is based on the Hungarian method.

* [The Hungarian method for the assignment problem](https://onlinelibrary.wiley.com/doi/10.1002/nav.3800020109)
* [Algorithms for the Assignment and Transportation Problems](http://www.jstor.org/stable/2098689)


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

