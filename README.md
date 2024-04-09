<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
<!--
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
-->



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/yymmt742/mobbrmsd">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

<h3 align="center">mobbrmsd</h3>

  <p align="center">
    molecular-oriented branch-and-bound for RMSD.
    <br />
    <a href="https://github.com/yymmt742/mobbrmsd"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/yymmt742/mobbrmsd">View Demo</a>
    ·
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
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#Reference">Reference</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

[![Product Name Screen Shot][product-screenshot]](https://github.com/yymmt742/mobbrmsd/blob/main/technical_notes.md)

Calculate mobbrmsd.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



### Built With

* [![cmake][cmake]][cmake-url]
* [![fortran][fortran-shield]][fortran-url]
* [![python][python-shield]][python-url]
* [![setuptools][setuptools-shield]][setuptools-url]
* [![scikit-build][skbuild-shield]][skbuild-url]

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

* gfortran >= 9.4.0
* cmake >= 3.9
* openmp
* blas/lapack

To use the Python interface, you additionally need the following:
* python >= 3.8
* pip
* setuptools >= 42
* numpy >= 1.21
* f2py
* skbuild
* wheel

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/yymmt742/mobbrmsd
   ```
2. Build fortran library
   ```sh
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=[Release, Debug, ...] -DCMAKE_PREFIX_PATH=your_local_cmake_path
   make install
   ```
   It doesn't work with MKL. Please use
   ```sh
   -DBLA_VENDOR=[OpenBLAS, ATLAS, ...]
   ```
   instead.
3. Build python interface.
   ```sh
   pip install .
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

Usage

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap

- [ ] Add Usage
- [ ] Compatible with compilers other than gnu
- [ ] Add detail documentation
- [ ] Add benchmarks

See the [open issues](https://github.com/yymmt742/mobbrmsd/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

YYMMT742 - yymmt@kuchem.kyoto-u.ac.jp

Project Link: [https://github.com/yymmt742/mobbrmsd](https://github.com/yymmt742/mobbrmsd)

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
[contributors-shield]: https://img.shields.io/github/contributors/yymmt742/mobbrmsd.svg?style=for-the-badge
[contributors-url]: https://github.com/yymmt742/mobbrmsd/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/yymmt742/mobbrmsd.svg?style=for-the-badge
[forks-url]: https://github.com/yymmt742/mobbrmsd/network/members
[stars-shield]: https://img.shields.io/github/stars/yymmt742/mobbrmsd.svg?style=for-the-badge
[stars-url]: https://github.com/yymmt742/mobbrmsd/stargazers
[issues-shield]: https://img.shields.io/github/issues/yymmt742/mobbrmsd.svg?style=for-the-badge
[issues-url]: https://github.com/yymmt742/mobbrmsd/issues
[license-shield]: https://img.shields.io/github/license/yymmt742/mobbrmsd.svg?style=for-the-badge
[license-url]: https://github.com/yymmt742/mobbrmsd/blob/master/LICENSE.txt
[product-screenshot]: images/screenshot.png
[cmake]: https://img.shields.io/badge/Cmake-064F8C?style=for-the-badge&logo=cmake&logoColor=EEEEEE
[cmake-url]: https://cmake.org/
[fortran-shield]: https://img.shields.io/badge/Fortran-734F96?style=for-the-badge&logo=fortran&logoColor=FFFFFF
[fortran-url]: https://fortran-lang.org/
[python-shield]: https://img.shields.io/badge/python-3776AB?style=for-the-badge&logo=python&logoColor=FFFFFF
[python-url]: https://www.python.org/
[setuptools-shield]: https://img.shields.io/badge/setuptools-3775A9?style=for-the-badge&logo=pypi&logoColor=FFFFFF
[setuptools-url]: https://pypi.org/project/setuptools/
[skbuild-shield]: https://img.shields.io/badge/skbuild-35495E?style=for-the-badge
[skbuild-url]: https://scikit-build.readthedocs.io/en/latest/

