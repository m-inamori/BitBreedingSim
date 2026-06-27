# BitBreedingSim

Fast breeding simulator package for R

## Installation

Install the package directly from GitHub.

```r
install.packages("remotes")  # if not already installed
remotes::install_github("m-inamori/BitBreedingSim", subdir = "R_package")
```

Alternatively, if you use `pak`:

```r
install.packages("pak")  # if not already installed
pak::pkg_install("m-inamori/BitBreedingSim/R_package")
```

## Building from source

Clone the repository and install the package using the standard R tools.

```bash
git clone https://github.com/m-inamori/BitBreedingSim.git
cd BitBreedingSim/R_package
R CMD INSTALL .
```
## Install (Windows)

On Windows, building from source requires Rtools.

1. Install Rtools (matching your R version)
2. Confirm availability:

```r
Sys.which("make")
```

3. Install dependencies:

```r
install.packages("Rcpp")
```

4. Install the package:

```r
install.packages("path/to/BitBreedingSim/R_package", repos = NULL, type = "source")
```


## Usage

```
> library(BitBreedingSim)
> info <- create_base_info(seed=2)
> add_trait_A(info, name="Trait1", mean=100, h2=0.6, sd=20, num_loci=10)
> origins <- create_origins(10, info, "orig_")
> prog <- cross_randomly(1000, origins, origins, "prog_")
> geno <- get_genotypes(prog)
> geno[1:4, 1:4]
       marker00000001 marker00000002 marker00000003 marker00000004
prog_1              0              0              0              0
prog_2              1             -1              0              0
prog_3              1              0              1              1
prog_4             -1              1              1              1
> pheno <- get_phenotypes(prog, 1)
i : 1 num_traits : 1
> pheno[1:10]
 [1] 116.26303  81.35491  80.49995 138.49164 141.23155 127.53185 128.23878
 [8]  95.75027  98.02536 110.55966
```

## License
MIT License

Copyright (c) 2024 Minoru Inamori

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
