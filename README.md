# AR(1) Biascorrection

Modified biascorrection by Sørbye

!The package is under active development!

## Installation

```{r}
devtools::install_github("https://github.com/michael-mueller-uibk/ar1MedianBiascorrection")
```

## Usage

```{r}
library("ar1MedianBiascorrection")
ar1_bias_corr(phi = 0.8, n = 10, method = "yw")
```
