
# PDSclassifier

<!-- badges: start -->
<!-- badges: end -->

PDSclassifier R package provides an ontolgoy-based molecular classification system for colorectal cancer, which can be applied to gene expression profiles to stratify into three PDS (Pathway-derived subtypes): PDS1, PDS2 and PDS3, with distinct molecular biology.

## Installation

You can install the development version of PDSclassifier like so:

``` r
devtools::install_github('sidmall/PDSclassifier')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PDSclassifier)
pds.calls <- PDSpredict(testData, species = 'human', threshold = 0.6)
```

