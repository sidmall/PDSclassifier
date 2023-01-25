
# PDSclassifier

<!-- badges: start -->
<!-- badges: end -->

*PDSclassifier* R package provides a pathway-based molecular classification system for colorectal cancer (CRC), which can be applied to gene expression profiles to stratify into three PDS (Pathway-Derived Subtype): PDS1, PDS2 and PDS3, with distinct molecular biology.

## Installation

You can install the development version of PDSclassifier like so:

``` r
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github('sidmall/PDSclassifier')
```

## Usage

An example where using the test dataset from the R package, PDS classification can be made with `PDSpredict()` function.
``` r
library(PDSclassifier)
pds.calls <- PDSpredict(testData, species = 'human', threshold = 0.6)
```

*PDSclassifier* can be applied to both human and mouse transcriptomic data with parameter:
`species = c("human", "mouse")`.

The default prediction probability `threshold = 0.6` has been set. It can be altered anywhere between 0 (less stringent) to 1 (very stringent). However, recommendation is be stay between 0.5-0.7 to retain enough samples without losing underlying biology that defines PDS.
