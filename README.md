# Krona charts from phyloseq objects in R(-Markdown)

This repo provides an R function, which allows creating [Krona charts](https://github.com/marbl/Krona/wiki)
from a [phyloseq object](https://github.com/joey711/phyloseq).

## Features

- The interactive charts can be embedded R-Markdown and thus directly explored from
  within RStudio
- If creating HTML output files using knitr (the Knit button in RStudio), the charts
  are embedded interactively as well. Embedding several Krona charts is in principle 
  possible, although it leads to pretty large HTML files.
- Embedding in PDF and DOCX (and possibly other formats) is possible, 
  currently this is solved with some Javascript hacks.
- Using a custom color scheme is enabled through a custom import script (`krona_import.pl`)

## How to use

An example R-Markdown document based on the GlobalPatterns dataset in phyloseq
is found [in this repository](phyloseq_example.Rmd), a compiled PDF snapshot is
[found here](phyloseq_example.pdf).

For more documentation refer to the [R script](embed_krona.R).

It is possible to supply an *install\_dir* argument to *plot\_krona*, which 
expects the path to the *embed\_krona* directory, and in which the KronaTools 
directory needs to be placed manually 
(download here: https://github.com/marbl/Krona/releases/latest).
This approach is required if using custom color schemes, the `krona_import.pl` script
requires access to the KronaTools directory.

Creating snapshots requires an up-to-date Krona version (tested with Krona version 2.8.1).
