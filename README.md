# Krona charts in R

This repo provides an R functions for creating [Krona charts](https://github.com/marbl/Krona/wiki) in R.

## Features

- The interactive charts can be explored directly within R-Markdown documents in RStudio and will be embedded into documents rendered to *HTML* webpages.
- For document formats like PDF, DOCX, ODT, etc., snapshot images will be inserted.
- Support for [phyloseq](https://github.com/joey711/phyloseq) objects.
- Flexible grouping into multiple datasets explored within the same chart
- Custom coloring schemes

## Example

The following snapshot image of a Krona chart gives an overview of the taxa in the GlobalPatterns dataset from [phyloseq](https://joey711.github.io/phyloseq/index.html) package. It was generated in the [example analysis](example.md).

![GlobalPatterns Krona chart](krona/global_patterns.png)

## How to use

An introduction [is found here](example.md). The original R-Markdown document used to generate this webpage is [example.Rmd](example.Rmd). The rendered PDF output of the same document [is found here](example.pdf). For a documentation currently refer to the [R script](embed_krona.R). 
