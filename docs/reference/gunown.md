# Gun Ownership Data Set

The data for this example are from a 2015 survey of 321 American voters
conducted by Smithson. Participants were asked to assign a number
representing their feeling about gun ownership from 0 (very negative) to
100 (very positive). They were also asked to choose their political
affiliation from four categories: Democrat, Republican, independent, and
no preference.

## Usage

``` r
gunown
```

## Format

A data frame with these columns:

- poliorient:

  Political affiliation by assignment to one of four categories

- demo:

  A binary variable, coded 1 = Democrat, 0 = any other category

- indep:

  A binary variable, coded 1 = independent, 0 = any other category

- nopref:

  A binary variable, coded 1 = no preference, 0 = any other category

- repub:

  A binary variable, coded 1 = Republican, 0 = any other category

- gun01:

  A numeric variable, linearly transformed from the 0-100 range to 0-1

## Source

Smithson, M. 2015. Unpublished dataset.
