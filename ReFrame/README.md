# ReFrame

To run *ReFrame*, place your input files (format details below) into the `in` folder of your working directory,
and save a file containing the reference relative peak heights (median of relative peak heights per marker, details below)
under the name `reflist.xlsx`.
Then proceed to run *ReFrame*

```

Rscript reframe.R

```

## Input File Format

The directory *ReFrame* is run from should contain a `reflist.xlsx` file containing reference relative peak heights for each marker of interest (details below),
an `in` directory containing the input data for each marker and an `out` directory to write outputs to.

### Definitions
|Name|Definition|
|----|----------|
| MP | main peak - the peak of the expected fragment length (K)|
| LX | minus-X peak - the peak of length K - X|
| RX | plus-X peak - the peak of length K + X|
| relative peak height | peak-height / sum of peak-heights|
| reference relative peak height | median of peak-height for a given marker|


### `reflist.xlsx`
`xlsx` format file containing a single sheet with the following layout:

| Gene ID | Run ID |  L4 | L3 | L2 | L1 | MP | R1 | R2 | R3 |
|----------|--------|----|----|----|----|----|----|----|----|
| gene-or-marker-id | *ignored* | median-relative-peak-height-L4| median-relative-peak-height-L3| median-relative-peak-height-L2| median-relative-peak-height-L1| median-relative-peak-height-MP| median-relative-peak-height-R1| median-relative-peak-height-R2| median-relative-peak-height-R3| median-relative-peak-height-R4|
|...|...|...|...|...|...|...|...|...|...|

This file contains *reference relative peak heights*, which can be computed by taking the median of relative peak heights from the data of a particular marker or gene of interest.
 
### Files in the `in` Folder

Inputs should be one `xlsx` format file per gene-of-interest with a sheet named "Heights" (case sensitive!) with the following layout:

| Tumor ID | Run ID | L4 | L3 | L2 | L1 | MP | R1 | R2 | R3 |
|----------|--------|----|----|----|----|----|----|----|----|
| *ignored*|*ignored*|L4-size| L3-size| L2-size| L1-size| MP-size| R1-size| R2-size| R3-size|
| tumor-identifier | run-identifier | peak-height L4| peak-height L3| peak-height L2| peak-height L1| main-peak-height| peak-height R1| peak-height R2| peak-height R3|
| tumor-identifier| run-identifier |...|...|...|...|...|...|...|...|

with one row in the format of row 3 for each tumor sample.

All other sheets in the `xlsx` file are ignored.

## Output File Format

Outputs for each marker of interest are written in separate files to the `out` folder. Files contain the original input data,
as well as a sheet "Result_of_algorithm" containing the processed outputs.
