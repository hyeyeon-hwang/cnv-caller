# cnv-caller
This copy number variation (CNV) caller detects large-scale variation in low coverage aligned whole genome bisulfite sequencing (WGBS) data. The CNV caller uses a fixed bin approach.

**cnv_bins.py** uses the fixed bin approach and takes the following input arguments:
- --controls: path to a directory of bam files of control samples
- --experimental: path to a directory of bam files of experimental samples
- --window: size of CNV bin
- --sexCallerOnly: either "yes" or "no" to indicate whether to execute the sex caller only

**cnv_plots_shiny_app.R** is an interactive Shiny app that plots the copy number value estimations of each chromosome bin and sample. Users may choose which chromosome and sample to display in the plots.
