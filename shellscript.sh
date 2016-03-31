#!/bin/sh
echo "Executing baobab.code.R script in the background"
Rscript --vanilla environ.R FALSE > environ.log 2>&1 &
# If you wan't to download the source data, do:
# Rscript --vanilla environ.R TRUE > environ.log 2>&1 &
echo "Check the progress with command 'tail -f environ.log'"
echo "Check the processor usage with command 'top'"
## End of script