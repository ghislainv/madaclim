#!/bin/sh
echo "Executing baobab.code.R script in the background"
Rscript --vanilla madaclim.R > madaclim.log 2>&1 &
echo "Check the progress with command 'tail -f baobab.log'"
echo "Check the processor usage with command 'top'"
## End of script