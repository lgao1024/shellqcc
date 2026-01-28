#!/bin/bash


# Extract Final Gibbs Free Energy from ORCA frequency calculation output files
grep "Final Gibbs free energy" */*.out | sed -E 's/[^/]+\/([^:]+)\.freq\.out:.*\.\.\.[[:space:]]+(-?[0-9]+\.[0-9]+)[[:space:]]+Eh/\1 \2/' > Gibbs.txt

