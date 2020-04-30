#!/bin/bash

for base in RC_Amy RC_Hippo RC_Neuron RC_NPC; do
	echo python calculate_phi.py outputs/${base}.distances.txt outputs/core_genes.hypothetical.for_phi.txt outputs/${base}.phi.txt --core_exclusions outputs/core_genes.known.for_distance.txt
	python calculate_phi.py outputs/${base}.distances.txt outputs/core_genes.hypothetical.for_phi.txt outputs/${base}.phi.txt --core_exclusions outputs/core_genes.known.for_distance.txt

done
