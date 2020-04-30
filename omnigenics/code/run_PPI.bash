inweb_mat="../../data/inweb_brain.ALL.txt"

dstr=`cat outputs/core_genes.known.for_distance.txt | tr '\n' ';' | tr '\t' '='`

Rscript extract_PPI_network_distance.R ../../data/inweb_brain.ALL.txt outputs/PPI.distances.txt $dstr

python calculate_phi.py outputs/PPI.distances.txt outputs/core_genes.hypothetical.for_phi.txt outputs/PPI.phi.txt --core_exclusions outputs/core_genes.known.for_distance.txt
