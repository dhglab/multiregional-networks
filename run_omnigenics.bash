set -e -x

cat omnigenics/data.tar.gz01 omnigenics/data.tar.gz02 omnigenics/data.tar.gz03 > omnigenics/data.tar.gz
tar xvzf omnigenics/data.tar.gz

mkdir -p omnigenics/outputs

if [ ! -e 'omnigenics/outputs/InWeb.distances.txt' ]; then
  python omnigenics/code/data_preprocessing.py
fi

python omnigenics/code/calculate_phi.py omnigenics/outputs/RC_Neuron.distances.txt omnigenics/outputs/core_genes.hypothetical.for_phi.txt omnigenics/outputs/RC_Neuron.phi.txt --core_exclusions omnigenics/outputs/core_genes.known.for_distance.txt

python omnigenics/code/calculate_phi.py omnigenics/outputs/WHOLE_BRAIN.distances.txt omnigenics/outputs/core_genes.hypothetical.for_phi.txt omnigenics/outputs/WHOLE_BRAIN.phi.txt --core_exclusions omnigenics/outputs/core_genes.known.for_distance.txt

python omnigenics/code/calculate_phi.py omnigenics/outputs/InWeb.distances.txt omnigenics/outputs/core_genes.hypothetical.for_phi.txt omnigenics/outputs/InWeb.phi.txt  --core_exclusions omnigenics/outputs/core_genes.known.for_distance.txt
