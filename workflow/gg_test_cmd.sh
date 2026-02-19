#!/usr/bin/env bash

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "" 1 1

########## Testing conda env base ##########

conda activate base
echo "$(date): Testing conda environment: base"

rpackages=( ape Rcpp edgeR ggtree ggimage phytools )
gg_test_r_packages "${rpackages[@]}"

pypackages=( Bio cython ete4 matplotlib numpy pandas scipy sqlalchemy )
gg_test_python_packages "${pypackages[@]}"

programs=(\
"amalgkit --version"
"AMAS.py -h"
"blastn -version"
"tblastn -version"
"cdskit -h"
"cdskit backalign -h"
"clipkit --version"
"csubst --version"
"diamond help"
"fastp --version"
"env OMPI_MCA_plm=isolated OMPI_MCA_plm_rsh_agent=/bin/false generax -h"
"hyphy --version"
"iqtree --version"
"kallisto version"
"jellyfish --version"
"mafft --version"
"cdskit maxalign -h"
"meme -version"
"nwkit -h"
"nwkit prune -h"
"orthofinder -h"
"env PYMOL_HEADLESS=1 QT_QPA_PLATFORM=offscreen pymol -cq -d quit"
"seqkit version"
"prefetch --version"
"provean"
"busco --version"
"augustus --species=help"
"metaeuk -h"
"prodigal -v"
"hmmsearch -h"
"TransDecoder.LongOrfs --version"
"trimal --version"
"Trinity --version")

gg_test_shell_commands "${programs[@]}"

conda deactivate


########## Testing conda env base (R stack) ##########

conda activate base
echo "$(date): Testing conda environment: base"

rpackages=(
ape
cowplot
ggimage
ggrepel
ggtree
igraph
l1ou
magick
nlme
phangorn
PhylogeneticEM
phytools
Rcpp
rkftools
Rphylopars
)
gg_test_r_packages "${rpackages[@]}"

conda deactivate

echo "$(date): Exiting Singularity environment"
