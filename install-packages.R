install.packages('pacman', repos='cran.wustl.edu')
install.packages('BiocManager', repos='cran.wustl.edu')
library(pacman)

pkgs <- scan('Rpackages.txt', what='%s')
for ( pkg in pkgs ) {
  pacman::p_install(pkg, character.only=T, update.bioconductor=T)
}

install.packages('WGCNA', repos='cran.wustl.edu')
