# wispack
 Rcpp implementation of warped-sigmoid Poisson-process mixed-effect models

# Installation 

1. Clone this git repo
2. Make "build_install.sh" executable by running "chmod +x build_install.sh" in bash terminal
3. In terminal, run "./build_install.sh" to build and install the package
4. To ensure a clean start, run: 

rm -f src/*.o src/*.so
rm -rf wispack.Rcheck
rm -f wispack_*.tar.gz
./build_install.sh