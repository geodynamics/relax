
# check out input.f90 for the correct version number

# make macports release source file
cd release/ports/
tar -cvjf relax3d-1.0.4.tar.bz2 relax3d-1.0.4

# upload relax3d-1.0.4.tar.bz2 to 
http://www.geodynamics.org/cig/software/relax/macports/relax3d-1.0.4.tar.gz

# checksum the tar ball and copy the number in the package info file
shasum -a 256 relax3d-1.0.4.tar.bz2
openssl rmd160 relax3d-1.0.4.tar.bz2

# copy the package.info file to fink
sudo cp relax.info /sw/fink/dists/local/main/finkinfo/

# validate the package
fink validate /sw/fink/dists/local/main/finkinfo/relax.info

# build the package
fink -m --build-as-nobody rebuild relax

# check out the content of the package
dpkg --contents /sw/fink/dists/local/main/binary-darwin-x86_64/relax_1.0.4-1_darwin-x86_64.deb 

# check the details of the fink package
dpkg --info /sw/fink/dists/local/main/binary-darwin-x86_64/relax_1.0.4-1_darwin-x86_64.deb
