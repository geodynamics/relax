
# make fink release source file
cd release/fink/
tar -cvzf relax-1.0.4.tar.gz relax-1.0.4

# upload relax-1.0.4.tar.gz to 
http://www.geodynamics.org/cig/software/relax/fink/relax-1.0.4.tar.gz

# checksum the tar ball and copy the number in the package info file
md5sum relax-1.0.4.tar.gz

# copy the package.info file to fink
sudo cp relax.info /sw/fink/dists/local/main/finkinfo/

# validate the package
fink validate /sw/fink/dists/local/main/finkinfo/relax.info

# build the package
fink -m --build-as-nobody rebuild relax

# check out the content of the package
dpkg --contents /sw/fink/dists/local/main/binary-darwin-x86_64/relax_1.0.4-1_darwin-x86_64.deb 

