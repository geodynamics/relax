
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
mkdir -p ~/ports/science/relax3d
cp science/relax3d/Portfile ~/ports/science/relax3d

# update the repository list
vi /opt/local/etc/macports/sources.conf

# and insert the lines
# 
#   file:///Users/whoami/ports

# update macports
sudo port -v selfupdate

# update the portfile
cd ~/ports
sudo portindex

# check for the package
port search relax3d

# install from the local repository
sudo port install relax3d@1.0.5
