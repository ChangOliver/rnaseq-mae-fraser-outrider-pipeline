apt update
apt install -y --no-install-recommends software-properties-common dirmngr
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt update
apt install -y --no-install-recommends r-base
apt install -y libssl-dev libcurl4-openssl-dev libxml2-dev libboost-all-dev
Rscript install.R

