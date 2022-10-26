sudo apt-get install libboost-all-dev libeigen3-dev libgmp-dev
cd /usr/include
sudo ln -sf eigen3/Eigen Eigen
sudo ln -sf eigen3/unsupported unsupported
# upgrade cmake to 3.22
apt remove cmake
cd /opt
wget https://github.com/Kitware/CMake/releases/download/v3.22.3/cmake-3.22.3-linux-x86_64.sh
chmod +x ./cmake-3.22.3-linux-x86_64.sh
bash ./cmake-3.22.3-linux-x86_64.sh
sudo ln -s /opt/cmake-3.22.3-linux-x86_64/bin/* /usr/bin
