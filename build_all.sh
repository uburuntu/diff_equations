# Usage
# ./build_all.sh
# or
# ./build_all.sh clean

cd ./laspack
make $1 > build.log
cd ..
make $1 >> build.log
cd ./eigen/arpack
make $1 >> build.log
cd ..
make $1 >> build.log
