pushd .
cd ..
if not exist build (
	mkdir build
)
cd build
if not exist ddmatch.sln (
	cmake ..
)
popd
