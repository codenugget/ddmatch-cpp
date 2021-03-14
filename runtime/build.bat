pushd .
cd ..
if not exist build (
	mkdir build
)
cd build
if not exist ddmatch.sln (
	cmake .. -D_CRT_SECURE_NO_WARNINGS=1
)
popd
