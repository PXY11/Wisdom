This project is for reimplementing the MATLAB code by C++,

which is used for pricing convertible bond


------------------------------------------------------Compile Command------------------------------------------------------

First change directory to the convertible_bond

Then copy the following commands to terminal to compile

[testJson.cpp]
g++ main/testJson.cpp ./lib_json/*.cpp -I ./ -o testJson

[testReadJsonFile.cpp]
g++ main/testReadJsonFile.cpp ./lib_json/*.cpp -I ./ -o testReadJsonFile

[testConvbond.cpp]
g++ main/testConvbond.cpp src/Parameter.cpp src/PDESolver.cpp lib_json/*.cpp -I ./ -I include/ -o testConvbond

g++ main/testConvbond.cpp  -o testConvbond

g++ testConvbond.cpp Parameter.cpp -I ./  -o testConvbond

g++ testConvbond.cpp Parameter.cpp PDESolver.cpp  -I ./  -o testConvbond