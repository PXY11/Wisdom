This project is for reimplementing the MATLAB code by C++,

which is used for pricing convertible bond


------------------------------------------------------Compile Command------------------------------------------------------

First change directory to the convertible_bond

Then copy the following commands to terminal to compile

g++ main/testJson.cpp ./lib_json/*.cpp -I ./ -o testJson

g++ main/testReadJsonFile.cpp ./lib_json/*.cpp -I ./ -o testReadJsonFile

g++ main/testCBP.cpp src/Parameter.cpp lib_json/*.cpp -I ./ -I include/ -o testCBP