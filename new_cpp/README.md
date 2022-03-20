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
g++ testConvbond.cpp Parameter.cpp PDESolver.cpp  -I ./eigen-3.4.0  -o3 -o testConvbond


eigen库使用方法参考

http://www.javashuo.com/article/p-fjkcsfmi-ey.html

https://eigen.tuxfamily.org/index.php?title=Main_Page

https://www.cnblogs.com/goingupeveryday/p/5699053.html

https://www.cnblogs.com/rainbow70626/p/8819119.html


g++ sp.cpp -I ./eigen-3.4.0 -o sp