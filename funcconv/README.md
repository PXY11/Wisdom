
编译命令

g++ src/main.cpp src/func1_initialize.cpp src/func1_terminate.cpp src/func1.cpp  -I ./ -I include/ -o func1main

g++ src/main.cpp src/*.cpp  -I ./ -I include/ -o funcconvmain