compile command example

[SimpleMCMain1.cpp]
g++ mains/SimpleMCMain1.cpp  source/Random1.cpp -I include/ -o SimpleMCMain -g -Og

[SimpleMCMain2.cpp]
g++  mains/SimpleMCMain2.cpp  source/Random1.cpp source/PayOff1.cpp source/SimpleMC.cpp -I include/ -o SimpleMCMain2 -g -Og

[SimpleMCMain3.cpp]
g++  mains/SimpleMCMain3.cpp  source/Random1.cpp source/PayOff2.cpp source/SimpleMC2.cpp -I include/ -o SimpleMCMain3 -g -Og