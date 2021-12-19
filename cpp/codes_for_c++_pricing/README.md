compile command example

[SimpleMCMain1.cpp]
g++ mains/SimpleMCMain1.cpp  source/Random1.cpp -I include/ -o SimpleMCMain -g -Og

[SimpleMCMain2.cpp]
g++  mains/SimpleMCMain2.cpp  source/Random1.cpp source/PayOff1.cpp source/SimpleMC.cpp -I include/ -o SimpleMCMain2 -g -Og

[SimpleMCMain3.cpp]
g++  mains/SimpleMCMain3.cpp  source/Random1.cpp source/PayOff2.cpp source/SimpleMC2.cpp -I include/ -o SimpleMCMain3 -g -Og

[SimpleMCMain4.cpp]
g++  mains/SimpleMCMain4.cpp  source/Random1.cpp source/PayOff2.cpp source/SimpleMC2.cpp -I include/ -o SimpleMCMain4 -g -Og

[SimpleMCMain5.cpp]
g++  mains/SimpleMCMain5.cpp  source/Random1.cpp source/PayOff2.cpp source/SimpleMC2.cpp source/DoubleDigital.cpp -I include/ -o SimpleMCMain5 -g -Og

[SimpleMCMain5.cpp]
g++  mains/SimpleMCMain5.cpp source/DoubleDigital.cpp -I include/ -o SimpleMCMain5 -g -Og