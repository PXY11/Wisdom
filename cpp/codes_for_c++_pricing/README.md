compile command exsample

[SimpleMCMain1.cpp]
g++ mains/SimpleMCMain1.cpp  source/Random1.cpp -I include/ -o SimpleMCMain -g -Og

[SimpleMCMain2.cpp]
g++  mains/SimpleMCMain2.cpp  source/PayOff1.cpp source/SimpleMC.cpp source/Random1.cpp -I include/ -o SimpleMCMain2 -g -Og

