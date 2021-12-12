#include <string>
#include <json/json.h>
#include <iostream>
#include <fstream>
#include <typeinfo> //判断类型
#include <vector>

using namespace std;

void readStrJson(); //从字符串中读取Json数据（简单Json格式)
//void readStrProJson(); //同样从字符串中读取Json数据（Json格式更加复杂）

int main(int argc,char** argv)
{
    cout<<"sample case:"<<endl;
    readStrJson();
    //cout<<"complex case:"<<endl;
    //readStrProJson();

    return 0;
}

void readStrJson()
{
    //字符串  -----  注意：使用""作为换行也可以，实际写成一行也可以"{\"praenomen\":\"Gaius\",\"nomen\":\"Julius\",\"cognomen\":\"Caezar\",\"born\":-100,\"died\":-44}"; 
    const char* str =   
      "{\"praenomen\":\"Gaius\","
      "\"nomen\":\"Julius\",\"cognomen\":\"Caezar\","
      "\"born\":-100,\"died\":-44,\"final\":0.66}";  

    /* 
    //json内容如下： key:val其中val为字符串、数字等基础类型
    { 
        "praenomen":"Gaius", 
        "nomen":"Julius", 
        "cognomen":"Caezar", 
        "born":-100, 
        "died":-44  
    } 
    */

    Json::Reader reader; //用于读取字符串
    Json::Value root; //用于存放结果到map中
    if(reader.parse(str,root))
    {
        string praenomen = root["praenomen"].asString();
        string nomen = root["nomen"].asString();
        string cognomen = root["cognomen"].asString();
        int born = root["born"].asInt();
        int died = root["died"].asInt();
        double final = root["final"].asDouble();
        //cout<<root["born"]<<typeid(root["born"]).name()<<endl; //封装为了一个对象N4Json5ValueE
        cout<<praenomen+" "+nomen+" "+cognomen<<" was born in year "<< born
        <<", died in year "<<died<<", final "<<final<<endl;
    }
}