#include"interpfun.h"
Interp::Interp()
{
 strMethod1 ="lgr";
 strMethod2 = "lg3";
 strMethod3 = "spl";
}

Interp::~Interp()
{
 
}

float Interp::lgr(float *x,float *y,int n,float t)
{
    int i,j,k,m;
    float z,s;
    z= 0.0f;
    if(n<1)
        return(z);
    if(n==1)
    {
        z=y[0];
        return(z);
    }

    if(n==2)
    {
        z=(y[0]*(t-x[1]) - y[1]*(t-x[0]))/(x[0]-x[1]);
        return(z);
    }

    if(n>1 && t>x[n-1])
    {
        z= y[n-1];
        return(z);
    }

    if(n>1 && t<x[0])
    {
        z= y[0];
        return(z);
    }

    i =0;
    while((x[i] <t) && (i<n)) 
    {
        i=i+1;
    }

    k = i-4;
    if(k<0) 
        k=0;
        m = i+3;
    if(m>n-1) 
        m = n-1;
    for( i = k;i<=m;i++)
    {
        s=1.0;
        for( j =k;j<=m;j++)
        {
            if(j!=i) 
                s= s*(t-x[j])/(x[i]-x[j]);
        }
        z = z+s*y[i];
    }
    return(z);
}

/*
版权声明：本文为CSDN博主「最大的素数」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
原文链接：https://blog.csdn.net/nneerr123/article/details/79296385
*/