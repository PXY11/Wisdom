#include"PDESolver.h"
#include<iostream>
#include<cmath>
#include <Eigen/Dense>  //https://www.cnblogs.com/rainbow70626/p/8819119.html 使用简介
#include <Eigen/LU>
#include <ctime>
using namespace Eigen;
using namespace std;


void PDESolver::set_iter(int iter_time)
{
    this->iteration = iter_time;
}

void PDESolver::show_iter()
{
    cout<<"iter time = "<<this->iteration<<endl;
}
double PDESolver::interpCB(MatrixXd x,MatrixXd y,double ind)
{
    //https://zhidao.baidu.com/question/619575742421839812.html MATLAB插值函数原理
    if(x.size()!=y.size()) 
    {
        cout<<"The two list have different length!"<<endl;
        return -1;
    }

    int i=0;
    while (i<ind) i++;
    double ratio = (ind - x(i))/(x(i+1) - x(i));
    double res = y(i) + (y(i+1) - y(i))*ratio;
    return res;
}

double PDESolver::solve()
{
    cout<<"PDESolver::solve() 开始执行"<<endl;
    //边界条件参数不在parameters中设置，在solve时再设置
    MatrixXd i(int(this->pde_param_ptr->Ns+1),1);  //i [1,2...,Ns+1]
    cout<<"i size = "<<i.size()<<endl;
    // MatrixXd i(3,1);
    for(int token=1;token<int(this->pde_param_ptr->Ns+2);token++) i(token-1)=token;
    // cout<<i<<endl;
    cout<<"i 初始化完成"<<endl;

    MatrixXd S(int(this->pde_param_ptr->Ns+1),1); // [0,h,2h...Ns*h]
    cout<<"S size = "<<S.size()<<endl;
    for(int token=0;token<int(this->pde_param_ptr->Ns+1);token++) S(token) = this->pde_param_ptr->ds*token;
    // cout<<S.transpose()<<endl;
    cout<<"S 初始化完成"<<endl;

    double Bc_T = this->pde_param_ptr->Bc_star + this->pde_param_ptr->final_coupon_rate*this->pde_param_ptr->F;
    double Bp_T = this->pde_param_ptr->Bp_star + this->pde_param_ptr->final_coupon_rate*this->pde_param_ptr->F;
    cout<<"Bc_T = "<<Bc_T<<endl;
    cout<<"Bp_T = "<<Bp_T<<endl;

    MatrixXd k_T(S.size(),1);
    cout<<"k_T size = "<<k_T.size()<<endl;
    for(int i=0;i<S.size();i++)
    {
        double tmp = this->pde_param_ptr->F/(this->pde_param_ptr->conversion_price-this->pde_param_ptr->q*ceil(this->pde_param_ptr->T-this->pde_param_ptr->dtime));
        k_T(i) = tmp;
    }
    // cout<<k_T.transpose()<<endl;
    cout<<"k_T 初始化完成"<<endl;

    MatrixXd u(S.size(),1);
    cout<<"u size = "<<u.size()<<endl;
    for (int i =0;i<S.size();i++)
    {
        double tmp_1 = k_T(i)*S(i)+this->pde_param_ptr->coupon_time_rate(0)*this->pde_param_ptr->F;
        double tmp_2 = this->pde_param_ptr->F+this->pde_param_ptr->final_coupon_rate*this->pde_param_ptr->F+this->pde_param_ptr->coupon_time_rate(0)*this->pde_param_ptr->F;
        double tmp_3 = max(tmp_1,tmp_2);
        u(i) = tmp_3;
        
    }
    // cout<<u.transpose()<<endl;
    cout<<"u 初始化完成"<<endl;

    MatrixXd B(int(this->pde_param_ptr->Ns+1),1);
    cout<<"B size = "<<B.size()<<endl;
    for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
    {
        double tmp = this->pde_param_ptr->F + this->pde_param_ptr->final_coupon_rate*this->pde_param_ptr->F;
        B(i) = tmp;
    }
    // cout<<B.transpose()<<endl;
    cout<<"B 初始化完成"<<endl;

    clock_t time_stt;
    int count = 0;
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int count4 = 0;
    int count5 = 0;
    int count6 = 0;
    int flag = 0;
    double t;
    double rate1;
    double rate2;
    double r;
    MatrixXd alpha(S.size()-2,1);
    MatrixXd beta(S.size()-2,1);
    MatrixXd alpha_MuMB(alpha.size()+2,1);
    MatrixXd beta_MuMB(beta.size()+2,1);
    MatrixXd gamma_Mu(alpha.size()+2,1);
    int height = int(this->pde_param_ptr->Ns)+1;
    int weight = int(this->pde_param_ptr->Ns)+1;
    MatrixXd Mu(height,weight);
    MatrixXd gamma_MB(alpha.size()+2,1);
    MatrixXd MB(height,weight);
    double AccI;
    double Bc_n;
    double Bp_n;
    MatrixXd k_n(S.size(),1);
    MatrixXd u_n;
    MatrixXd Muu(k_n.size(),1);
    MatrixXd P1_old(height,weight);
    MatrixXd P2_old(height,weight);
    MatrixXd u_old;
    MatrixXd P1(height,weight);
    MatrixXd P2(height,weight);
    for(int n=1;n<2;n++)
    {
        if(flag==1) break;
        t = n*this->pde_param_ptr->dt;
        cout<<"t = "<<t<<endl;
        rate1 = this->interpCB(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,n*this->pde_param_ptr->dt);
        cout<<"rate1 =  "<<rate1<<endl;
        rate2 = this->interpCB(this->pde_param_ptr->rate_T,this->pde_param_ptr->risk_free_rate,(n-1)*this->pde_param_ptr->dt);
        cout<<"rate2 = "<<rate2<<endl;
        r = -(log(rate2/rate1)/this->pde_param_ptr->dt);
        cout<<"r = "<<r<<endl;

        if ((ceil(t-this->pde_param_ptr->dtime) - ceil(t-this->pde_param_ptr->dt-this->pde_param_ptr->dtime)) == 1)
        {//此段代码对应MATLAB源码中的84至88行  第一次循环中此判断条件未成立，故还未检查此段代码执行的正确性
            count6++;
            // vector<double> tmpB = B;
            for(int i=0;i<B.size();i++) B(i) = this->interpCB(S,B,(1-this->pde_param_ptr->q)*S(i));
            for(int i=0;i<u.size();i++) u(i) = this->interpCB(S,u,(1-this->pde_param_ptr->q)*S(i));

        }

        for(int i=1;i<S.size()-1;i++)//循环赋值给alpha 在MATLAB源代码92行中是以向量化操作实现的
        { //alpha shape 2399 x 1
            double tmp1 = (this->pde_param_ptr->sigma*this->pde_param_ptr->sigma*S(i)*S(i))/(2*this->pde_param_ptr->ds*this->pde_param_ptr->ds);
            double tmp2 = (r+this->pde_param_ptr->p*this->pde_param_ptr->eta-this->pde_param_ptr->q)*S(i)/(2*this->pde_param_ptr->ds);
            double tmp3 = this->pde_param_ptr->dt;
            alpha(i-1) = ((tmp1-tmp2)*tmp3);
        }
        // cout<<alpha.transpose()<<endl;
        cout<<"alpha size = "<<alpha.size()<<endl;
        cout<<"alpha 计算完毕"<<endl;

        for(int i=1;i<S.size()-1;i++)//beta 在MATLAB源代码92行中是以向量化操作实现的
        { //beta shape 2399 x 1
            double tmp1 = (this->pde_param_ptr->sigma*this->pde_param_ptr->sigma*S(i)*S(i))/(2*this->pde_param_ptr->ds*this->pde_param_ptr->ds);
            double tmp2 = (r+this->pde_param_ptr->p*this->pde_param_ptr->eta-this->pde_param_ptr->q)*S(i)/(2*this->pde_param_ptr->ds);
            double tmp3 = this->pde_param_ptr->dt;
            beta(i-1) = ((tmp1+tmp2)*tmp3);
        }
        // cout<<beta.transpose()<<endl;
        cout<<"beta size = "<<beta.size()<<endl;
        cout<<"beta 计算完毕"<<endl;        

        for(int i=0;i<alpha.size();i++) alpha_MuMB(i) = alpha(i);
        alpha_MuMB(alpha_MuMB.size()-2) = 0;
        alpha_MuMB(alpha_MuMB.size()-1) = 0;
        // cout<<alpha_MuMB.transpose()<<endl;
        cout<<"alpha_MuMB size = "<<alpha_MuMB.size()<<endl;
        cout<<"alpha_MuMB 计算完毕"<<endl;        

        beta_MuMB(0) = 0;
        beta_MuMB(1) = 0;
        for(int i=0;i<beta.size();i++) beta_MuMB(i+2) = beta(i);
        // cout<<beta_MuMB.transpose()<<endl;
        cout<<"beta_MuMB size = "<<beta_MuMB.size()<<endl;
        cout<<"beta_MuMB 计算完毕"<<endl; 

        gamma_Mu(0) = -(r+this->pde_param_ptr->p)*this->pde_param_ptr->dt;
        for(int i=0;i<alpha.size();i++)
        {
            double tmp = -(alpha(i) + beta(i) + (r+this->pde_param_ptr->p)*this->pde_param_ptr->dt);
            gamma_Mu(i+1) = tmp;
        }
        // gamma_Mu.middleRows(1,2400) = -(alpha+beta+(r+this->pde_param_ptr->p*(1-this-> pde_param_ptr-> R))*this->pde_param_ptr->dt);
        // cout<<gamma_Mu.transpose()<<endl;
        cout<<"gamma_Mu size = "<<gamma_Mu.size()<<endl;
        cout<<"gamma_Mu 计算完毕"<<endl; 


        for(int i=0;i<height;i++)
        {
            Mu(i,i) = gamma_Mu(i);
        }
        for(int i=1;i<height;i++)
        {
            Mu(i-1,i) = beta_MuMB(i);
        }
        for(int i=1;i<height;i++)
        {
            Mu(i,i-1) = alpha_MuMB(i-1);
        }
        // cout<<Mu(1,1)<<endl;
        cout<<"Mu size = "<<Mu.size()<<endl;
        cout<<"Mu 计算完毕"<<endl; 

        gamma_MB(0) = -(r+this->pde_param_ptr->p*(1-this->pde_param_ptr->R))*this->pde_param_ptr->dt;
        for(int i=0;i<alpha.size();i++)
        {
            double tmp = -(alpha(i) + beta(i) + (r+this->pde_param_ptr->p*(1-this->pde_param_ptr->R))*this->pde_param_ptr->dt);
            gamma_MB(i+1) = tmp;
        }
        // gamma_MB.middleRows(1,2400) = -(alpha+beta+(r+this->pde_param_ptr->p*(1-this-> pde_param_ptr-> R))*this->pde_param_ptr->dt);
        // cout<<gamma_MB.transpose()<<endl;
        cout<<"gamma_MB size = "<<gamma_MB.size()<<endl;
        cout<<"gamma_MB 计算完毕"<<endl; 

        for(int i=0;i<height;i++)
        {
            MB(i,i) = gamma_MB(i);
        }
        for(int i=1;i<height;i++)
        {
            MB(i-1,i) = beta_MuMB(i);
        }
        for(int i=1;i<height;i++)
        {
            MB(i,i-1) = alpha_MuMB(i-1);
        }
        // cout<<MB(1,0)<<endl;
        cout<<"MB size = "<<MB.size()<<endl;
        cout<<"MB 计算完毕"<<endl; 

        AccI = (ceil(t) - t)*this->pde_param_ptr->coupon_time_rate(int(ceil(t)))*this->pde_param_ptr->F;
        cout<<"AccI = "<<AccI<<endl;

        if(this->pde_param_ptr->no_call_time(0)<=t&&t<=this->pde_param_ptr->no_call_time(1))
            Bc_n = 1000;
        else 
            Bc_n=this->pde_param_ptr->Bc_star+AccI;
        cout<<"Bc_n = "<<Bc_n<<endl;
        if(this->pde_param_ptr->no_put_time(0)<=t&&t<=this->pde_param_ptr->no_put_time(1))
            Bp_n=0.01;
        else
            Bp_n=this->pde_param_ptr->Bp_star+AccI;
        cout<<"Bp_n = "<<Bp_n<<endl;

        if(this->pde_param_ptr->no_convert_time(0)<=t&&t<=this->pde_param_ptr->no_convert_time(1))
        {
            for(int i=0;i<S.size();i++) 
            {   
                 k_n(i) = (0.01*1); 
            }
        }
        else
        {
            for(int i=0;i<S.size();i++)
            {
                double tmp1 = this->pde_param_ptr->F;
                double tmp2 = this->pde_param_ptr->conversion_price-this->pde_param_ptr->q*S(i)*(ceil(this->pde_param_ptr->T-this->pde_param_ptr->dtime)-ceil(t-this->pde_param_ptr->dtime));
                k_n(i) = (tmp1/tmp2);
            }
        }
        // cout<<"k_n = "<<k_n.transpose()<<endl;
        cout<<"k_n size = "<<k_n.size()<<endl;
        cout<<"k_n 计算完毕"<<endl; 

        u_n = u;
        // cout<<"u_n = "<<u_n.transpose()<<endl;
        cout<<"u_n size = "<<u_n.size()<<endl;
        cout<<"u_n 计算完毕"<<endl; 

        MatrixXd tmp1 = (MatrixXd::Identity(this->pde_param_ptr->Ns+1,this->pde_param_ptr->Ns+1) - (this->pde_param_ptr->theta)*MB);
        MatrixXd tmp2 = (MatrixXd::Identity(this->pde_param_ptr->Ns+1,this->pde_param_ptr->Ns+1) - (1-this->pde_param_ptr->theta)*MB)*B;
        // time_stt = clock();
        // MatrixXd tmp3 = tmp1.inverse(); //求逆矩阵非常耗时间
        cout<<"tmp1 size="<<tmp1.size()<<endl;
        cout<<"tmp2 size="<<tmp2.size()<<endl;
        // B = tmp3*tmp2;
        // cout << "time use in LU composition is       " << 1000 * (clock() - time_stt) / (double)
        //     CLOCKS_PER_SEC << "ms" << endl;
        time_stt = clock();
        B = tmp1.partialPivLu().solve(tmp2); //用LU分解替代直接求逆 方法参考 http://www.javashuo.com/article/p-fjkcsfmi-ey.html
        // cout<<B.transpose()<<endl;
        cout << "time use in LU composition is       " << 1000 * (clock() - time_stt) / (double)
            CLOCKS_PER_SEC << "ms" << endl;
        cout<<"B size = "<<B.size()<<endl;
        cout<<"B 计算完毕"<<endl;

        for(int i=0;i<B.size();i++)
        {
            if(B(i)>Bc_n) B(i) = Bc_n;
        }


        for(int i=0;i<S.size();i++)
        {
            double tmp1 = max(k_n(i)*(1-this->pde_param_ptr->eta)*S(i),this->pde_param_ptr->R*B(i));
            double tmp2 = this->pde_param_ptr->p*this->pde_param_ptr->dt;
            Muu(i) = (tmp1*tmp2);
        }
        // cout<<Muu.transpose()<<endl;
        cout<<"Muu size = "<<Muu.size()<<endl;
        cout<<"Muu 计算完毕"<<endl;

        tmp1 = (MatrixXd::Identity(this->pde_param_ptr->Ns+1,this->pde_param_ptr->Ns+1) - (this->pde_param_ptr->theta)*Mu);
        tmp2 = (MatrixXd::Identity(this->pde_param_ptr->Ns+1,this->pde_param_ptr->Ns+1) - (1-this->pde_param_ptr->theta)*Mu)*u_n;
        for(int i=0;i<this->pde_param_ptr->Ns;i++) tmp2(i) = tmp2(i) + Muu(i);
        tmp2(tmp2.size()-1) = tmp2(tmp2.size()-1) + 0;
        // tmp3 = tmp1.transpose();
        cout<<"tmp1 size="<<tmp1.size()<<endl;
        cout<<"tmp2 size="<<tmp2.size()<<endl;
        time_stt = clock();
        u = tmp1.partialPivLu().solve(tmp2); //用LU分解替代直接求逆
        cout << "time use in LU composition is       " << 1000 * (clock() - time_stt) / (double)
            CLOCKS_PER_SEC << "ms" << endl;
        // cout<<u.transpose()<<endl;
        cout<<"u size = "<<u.size()<<endl;
        cout<<"u 计算完毕"<<endl;

        // for(int i=0;i<u.size();i++)
        
        
        MatrixXd tmp111 = k_n.array()*S.array();
        MatrixXd tmp11 = tmp111.array().max(Bc_n);
        tmp1 = tmp11.array().min(u.array());
        MatrixXd tmp22 = k_n.array()*S.array();
        tmp2 = tmp22.array().max(Bp_n);
        u = tmp1.array().max(tmp2.array());
        // cout<<"u -> "<<u.transpose()<<endl;
        for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
        {
            double tmp;
            if(u(i)<=tmp2(i))
                tmp = 1;
            else
                tmp = 0;
            P1_old(i,i) = tmp;
        }
        P1_old =  this->pde_param_ptr->rhopenltyput * P1_old;
        // cout<<"P1_old \n"<<P1_old.transpose()<<endl;
        for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
        {
            double tmp;
            if(u(i)>=tmp11(i))
                tmp = 1;
            else
                tmp = 0;
            P2_old(i,i) = tmp;
        }
        P2_old = -this->pde_param_ptr->rhopenltycall * P2_old;
        // cout<<"P2_old \n"<<P2_old.transpose()<<endl;
        while(1)
        {
            count++;
            u_old = u;
            cout<<"u_old size = "<<u_old.size();
            for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
            {
                double tmp;
                if(u(i)<=tmp2(i))
                    tmp = 1;
                else
                    tmp = 0;
                P1(i,i) = tmp;
            }
            cout<<"complete"<<endl;
            P1 = this->pde_param_ptr->rhopenltyput * P1;
            for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
            {
                double tmp;
                if(u(i)>=tmp11(i))
                    tmp = 1;
                else
                    tmp = 0;
                P2_old(i,i) = tmp;
            }
            P2 = -this->pde_param_ptr->rhopenltycall * P2;
            MatrixXd tmpu1(height,weight); 
            // cout<<"speye size = "<<(MatrixXd::Identity(this->pde_param_ptr->Ns+1,this->pde_param_ptr->Ns+1)).size()<<endl;
            // cout<<"Mu size = "<<Mu.size()<<endl;
            // cout<<"P1 size = "<<P1.size()<<endl;
            // cout<<"P2 size = "<<P2.size()<<endl;
            tmpu1 = (MatrixXd::Identity(this->pde_param_ptr->Ns+1,this->pde_param_ptr->Ns+1) - (this->pde_param_ptr->theta)*Mu + P1- P2 );
            // cout<<tmpu1<<endl;

            MatrixXd tmpu2(height,weight); 
            tmpu2 = (MatrixXd::Identity(this->pde_param_ptr->Ns+1,this->pde_param_ptr->Ns+1) - (1-this->pde_param_ptr->theta)*Mu)*u_n;
            for(int i=0;i<this->pde_param_ptr->Ns;i++) tmpu2(i) = tmpu2(i) + Muu(i);
            tmpu2(tmpu2.size()-1) = tmpu2(tmpu2.size()-1) + 0;
            
            MatrixXd tmpu3(height,weight);
            cout<<"P1 size = "<<P1.size()<<endl;
            tmpu3 = (k_n.array()*S.array()).array().max(Bp_n);
            for(int i=0;i<height;i++)
            {
                if(P1(i,i)==1)
                    tmpu3(i) *= 1;
                else
                    tmpu3(i) *= 0;
            }
            cout<<"tmpu3 size = "<<tmpu3.size()<<endl;
            // cout<<"$$$-> "<<tmpu3.transpose()*tmpu3<<endl;

            MatrixXd tmpu4(height,weight);
            tmpu4 = (k_n.array()*S.array()).array().min(Bc_n);
            for(int i=0;i<height;i++)
            {
                if(P2(i,i)==1)
                    tmpu4(i) *= 1;
                else
                    tmpu4(i) *= 0;
            }

            u = tmpu1.partialPivLu().solve(tmpu2) + tmpu3 - tmpu4 ;
            // cout<<"u -> \n"<<u.transpose()<<endl;
            cout<<"u size = "<<u.size()<<endl;
            for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
            {
                double tmp;
                if(u(i)<=tmp2(i))
                    tmp = 1;
                else
                    tmp = 0;
                P1(i,i) = tmp;
            }
            cout<<"P1 size = "<<P1.size()<<endl;
            P1 = this->pde_param_ptr->rhopenltyput * P1;
            for(int i=0;i<this->pde_param_ptr->Ns+1;i++)
            {
                double tmp;
                if(u(i)>=tmp11(i))
                    tmp = 1;
                else
                    tmp = 0;
                P2_old(i,i) = tmp;
            }
            P2 = -this->pde_param_ptr->rhopenltycall * P2;
            cout<<"P2 size = "<<P2.size()<<endl;

            

            bool sameP1=1;
            for(int i=0;i<height;i++)
            {
                for(int j =0;j<weight;j++)
                {
                    if (P1(i,j)!=P1_old(i,j))
                    {
                        sameP1 = 0;
                        break;
                    }
                }
            }
            bool sameP2=1;
            for(int i=0;i<height;i++)
            {
                for(int j =0;j<weight;j++)
                {
                    if (P2(i,j)!=P2_old(i,j))
                    {   
                        sameP2 = 0;
                        break;
                    }
                }
            }

            bool u_flag = 1;
            MatrixXd tmpuflag = (u-u_old).array().abs().array()/u.array().abs().array().max(1).array();
            double maxval = tmpuflag(0);
            for(int i=0;i<tmpuflag.size();i++) 
            {
                if(tmpuflag(i)>maxval) maxval = tmpuflag(i);
            }
            if (maxval < 1/0.0001) u_flag = 0;
            cout<<"tmpuflag size = "<<tmpuflag.size()<<endl;
            if(sameP1&&sameP2) 
            {
                count1++;  
                break;
            }
            else if (u_flag)
            {
                count2++;
                break;
            }
            else if(count > this->pde_param_ptr->Nt+40)
            {
                count4++;
                cout<<"Probably does not converge, current n is "<<n<<endl;
                flag = 1;
                break;
            }
            else
            { 
                count3++;
            }
            P1_old = P1;
            P2_old = P2;
            

        }
        cout<<"while loop OVER"<<endl;
    }

    double Uindex = this->pde_param_ptr->Ns*this->pde_param_ptr->S0/S(S.size()-1)+1;
    cout<<"Uindex = "<<Uindex<<endl;
    double U;
    if(Uindex >1)
        U = u(int(ceil(Uindex)))*(Uindex-floor(Uindex+0.001))+u(int(floor(Uindex)))*(ceil(Uindex+0.001)-Uindex);
    else 
        U = u(int(ceil(Uindex)));
    return U;
}