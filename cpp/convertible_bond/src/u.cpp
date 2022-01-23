u=max(min(u,max(Bc_n,(k_n.*S+((ceil(t+dt)-ceil(t))==1).*coupon_time_rate(ceil(t+dt+1e-6))*F))),...
        max(Bp_n,(k_n.*S+((ceil(t+dt)-ceil(t))==1).*coupon_time_rate(ceil(t+dt+1e-6))*F)));


tmp_1 = max(Bp_n,(k_n[i]*this->pde_param_ptr->S+((ceil(this->pde_param_ptr->t+this->pde_param_ptr->dt)-ceil(this->pde_param_ptr->t))==1)*this->pde_param_ptr->coupon_time_rate[ceil(this->pde_param_ptr->t+this->pde_param_ptr->dt+1e-6)]*this->pde_param_ptr->F));
max(Bp_n,(k_n.*S+((ceil(t+dt)-ceil(t))==1).*coupon_time_rate(ceil(t+dt+1e-6))*F))

tmp_2 = min(u[i],max(Bc_n,(k_n[i]*this->pde_param_ptr->S+((ceil(this->pde_param_ptr->t+this->pde_param_ptr->dt)-ceil(this->pde_param_ptr->t))==1)*this->pde_param_ptr->coupon_time_rate[ceil(this->pde_param_ptr->t+this->pde_param_ptr->dt+1e-6)]*this->pde_param_ptr->F)));
min(u,max(Bc_n,(k_n.*S+((ceil(t+dt)-ceil(t))==1).*coupon_time_rate(ceil(t+dt+1e-6))*F)))

tmp_3 = max(tmp_2,tmp_1);
u[i] = tmp_3;