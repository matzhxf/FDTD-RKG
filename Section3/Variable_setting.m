function [x_left,x_right,y_left,y_right,t_begin,t_end, ...,
    phi_0,phi_1,gamma,lambda,h1,h2,tau,J,K,N] = Variable_setting(num)
%VARIABLE_SETTING RKG方程变量预设置
%   exact_solution: T=1.5, gamma=0.2, h=1/64, tau=1/1000,  BS4
%   exact_solution1: T=1, gamma=0.2, h=1/64, tau=1/4000,  BS4 
%   exact_solution3: T=1, gamma=1, h=1/64, tau=1/4000,  BS4 

[x_left,x_right]=deal(-4,4);  
[y_left,y_right]=deal(-4,4);
[t_begin,t_end]=deal(0,2); %% 

phi_0=@(x,y) ( exp(-2.*(x.^2+y.^2)) ); 
phi_1=@(x,y) ( (1+1i).*exp(-2.*(x.^2+y.^2)) );

Gamma=[0.2,0.3,0.4,0.4]; %%
Space_step=[1/16,1/16,1/16,1/16]; %%
Time_step=[1/64,1/64,1/64,1/256]; %% 

[lambda,gamma]=deal(1,Gamma(num));
[h1,h2,tau]=deal(Space_step(num),Space_step(num),Time_step(num));
[J,K,N]=deal(ceil((x_right-x_left)/h1),ceil((y_right-y_left)/h2), ...,
    ceil((t_end-t_begin)/tau));
end

