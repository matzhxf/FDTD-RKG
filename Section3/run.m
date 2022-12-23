clc
clear

num=1;
plot_type = 2;  % plot_mod=其它,不画图; =1,质量能量图; =2,动态图; =3,最大模范数变化图
line_type = '-'; 

[x_left,x_right,y_left,y_right,t_begin,t_end, ...,
    phi_0,phi_1,gamma,lambda,h1,h2,tau,J,K,N] = Variable_setting(num);
%======== 计算psi_N ========
tic;
% [psi_N] = generate_exact_scheme(num,plot_type,line_type,'LCFD','\mathcal{L}_z^{2,I}');
% [psi_N] = BS4_scheme(num,plot_type,line_type,'LCFD','\mathcal{L}_z^{2,I}');
% [psi_N] = BS4_scheme(num,plot_type,line_type,'NCFD-II','\mathcal{L}_z^{2,I}');
% [psi_N] = BS4_scheme(num,plot_type,line_type,'RCFD','\mathcal{L}_z^{2,I}');

% [psi_N] = BS2_scheme(num,plot_type,line_type,'LCFD','L_z^{2,I}');
% [psi_N] = BS2_scheme(num,plot_type,line_type,'NCFD-I','L_z^{2,I}');
% [psi_N] = BS2_scheme(num,plot_type,line_type,'NCFD-II','L_z^{2,I}');
% [psi_N] = BS2_scheme(num,plot_type,line_type,'RCFD','L_z^{2,I}');

[psi_N] = BS2_scheme(num,plot_type,line_type,'LCFD','\hat{L}_z^{2,h}');
% [psi_N] = BS2_scheme(num,plot_type,line_type,'NCFD-I','\hat{L}_z^{2,h}');
% [psi_N] = BS2_scheme(num,plot_type,line_type,'NCFD-II','\hat{L}_z^{2,h}');
% [psi_N] = BS2_scheme(num,plot_type,line_type,'RCFD','\hat{L}_z^{2,h}');
toc;
% ======== 储存数值精确解 ========
% exact_solution1 = psi_N;
% name = 'exact_solution1'; 
% save(name,name);