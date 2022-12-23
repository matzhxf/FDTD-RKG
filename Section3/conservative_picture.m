clc;
clear;

plot_type = 1;  % plot_mod=其它,不画图; =1,质量能量图; =2,动态图; =3,最大模范数变化图
Scheme_type = {'LCFD','NCFD-I'};

fig = 1;
for num=1:1
    [psi_N]=BS4_scheme(num,plot_type,'-',Scheme_type{num},'\mathcal{L}_z^{2,I}');
    saveas(fig,['C:\Users\Mr.Wang\Desktop\RKG-FD\matlab code\Picture_fig\Section3\BS4\', ...,
        Scheme_type{num},'_mass']);
    saveas(fig+1,['C:\Users\Mr.Wang\Desktop\RKG-FD\matlab code\Picture_fig\Section3\BS4\', ...,
        Scheme_type{num},'_energy']);
    fig = fig+2;
end
close all;
