clc;
clear;

plot_type = 3;  % plot_mod=其它,不画图; =1,质量能量图; =2,动态图; =3,最大模范数变化图

Picture_mod={'-k',':.b','--r','-.g'};

figure;
for num=1:4
    hold on;
    [psi_N]=BS2_scheme(num,plot_type,Picture_mod{num},'LCFD','\hat{L}_z^{2,h}');
    box on;
end  
legend('$\gamma=0.2,\tau=h/4=1/64$','$\gamma=0.3,\tau=h/4=1/64$', ...,
    '$\gamma=0.4,\tau=h/4=1/64$','$\gamma=0.4,\tau=h^2=1/256$','Interpreter','Latex','FontSize',14);
saveas(gcf,'C:\Users\Mr.Wang\Desktop\RKG-FD\matlab code\Picture_fig\Section3\BUS2\LCFD_infty_norm');

figure;
for num=1:4
    hold on;
    [psi_N]=BS2_scheme(num,plot_type,Picture_mod{num},'NCFD-I','\hat{L}_z^{2,h}');
    box on;
end  
legend('$\gamma=0.2,\tau=h/4=1/64$','$\gamma=0.3,\tau=h/4=1/64$', ...,
    '$\gamma=0.4,\tau=h/4=1/64$','$\gamma=0.4,\tau=h^2=1/256$','Interpreter','Latex','FontSize',14);
saveas(gcf,'C:\Users\Mr.Wang\Desktop\RKG-FD\matlab code\Picture_fig\Section3\BUS2\NCFD-I_infty_norm');

figure;
for num=1:4
    hold on;
    [psi_N]=BS2_scheme(num,plot_type,Picture_mod{num},'RCFD','\hat{L}_z^{2,h}');
    box on;
end  
legend('$\gamma=0.2,\tau=h/4=1/64$','$\gamma=0.3,\tau=h/4=1/64$', ...,
    '$\gamma=0.4,\tau=h/4=1/64$','$\gamma=0.4,\tau=h^2=1/256$','Interpreter','Latex','FontSize',14);
saveas(gcf,'C:\Users\Mr.Wang\Desktop\RKG-FD\matlab code\Picture_fig\Section3\BUS2\RCFD_infty_norm');

close all;