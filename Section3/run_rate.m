clc
clear
%============== RNKG方程变量预设置 ===============
plot_type = 0;  % plot_mod=其它,不画图; =1,质量能量图; =2,动态图; =3,最大模范数变化图
Scheme_type = {'LCFD','NCFD-I'};
%======== 计算psi_N ========
load exact_solution3.mat;
exact_solution = exact_solution3; 
err_norm = zeros(1,4);
rate = zeros(1,4);
exact_solution_h = 1/64;

% disp('=========  BS2  =========');
% %% 空间收敛阶
% for scheme_num = 2:2
%     for num=1:4
%         %======== 变量声明 ========
%         [x_left,x_right,y_left,y_right,t_begin,t_end, ...,
%             phi_0,phi_1,gamma,lambda,h1,h2,tau,J,K,N] = Variable_setting(num);
% 
%         [psi_N] = BS2_scheme(num,plot_type,'-', ...,
%             Scheme_type{scheme_num},'L_z^{2,I}');
% 
%         err_point = zeros(1,J-1);
%         for j=1:J-1
%             exact_solution_j=floor(j.*h1./exact_solution_h);
%             theta_x=(j.*h1-exact_solution_j.*exact_solution_h)./exact_solution_h;
%             for k=1:K-1
%                 exact_solution_k=floor(k.*h2./exact_solution_h);
%                 theta_y=(k.*h2-exact_solution_k.*exact_solution_h)./exact_solution_h;
%                 err_point(j,k)=abs(psi_N(j+1,k+1)-(1-theta_x).*(1-theta_y) ...,
%                     .*exact_solution(exact_solution_j+1, ...,
%                     exact_solution_k+1)-(1-theta_x).*theta_y.* ...,
%                     exact_solution(exact_solution_j+1, ...,
%                     exact_solution_k+2)-(1-theta_y).*theta_x.* ...,
%                     exact_solution(exact_solution_j+2, ...,
%                     exact_solution_k+1)-theta_x.*theta_y.* ...,
%                     exact_solution(exact_solution_j+2, ...,
%                     exact_solution_k+2));
%             end
%         end
%         err_norm(1,num)=max(max(err_point));
%         if num>1
%             rate(num)=log(err_norm(1,num-1)./err_norm(1,num))./log(2);
%         end
%     end
%     disp([Scheme_type{scheme_num},'——spatial rate']);
%     disp(['err_norm= &',num2str(err_norm(1,1),' %.4E'), ...,
%         ' &',num2str(err_norm(1,2),' %.4E'), ...,
%         ' &',num2str(err_norm(1,3),' %.4E'), ...,
%         ' &',num2str(err_norm(1,4),' %.4E')]);
%     disp(['rate=     &--', ...,
%         ' &',num2str(rate(1,2),' %.2f'), ...,
%         ' &',num2str(rate(1,3),' %.2f'), ...,
%         ' &',num2str(rate(1,4),' %.2f')]);
%     disp('==================');
% end

disp('=========  BS4-\mathcal{L}_z^{2,I}  =========');
%% 空间收敛阶
for scheme_num = 1:1
    for num=1:1
        %======== 变量声明 ========
        [x_left,x_right,y_left,y_right,t_begin,t_end, ...,
            phi_0,phi_1,gamma,lambda,h1,h2,tau,J,K,N] = Variable_setting(num);

        [psi_N] = BS4_scheme(num,plot_type,'-', ...,
            Scheme_type{scheme_num},'\mathcal{L}_z^{2,I}');

        err_point = zeros(1,J-1);
        for j=1:J-1
            exact_solution_j=floor(j.*h1./exact_solution_h);
            theta_x=(j.*h1-exact_solution_j.*exact_solution_h)./exact_solution_h;
            for k=1:K-1
                exact_solution_k=floor(k.*h2./exact_solution_h);
                theta_y=(k.*h2-exact_solution_k.*exact_solution_h)./exact_solution_h;
                err_point(j,k)=abs(psi_N(j+1,k+1)-(1-theta_x).*(1-theta_y) ...,
                    .*exact_solution(exact_solution_j+1, ...,
                    exact_solution_k+1)-(1-theta_x).*theta_y.* ...,
                    exact_solution(exact_solution_j+1, ...,
                    exact_solution_k+2)-(1-theta_y).*theta_x.* ...,
                    exact_solution(exact_solution_j+2, ...,
                    exact_solution_k+1)-theta_x.*theta_y.* ...,
                    exact_solution(exact_solution_j+2, ...,
                    exact_solution_k+2));
            end
        end
        err_norm(1,num)=max(max(err_point));
        if num>1
            rate(num)=log(err_norm(1,num-1)./err_norm(1,num))./log(2);
        end
    end
    disp([Scheme_type{scheme_num},'——spatial rate']);
    disp(['err_norm= &',num2str(err_norm(1,1),' %.4E'), ...,
        ' &',num2str(err_norm(1,2),' %.4E'), ...,
        ' &',num2str(err_norm(1,3),' %.4E'), ...,
        ' &',num2str(err_norm(1,4),' %.4E')]);
    disp(['rate=     &--', ...,
        ' &',num2str(rate(1,2),' %.2f'), ...,
        ' &',num2str(rate(1,3),' %.2f'), ...,
        ' &',num2str(rate(1,4),' %.2f')]);
    disp('==================');
end


