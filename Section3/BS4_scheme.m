function [psi_N] = BS4_scheme(num,plot_type,line_type,scheme_type, ...,
    Quadratic_rotation_term)
%FD_SCHEME 此处显示有关此函数的摘要
%   此处显示详细说明
%% 报错
if ismember(scheme_type,{'LCFD','NCFD-I','NCFD-II','RCFD'}) == 0
    error('scheme_type 输入错误');
elseif ismember(Quadratic_rotation_term,{'\mathcal{L}_z^{2,I}', ...,
        '\mathcal{L}_z^{2,II}','\mathcal{L}_z^{2,III}'}) == 0
    error('Quadratic_rotation_term 输入错误');
end
%% 变量准备
tol = 1e-10;
iteration_max = 100;

[x_left,x_right,y_left,y_right,t_begin,t_end, ...,
    phi_0,phi_1,gamma,lambda,h1,h2,tau,J,K,N] = Variable_setting(num);
[scheme_point] = point_type_judgment(Quadratic_rotation_term);

xh_inner = x_left+h1 : h1 : x_left+(J-1)*h1;
yh_inner = y_left+h2 : h2 : y_left+(K-1)*h2;
t_all_point = t_begin : tau : t_end;
[Yh_Inner,Xh_Inner] = meshgrid(yh_inner,xh_inner);
[r1,r2] = deal(tau./h1.^2,tau./h2.^2);

[xh,yh] = deal(xh_inner,yh_inner);
for m = 1 : (scheme_point-1)/2
    xh=[xh(1)-h1,xh,xh(end)+h1];
    yh=[yh(1)-h2,yh,yh(end)+h2];
end
[Yh,Xh] = meshgrid(yh,xh);
[j1,jend,k1,kend] = deal(m+1,m+(J-1),m+1,m+(K-1));
[K_Index,J_Index] = meshgrid((1-m):1:K+(m-1),(1-m):1:J+(m-1));
[number_of_x_points,number_of_y_points] = size(Xh);

%% 生成稀疏矩阵
if ismember(scheme_type,{'LCFD','NCFD-I'}) == 1
    A = 2.*generate_A(1)-generate_A(2)+generate_A(3) ...,
        -2.*gamma.*generate_A(4)-tau.*gamma.^2.*generate_A(5);
elseif ismember(scheme_type,{'NCFD-II','RCFD'}) == 1
    A = 4.*generate_A(1)-generate_A(2)+generate_A(3) ...,
        -4.*gamma.*generate_A(4)-tau.*gamma.^2.*generate_A(5);
end

%% 数值解计算
% 数值解分层预设
psin_minus = zeros(number_of_x_points,number_of_y_points);
psin = zeros(number_of_x_points,number_of_y_points);
phi1 = zeros(number_of_x_points,number_of_y_points);

Infinity_norm = zeros(1,N+1);
if plot_type == 1   %=========== 守恒律图
    if strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,I}') == 1 
        if ismember(scheme_type,{'LCFD','NCFD-I','NCFD-II'}) == 1
            [Mh,Error_Mh,Eh,Error_Eh] = deal(zeros(1,N),zeros(1,N), ...,
                zeros(1,N),zeros(1,N));
        elseif ismember(scheme_type,{'RCFD'}) == 1
            [Mh,Error_Mh,Eh,Error_Eh] = deal(zeros(1,N+1),zeros(1,N+1), ...,
                zeros(1,N+1),zeros(1,N+1));
        end
    else 
        disp('该方法没有守恒律');
    end
end


% n=1时
psin_minus(j1:jend,k1:kend) = feval(phi_0,Xh_Inner,Yh_Inner); % 载入第0层
Infinity_norm(1,1) = max(max(abs(psin_minus)));
phi1(j1:jend,k1:kend) = feval(phi_1,Xh_Inner,Yh_Inner);

% 计算第1层(此时n=1)
if ismember(scheme_type,{'LCFD','NCFD-I','NCFD-II'}) == 1
    psin(j1:jend,k1:kend) = ...,
        psin_minus(j1:jend,k1:kend)+tau.*phi1(j1:jend,k1:kend) ...,
        +0.5.*tau.^2.*(Delta_h(psin_minus)-(1+lambda ...,
        .*abs(psin_minus(j1:jend,k1:kend)).^2).*psin_minus(j1:jend,k1:kend) ...,
        +2i.*gamma.*L_zh(phi1)+gamma.^2.*L_zh_2(psin_minus));   
elseif ismember(scheme_type,{'RCFD'}) == 1
    uh_s = zeros(number_of_x_points,number_of_y_points,2);
    R = zeros(number_of_x_points,number_of_y_points);
    uh_s(j1:jend,k1:kend,1) = psin_minus(j1:jend,k1:kend); % 载入迭代初值
    for s = 1:iteration_max
        R(j1:jend,k1:kend) =  ...,
            4./tau.*psin_minus(j1:jend,k1:kend)+4.*phi1(j1:jend,k1:kend) ...,
            +tau.*Delta_h(psin_minus)-tau.*(1 ...,
            +0.5.*lambda.*(abs(psin_minus(j1:jend,k1:kend)).^2 ...,
            +abs(uh_s(j1:jend,k1:kend,1)).^2)).*(psin_minus(j1:jend,k1:kend)) ...,
            -4i.*gamma.*L_zh(psin_minus)+tau.*gamma.^2.*L_zh_2(psin_minus);
        R_vector = reshape(R,[],1);
        A1 = A;
        for j = 1:J-1
            for k = 1:K-1
                A1(Index(j,k),Index(j,k)) = A(Index(j,k),Index(j,k)) ...,
                    +0.5.*tau.*lambda.*(abs(psin_minus(j1+j-1,k1+k-1)).^2 ...,
                    +abs(uh_s(j1+j-1,k1+k-1,1)).^2);
            end
        end
        uh_s(:,:,2) = reshape(A1\R_vector,number_of_x_points,number_of_y_points);
        err = max(max(abs(uh_s(:,:,2)-uh_s(:,:,1))));
        uh_s(:,:,1) = uh_s(:,:,2);
        if err < tol
            break
        end
        if s == iteration_max
            warning_text = '计算第1层时迭代次数过多';
            warning(warning_text);
        end
    end
    psin = uh_s(:,:,2);
end
Infinity_norm(1,2) = max(max(abs(psin)));

if plot_type == 1 &&  ...,    %=========== 守恒律图
        strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,I}') == 1       
    if ismember(scheme_type,{'LCFD','NCFD-I','NCFD-II'}) == 1
        [Mh(1,1),Eh(1,1)] = Comput_Mass_Energy(psin_minus,psin);
    elseif ismember(scheme_type,{'RCFD'}) == 1
        u = phi1;
        [Mh(1,1),Eh(1,1)] = Comput_Mass_Energy(psin,psin_minus);
        u = 2./tau.*(psin-psin_minus)-u;
        [Mh(1,2),Eh(1,2)] = Comput_Mass_Energy(psin_minus,psin);
        Error_Mh(1,2) = abs(Mh(1,2)-Mh(1,1))./Mh(1,1);
        Error_Eh(1,2) = abs(Eh(1,2)-Eh(1,1))./Eh(1,1);
    end
elseif plot_type == 2
    figure;
    Dynal=moviein(41);   %=========== 动态图
    split_N=ceil(N./40);
    split_n=1;
    surf(Xh,Yh,abs(psin_minus));
    x_tick=x_left:2:x_right; y_tick=y_left:2:y_right;    
    axis([x_left x_right y_left y_right 0 1.5]);
    set(gca,'XTick',x_tick,'YTick',y_tick);
    set(gca,'Fontsize',18,'position',[0.12,0.17,0.82,0.75]);
    xlabel('\itx','Fontname','Times New Roman','FontSize',20); 
    ylabel('\ity','Fontname','Times New Roman','FontSize',20);
    colorbar; shading interp;  view(0,90);
    title('$t=0$','Interpreter','Latex','FontSize',20);
    Dynal(:,1)=getframe;

end

% 预设循环变量
uh_s = zeros(number_of_x_points,number_of_y_points,2);
R = zeros(number_of_x_points,number_of_y_points);

% 计算第n+1层(n>=1) ========================
for n=1:N-1
    uh_s(j1:jend,k1:kend,1) =  ...,
        2.*psin(j1:jend,k1:kend)-psin_minus(j1:jend,k1:kend); % 载入迭代初值
    for s = 1:iteration_max
        R(j1:jend,k1:kend) = computer_R(scheme_type);
        R_vector = reshape(R,[],1);
        A1 = A+generate_A(6);
        uh_s(:,:,2) = reshape(A1\R_vector,number_of_x_points,number_of_y_points);
        err = max(max(abs(uh_s(:,:,2)-uh_s(:,:,1))));
        uh_s(:,:,1) = uh_s(:,:,2);
        if err < tol
            break
        end
        if s == iteration_max
            warning_text = ['计算第',num2str(n+1),'层时迭代次数过多'];
            warning(warning_text);
        end
    end
    psin_add = uh_s(:,:,2);
    Infinity_norm(1,n+2) = max(max(abs(psin_add)));

    if Infinity_norm(1,n+2)>10
        warning('该方法产生了数值爆破');
        Infinity_norm(1,n+2:end)=Inf;
        break;
    end

    psin_minus = psin;
    psin = psin_add;

    if plot_type == 1 &&  ...,    %=========== 守恒律图
            strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,I}') == 1
        if ismember(scheme_type,{'LCFD','NCFD-I','NCFD-II'}) == 1
            [Mh(1,n+1),Eh(1,n+1)] = Comput_Mass_Energy(psin_minus,psin);
            Error_Mh(1,n+1) = abs(Mh(1,n+1)-Mh(1,1))./Mh(1,1);
            Error_Eh(1,n+1) = abs(Eh(1,n+1)-Eh(1,1))./Eh(1,1);
        elseif ismember(scheme_type,{'RCFD'}) == 1
            u = 2./tau.*(psin-psin_minus)-u;
            [Mh(1,n+2),Eh(1,n+2)] = Comput_Mass_Energy(psin_minus,psin);
            Error_Mh(1,n+2) = abs(Mh(1,n+2)-Mh(1,1))./Mh(1,1);
            Error_Eh(1,n+2) = abs(Eh(1,n+2)-Eh(1,1))./Eh(1,1);
        end
    elseif plot_type == 2 && mod(n+1,split_N)==0  %===========画动态图
        split_n=split_n+1;
        surf(Xh,Yh,abs(psin));
        x_tick=x_left:2:x_right; y_tick=y_left:2:y_right;
        axis([x_left x_right y_left y_right 0 1.5]);
        set(gca,'XTick',x_tick,'YTick',y_tick);
        set(gca,'Fontsize',18,'position',[0.12,0.17,0.82,0.75]);
        xlabel('\itx','Fontname','Times New Roman','FontSize',20);
        ylabel('\ity','Fontname','Times New Roman','FontSize',20);
        colorbar; shading interp;  view(0,90);
        title(['$t=',num2str((n+1).*tau),'$'], ...,
            'Interpreter','Latex','FontSize',20);
        Dynal(:,split_n)=getframe;
    end
end
psi_N=psin(j1-1:jend+1,k1-1:kend+1);

%% 画图 =================================================================================
if plot_type == 3   %=========== 最大模范数变化图
    plot(t_all_point,Infinity_norm,line_type,'LineWidth',2);
    % 坐标轴刻度设置
    axis([0 t_end 0 5]);
    y_tick=0:1:5;
    axis([t_begin t_end y_tick(1) y_tick(end)]);
    set(gca,'YTick',y_tick);
    set(gca,'Fontsize',18,'position',[0.13,0.17,0.8,0.73]);
    xlabel('$t$','Interpreter','Latex','FontSize',20); 
    ylabel('$\|\psi^{n}\|_{L^{\infty}}$','Interpreter','Latex','FontSize',20);
    title(scheme_type,'FontSize',20);

elseif plot_type == 1 &&  ...,    %=========== 守恒律图
        strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,I}') == 1
    if ismember(scheme_type,{'LCFD','NCFD-I','NCFD-II'}) == 1
        t_all_point(1)=[];
    end
    % 质量误差图像
    figure
    plot(t_all_point,Error_Mh);
    set(gca,'Fontsize',18,'position',[0.15,0.17,0.8,0.75]);
    xlabel('\itt','Fontname','Times New Roman','FontSize',20);
    ylabel('$|Q^{n}-Q^{0}|/Q^{0}$','Interpreter','Latex','FontSize',20);
    title(['BS4-',scheme_type,' with $',Quadratic_rotation_term,'$'],'Interpreter','Latex','FontSize',20);

    % 能量误差图像
    figure
    plot(t_all_point,Error_Eh);
    set(gca,'Fontsize',18,'position',[0.15,0.17,0.8,0.75]);
    xlabel('\itt','Fontname','Times New Roman','FontSize',20);
    ylabel('$|E^{n}-E^{0}|/E^{0}$','Interpreter','Latex','FontSize',20);
    title(['BS4-',scheme_type,' with $',Quadratic_rotation_term,'$'],'Interpreter','Latex','FontSize',20);
end



%% 格式点数类型判断函数 =================================================================================
    function [method_point] = point_type_judgment(coe1)
        if strcmp(coe1,'\mathcal{L}_z^{2,I}') == 1 
            method_point = 9;
        elseif strcmp(coe1,'\mathcal{L}_z^{2,II}') == 1
            method_point = 5;
        elseif strcmp(coe1,'\mathcal{L}_z^{2,III}') == 1
            method_point = 9;
        end
    end

%% 计算Delta_h =================================================================================
    function [result] = Delta_h(uh)   
        result=1./(12.*h1.^2).*( -1.*uh(j1-2:jend-2,k1:kend) ...,
            +16.*uh(j1-1:jend-1,k1:kend)-30.*uh(j1:jend,k1:kend) ...,
            +16.*uh(j1+1:jend+1,k1:kend)-uh(j1+2:jend+2,k1:kend) ) ...,
             ...,
            +1./(12.*h2.^2).*( -1.*uh(j1:jend,k1-2:kend-2) ...,
            +16.*uh(j1:jend,k1-1:kend-1)-30.*uh(j1:jend,k1:kend) ...,
            +16.*uh(j1:jend,k1+1:kend+1)-uh(j1:jend,k1+2:kend+2) );
    end

%% 计算L_zh =================================================================================
    function [result] = L_zh(uh)
        result=Xh_Inner./(12.*h2).*( uh(j1:jend,k1-2:kend-2) ...,
            -8.*uh(j1:jend,k1-1:kend-1)+8.*uh(j1:jend,k1+1:kend+1) ...,
            -uh(j1:jend,k1+2:kend+2) ) ...,
             ...,
            -Yh_Inner./(12.*h1).*( uh(j1-2:jend-2,k1:kend) ...,
            -8.*uh(j1-1:jend-1,k1:kend)+8.*uh(j1+1:jend+1,k1:kend) ...,
            -uh(j1+2:jend+2,k1:kend) ) ;
        result=-1i.*result;
    end

%% 计算L_zh_2 =================================================================================
    function [result] = L_zh_2(uh)
        if strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,I}') == 1 
            result=(-1./(144.*h2.^2)).*Xh_Inner.^2.*( uh(j1:jend,k1-4:kend-4) ...,
                -16.*uh(j1:jend,k1-3:kend-3)+64.*uh(j1:jend,k1-2:kend-2) ...，
                +16.*uh(j1:jend,k1-1:kend-1)-130.*uh(j1:jend,k1:kend) ...，
                +16.*uh(j1:jend,k1+1:kend+1)+64.*uh(j1:jend,k1+2:kend+2) ...,
                -16.*uh(j1:jend,k1+3:kend+3)+uh(j1:jend,k1+4:kend+4) ) ...,
                ...,
                +Yh_Inner./(144.*h1.*h2).*( Xh(j1-2:jend-2,k1:kend).* ...,
                (uh(j1-2:jend-2,k1-2:kend-2)-8.*uh(j1-2:jend-2,k1-1:kend-1) ...,
                +8.*uh(j1-2:jend-2,k1+1:kend+1)-uh(j1-2:jend-2,k1+2:kend+2)) ...,
                -8.*Xh(j1-1:jend-1,k1:kend).* ...,
                (uh(j1-1:jend-1,k1-2:kend-2)-8.*uh(j1-1:jend-1,k1-1:kend-1) ...,
                +8.*uh(j1-1:jend-1,k1+1:kend+1)-uh(j1-1:jend-1,k1+2:kend+2)) ...,
                +8.*Xh(j1+1:jend+1,k1:kend).* ...,
                (uh(j1+1:jend+1,k1-2:kend-2)-8.*uh(j1+1:jend+1,k1-1:kend-1) ...,
                +8.*uh(j1+1:jend+1,k1+1:kend+1)-uh(j1+1:jend+1,k1+2:kend+2)) ...,
                -Xh(j1+2:jend+2,k1:kend).* ...,
                (uh(j1+2:jend+2,k1-2:kend-2)-8.*uh(j1+2:jend+2,k1-1:kend-1) ...,
                +8.*uh(j1+2:jend+2,k1+1:kend+1)-uh(j1+2:jend+2,k1+2:kend+2)) )...,
                ...,
                +Xh_Inner./(144.*h1.*h2).*( Yh(j1:jend,k1-2:kend-2).* ...,
                (uh(j1-2:jend-2,k1-2:kend-2)-8.*uh(j1-1:jend-1,k1-2:kend-2) ...,
                +8.*uh(j1+1:jend+1,k1-2:kend-2)-uh(j1+2:jend+2,k1-2:kend-2)) ...,
                -8.*Yh(j1:jend,k1-1:kend-1).* ...,
                (uh(j1-2:jend-2,k1-1:kend-1)-8.*uh(j1-1:jend-1,k1-1:kend-1) ...,
                +8.*uh(j1+1:jend+1,k1-1:kend-1)-uh(j1+2:jend+2,k1-1:kend-1)) ...,
                +8.*Yh(j1:jend,k1+1:kend+1).* ...,
                (uh(j1-2:jend-2,k1+1:kend+1)-8.*uh(j1-1:jend-1,k1+1:kend+1) ...,
                +8.*uh(j1+1:jend+1,k1+1:kend+1)-uh(j1+2:jend+2,k1+1:kend+1)) ...,
                -Yh(j1:jend,k1+2:kend+2).* ...,
                (uh(j1-2:jend-2,k1+2:kend+2)-8.*uh(j1-1:jend-1,k1+2:kend+2) ...,
                +8.*uh(j1+1:jend+1,k1+2:kend+2)-uh(j1+2:jend+2,k1+2:kend+2)) )...,
                ...,
                -(1./(144.*h1.^2)).*Yh_Inner.^2.*( uh(j1-4:jend-4,k1:kend) ...,
                -16.*uh(j1-3:jend-3,k1:kend)+64.*uh(j1-2:jend-2,k1:kend) ...,
                +16.*uh(j1-1:jend-1,k1:kend)-130.*uh(j1:jend,k1:kend) ...,
                +16.*uh(j1+1:jend+1,k1:kend)+64.*uh(j1+2:jend+2,k1:kend) ...,
                -16.*uh(j1+3:jend+3,k1:kend)+uh(j1+4:jend+4,k1:kend) ) ;
            
        elseif strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,II}') == 1 
            result=(-1./(12.*h2.^2)).*Xh_Inner.^2.*( -1.*uh(j1:jend,k1-2:kend-2) ...,
                +16.*uh(j1:jend,k1-1:kend-1)-30.*uh(j1:jend,k1:kend) ...，
                +16.*uh(j1:jend,k1+1:kend+1)-uh(j1:jend,k1+2:kend+2) ) ...,
                ...,
                -(1./(12.*h1.^2)).*Yh_Inner.^2.*( -1.*uh(j1-2:jend-2,k1:kend) ...,
                +16.*uh(j1-1:jend-1,k1:kend)-30.*uh(j1:jend,k1:kend) ...,
                +16.*uh(j1+1:jend+1,k1:kend)-uh(j1+2:jend+2,k1:kend) ) ...,
                ...,
                +Xh_Inner./(12.*h1).*(uh(j1-2:jend-2,k1:kend)-8.*uh(j1-1:jend-1,k1:kend) ...,
                +8.*uh(j1+1:jend+1,k1:kend)-uh(j1+2:jend+2,k1:kend)) ...,
                ...,
                +Yh_Inner./(12.*h2).*(uh(j1:jend,k1-2:kend-2)-8.*uh(j1:jend,k1-1:kend-1) ...，
                +8.*uh(j1:jend,k1+1:kend+1)-uh(j1:jend,k1+2:kend+2)) ...,
                ...,
                +Xh_Inner.*Yh_Inner./(72.*h1.*h2).*(  ...,
                (uh(j1-2:jend-2,k1-2:kend-2)-8.*uh(j1-2:jend-2,k1-1:kend-1) ...,
                +8.*uh(j1-2:jend-2,k1+1:kend+1)-uh(j1-2:jend-2,k1+2:kend+2)) ...,
                -8.*(uh(j1-1:jend-1,k1-2:kend-2)-8.*uh(j1-1:jend-1,k1-1:kend-1) ...,
                +8.*uh(j1-1:jend-1,k1+1:kend+1)-uh(j1-1:jend-1,k1+2:kend+2)) ...,
                +8.*(uh(j1+1:jend+1,k1-2:kend-2)-8.*uh(j1+1:jend+1,k1-1:kend-1) ...,
                +8.*uh(j1+1:jend+1,k1+1:kend+1)-uh(j1+1:jend+1,k1+2:kend+2)) ...,
                -(uh(j1+2:jend+2,k1-2:kend-2)-8.*uh(j1+2:jend+2,k1-1:kend-1) ...,
                +8.*uh(j1+2:jend+2,k1+1:kend+1)-uh(j1+2:jend+2,k1+2:kend+2)) ) ;
            
        elseif strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,III}') == 1
            result=(-1./(144.*h2.^2)).*Xh_Inner.^2.*( uh(j1:jend,k1-4:kend-4) ...,
                -16.*uh(j1:jend,k1-3:kend-3)+64.*uh(j1:jend,k1-2:kend-2) ...，
                +16.*uh(j1:jend,k1-1:kend-1)-130.*uh(j1:jend,k1:kend) ...，
                +16.*uh(j1:jend,k1+1:kend+1)+64.*uh(j1:jend,k1+2:kend+2) ...,
                -16.*uh(j1:jend,k1+3:kend+3)+uh(j1:jend,k1+4:kend+4) ) ...,
                ...,
                -(1./(144.*h1.^2)).*Yh_Inner.^2.*( uh(j1-4:jend-4,k1:kend) ...,
                -16.*uh(j1-3:jend-3,k1:kend)+64.*uh(j1-2:jend-2,k1:kend) ...,
                +16.*uh(j1-1:jend-1,k1:kend)-130.*uh(j1:jend,k1:kend) ...,
                +16.*uh(j1+1:jend+1,k1:kend)+64.*uh(j1+2:jend+2,k1:kend) ...,
                -16.*uh(j1+3:jend+3,k1:kend)+uh(j1+4:jend+4,k1:kend) ) ...,
                ...,
                +Xh_Inner./(12.*h1).*(uh(j1-2:jend-2,k1:kend)-8.*uh(j1-1:jend-1,k1:kend) ...,
                +8.*uh(j1+1:jend+1,k1:kend)-uh(j1+2:jend+2,k1:kend)) ...,
                ...,
                +Yh_Inner./(12.*h2).*(uh(j1:jend,k1-2:kend-2)-8.*uh(j1:jend,k1-1:kend-1) ...，
                +8.*uh(j1:jend,k1+1:kend+1)-uh(j1:jend,k1+2:kend+2)) ...,
                ...,
                +Xh_Inner.*Yh_Inner./(72.*h1.*h2).*(  ...,
                (uh(j1-2:jend-2,k1-2:kend-2)-8.*uh(j1-2:jend-2,k1-1:kend-1) ...,
                +8.*uh(j1-2:jend-2,k1+1:kend+1)-uh(j1-2:jend-2,k1+2:kend+2)) ...,
                -8.*(uh(j1-1:jend-1,k1-2:kend-2)-8.*uh(j1-1:jend-1,k1-1:kend-1) ...,
                +8.*uh(j1-1:jend-1,k1+1:kend+1)-uh(j1-1:jend-1,k1+2:kend+2)) ...,
                +8.*(uh(j1+1:jend+1,k1-2:kend-2)-8.*uh(j1+1:jend+1,k1-1:kend-1) ...,
                +8.*uh(j1+1:jend+1,k1+1:kend+1)-uh(j1+1:jend+1,k1+2:kend+2)) ...,
                -(uh(j1+2:jend+2,k1-2:kend-2)-8.*uh(j1+2:jend+2,k1-1:kend-1) ...,
                +8.*uh(j1+2:jend+2,k1+1:kend+1)-uh(j1+2:jend+2,k1+2:kend+2)) ) ;
        end
    end

%% 点坐标索引 =================================================================================
    function [result] = Index(j,k)
        result=(k+(m-1)).*(number_of_x_points)+j+m;
    end

%% 稀疏矩阵生成 =================================================================================
    function [result] = generate_A(coe2)
        if coe2 == 1  % 1/tau.*uh的系数矩阵
            Ar=Index(J_Index,K_Index);
            Ac=Index(J_Index,K_Index);
            Anum=ones(number_of_x_points,number_of_y_points);
            Anum(j1:jend,k1:kend)=1./tau.*Anum(j1:jend,k1:kend);
            result=sparse(Ar,Ac,Anum,number_of_x_points.*number_of_y_points, ...,
                number_of_x_points.*number_of_y_points);

        elseif coe2 == 2  % tau.*Delta_h(uh)的系数矩阵
            Ar=ones(number_of_x_points,number_of_y_points,9);
            Ac=ones(number_of_x_points,number_of_y_points,9);
            Anum=zeros(number_of_x_points,number_of_y_points,9);
            
            Point_index=[ 0,0;  -1,0;   0,-1;   1,0;   0,1;
                -2,0;   0,-2;   2,0;   0,2];
            Auxiliary_metrix=[-30*r1/12-30*r2/12, ...,
                16*r1/12, 16*r2/12, 16*r1/12, 16*r2/12, ...,
                (-1)*r1/12, (-1)*r2/12, (-1)*r1/12, (-1)*r2/12];

            Ar(:,:,1)=Index(J_Index,K_Index);
            Ac(:,:,1)=Index(J_Index,K_Index);
            Anum(:,:,1)=ones(number_of_x_points,number_of_y_points);
            for coe3 = 1:9
                Ar(j1:jend,k1:kend,coe3) = ...，
                    Index(J_Index(j1:jend,k1:kend),K_Index(j1:jend,k1:kend));
                Ac(j1:jend,k1:kend,coe3) = ...，
                    Index(J_Index(j1+Point_index(coe3,1):jend+Point_index(coe3,1),k1:kend), ...,
                    K_Index(j1:jend,k1+Point_index(coe3,2):kend+Point_index(coe3,2)));
                Anum(j1:jend,k1:kend,coe3) = Auxiliary_metrix(coe3);
            end
            Ar_new=reshape(Ar,number_of_x_points,[]);
            Ac_new=reshape(Ac,number_of_x_points,[]);
            Anum_new=reshape(Anum,number_of_x_points,[]);
            result=sparse(Ar_new,Ac_new,Anum_new,number_of_x_points.*number_of_y_points, ...,
                number_of_x_points.*number_of_y_points);

        elseif coe2 == 3 % tau.*uh的系数矩阵
            Ar=Index(J_Index,K_Index);
            Ac=Index(J_Index,K_Index);
            Anum=ones(number_of_x_points,number_of_y_points);
            Anum(j1:jend,k1:kend)=tau.*Anum(j1:jend,k1:kend);
            result=sparse(Ar,Ac,Anum,number_of_x_points.*number_of_y_points, ...,
                number_of_x_points.*number_of_y_points);            

        elseif coe2 == 4 % 1i.*L_zh(uh)的系数矩阵
            Ar=ones(number_of_x_points,number_of_y_points,8);
            Ac=ones(number_of_x_points,number_of_y_points,8);
            Anum=zeros(number_of_x_points,number_of_y_points,8);
            
            Point_index=[-1,0;   0,-1;   1,0;   0,1;
                -2,0;   0,-2;   2,0;   0,2];
            Auxiliary_metrix=[8/(12*h1), (-8)/(12*h2), (-8)/(12*h1), 8/(12*h2), ...,
                (-1)/(12*h1), 1/(12*h2), 1/(12*h1), (-1)/(12*h2)];

            Ar(:,:,1)=Index(J_Index,K_Index);
            Ac(:,:,1)=Index(J_Index,K_Index);
            Anum(:,:,1)=ones(number_of_x_points,number_of_y_points);
            for coe3 = 1:8
                Ar(j1:jend,k1:kend,coe3) = ...，
                    Index(J_Index(j1:jend,k1:kend),K_Index(j1:jend,k1:kend));
                Ac(j1:jend,k1:kend,coe3) = ...，
                    Index(J_Index(j1+Point_index(coe3,1):jend+Point_index(coe3,1),k1:kend), ...,
                    K_Index(j1:jend,k1+Point_index(coe3,2):kend+Point_index(coe3,2)));
                if Point_index(coe3,1)~=0 && Point_index(coe3,2)==0
                    Anum(j1:jend,k1:kend,coe3) = Yh_Inner.*Auxiliary_metrix(coe3);
                else
                    Anum(j1:jend,k1:kend,coe3) = Xh_Inner.*Auxiliary_metrix(coe3);
                end
            end
            Ar_new=reshape(Ar,number_of_x_points,[]);
            Ac_new=reshape(Ac,number_of_x_points,[]);
            Anum_new=reshape(Anum,number_of_x_points,[]);
            result=sparse(Ar_new,Ac_new,Anum_new,number_of_x_points.*number_of_y_points, ...,
                number_of_x_points.*number_of_y_points);            

        elseif coe2 == 5 % L_zh_2(uh)的系数矩阵
            if strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,I}') == 1
                Ar=ones(number_of_x_points,number_of_y_points,33);
                Ac=ones(number_of_x_points,number_of_y_points,33);
                Anum=zeros(number_of_x_points,number_of_y_points,33);

                Point_index=[0,0; ...,
                    -1,0; -1,-1;  0,-1;  1,-1;   1,0;   1,1;   0,1;  -1,1; ...,
                    -2,0; -2,-1; -2,-2; -1,-2;  0,-2;  1,-2;  2,-2;  2,-1; ...,
                     2,0;   2,1;   2,2;   1,2;   0,2;  -1,2;  -2,2;  -2,1; ...,
                    -3,0;  0,-3;   3,0;   0,3;  -4,0;  0,-4;   4,0;   0,4];
                Auxiliary_metrix=[130, ...,
                    -16, 64,-16,-64,-16, 64,-16,-64 ...,
                    -64, -8,  1, -8,-64,  8, -1,  8, ...,
                    -64, -8,  1, -8,-64,  8, -1,  8, ...,
                     16, 16, 16, 16, -1, -1, -1, -1]./144;

                Ar(:,:,1) = Index(J_Index,K_Index);
                Ac(:,:,1) = Index(J_Index,K_Index);
                Anum(:,:,1) = ones(number_of_x_points,number_of_y_points);
                Anum(j1:jend,k1:kend,1) = Auxiliary_metrix(1) ...,
                    .*(Xh_Inner.^2./h2.^2+Yh_Inner.^2./h1.^2);
                for coe3 = 2:33
                    Ar(j1:jend,k1:kend,coe3) = ...，
                        Index(J_Index(j1:jend,k1:kend),K_Index(j1:jend,k1:kend));
                    Ac(j1:jend,k1:kend,coe3) = ...，
                        Index(J_Index(j1+Point_index(coe3,1):jend+Point_index(coe3,1),k1:kend), ...,
                        K_Index(j1:jend,k1+Point_index(coe3,2):kend+Point_index(coe3,2)));
                    
                    if Point_index(coe3,1)~=0 && Point_index(coe3,2)~=0
                        Anum(j1:jend,k1:kend,coe3) ...，
                            =Auxiliary_metrix(coe3)./(h1.*h2).*(Yh_Inner.* ...,
                            Xh(j1+Point_index(coe3,1):jend+Point_index(coe3,1),k1:kend) ...,
                            +Xh_Inner.*Yh(j1:jend,k1+Point_index(coe3,2):kend+Point_index(coe3,2))); 

                    elseif abs(Point_index(coe3,1))>=3 && abs(Point_index(coe3,2))==0
                        Anum(j1:jend,k1:kend,coe3) ...,
                            =Auxiliary_metrix(coe3)./(h1.^2).*Yh_Inner.^2;

                    elseif abs(Point_index(coe3,1))==0 && abs(Point_index(coe3,2))>=3
                        Anum(j1:jend,k1:kend,coe3) ...,
                            =Auxiliary_metrix(coe3)./(h2.^2).*Xh_Inner.^2;

                    elseif (abs(Point_index(coe3,1))==1 || abs(Point_index(coe3,1))==2) ...,
                            && abs(Point_index(coe3,2))==0
                        Anum(j1:jend,k1:kend,coe3) ...,
                            =Auxiliary_metrix(coe3)./(h1.^2).*Yh_Inner.^2;

                    elseif (abs(Point_index(coe3,2))==1 || abs(Point_index(coe3,2))==2) ...,
                            && abs(Point_index(coe3,1))==0
                        Anum(j1:jend,k1:kend,coe3) ...,
                            =Auxiliary_metrix(coe3)./(h2.^2).*Xh_Inner.^2;
                    end
                end  
                
            elseif strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,II}') == 1
                Ar=ones(number_of_x_points,number_of_y_points,25);
                Ac=ones(number_of_x_points,number_of_y_points,25);
                Anum=zeros(number_of_x_points,number_of_y_points,25);

                Point_index=[ 0,0; ...,
                    -1,0; -1,-1;  0,-1;  1,-1;   1,0;   1,1;   0,1;  -1,1; ...,
                    -2,0; -2,-1; -2,-2; -1,-2;  0,-2;  1,-2;  2,-2;  2,-1; ...,
                     2,0;   2,1;   2,2;   1,2;   0,2;  -1,2;  -2,2;  -2,1];
                Auxiliary_metrix1=[30/12, ...,
                    -16/12, 64/72, -16/12, -64/72,  -16/12, 64/72, -16/12, -64/72, ...,
                    1/12,   -8/72,   1/72,  -8/72,    1/12,  8/72,  -1/72,   8/72, ...,
                    1/12,   -8/72,   1/72,  -8/72,    1/12,  8/72,  -1/72,   8/72];
                Auxiliary_metrix2=[0, ...,
                    -8/12, 0, -8/12, 0,  8/12, 0, 8/12, 0, ...,
                     1/12, 0,     0, 0,  1/12, 0,    0, 0, ...,
                    -1/12, 0,     0, 0, -1/12, 0,    0, 0];

                Ar(:,:,1) = Index(J_Index,K_Index);
                Ac(:,:,1) = Index(J_Index,K_Index);
                Anum(:,:,1) = ones(number_of_x_points,number_of_y_points);
                Anum(j1:jend,k1:kend,1) = Auxiliary_metrix1(1) ...,
                    .*(Xh_Inner.^2./h2.^2+Yh_Inner.^2./h1.^2);
                for coe3 = 2:25
                    Ar(j1:jend,k1:kend,coe3) = ...，
                        Index(J_Index(j1:jend,k1:kend),K_Index(j1:jend,k1:kend));
                    Ac(j1:jend,k1:kend,coe3) = ...，
                        Index(J_Index(j1+Point_index(coe3,1):jend+Point_index(coe3,1),k1:kend), ...,
                        K_Index(j1:jend,k1+Point_index(coe3,2):kend+Point_index(coe3,2)));
                    
                    if Point_index(coe3,1)~=0 && Point_index(coe3,2)~=0
                        Anum(j1:jend,k1:kend,coe3) ...，
                            =Auxiliary_metrix1(coe3)./(h1.*h2).*Yh_Inner.*Xh_Inner; 

                    elseif (abs(Point_index(coe3,1))==1 || abs(Point_index(coe3,1))==2) ...,
                            && abs(Point_index(coe3,2))==0
                        Anum(j1:jend,k1:kend,coe3) ...,
                            =Auxiliary_metrix1(coe3)./(h1.^2).*Yh_Inner.^2 ...,
                            +Auxiliary_metrix2(coe3)./h1.*Xh_Inner;

                    elseif (abs(Point_index(coe3,2))==1 || abs(Point_index(coe3,2))==2) ...,
                            && abs(Point_index(coe3,1))==0
                        Anum(j1:jend,k1:kend,coe3) ...,
                            =Auxiliary_metrix1(coe3)./(h2.^2).*Xh_Inner.^2 ...,
                            +Auxiliary_metrix2(coe3)./h2.*Yh_Inner;
                    end
                end
                
            elseif strcmp(Quadratic_rotation_term,'\mathcal{L}_z^{2,III}') == 1
                Ar=ones(number_of_x_points,number_of_y_points,33);
                Ac=ones(number_of_x_points,number_of_y_points,33);
                Anum=zeros(number_of_x_points,number_of_y_points,33);

                Point_index=[ 0,0; ...,
                    -1,0; -1,-1;  0,-1;  1,-1;   1,0;   1,1;   0,1;  -1,1; ...,
                    -2,0; -2,-1; -2,-2; -1,-2;  0,-2;  1,-2;  2,-2;  2,-1; ...,
                     2,0;   2,1;   2,2;   1,2;   0,2;  -1,2;  -2,2;  -2,1; ...,
                    -3,0;  0,-3;   3,0;   0,3;  -4,0;  0,-4;   4,0;   0,4];
                Auxiliary_metrix1=[130/144, ...,
                     -16/144,  64/72, -16/144, -64/72,  -16/144,  64/72, -16/144, -64/72, ...,
                     -64/144,  -8/72,    1/72,  -8/72,  -64/144,   8/72,   -1/72,   8/72, ...,
                     -64/144,  -8/72,    1/72,  -8/72,  -64/144,   8/72,   -1/72,   8/72, ...,
                      16/144, 16/144,  16/144, 16/144,   -1/144, -1/144,  -1/144, -1/144];
                Auxiliary_metrix2=[0, ...,
                    -8/12, 0, -8/12, 0,  8/12, 0, 8/12, 0, ...,
                     1/12, 0,     0, 0,  1/12, 0,    0, 0, ...,
                    -1/12, 0,     0, 0, -1/12, 0,    0, 0, ...,
                        0, 0,     0, 0,     0, 0,    0, 0];

                Ar(:,:,1) = Index(J_Index,K_Index);
                Ac(:,:,1) = Index(J_Index,K_Index);
                Anum(:,:,1) = ones(number_of_x_points,number_of_y_points);
                Anum(j1:jend,k1:kend,1) = Auxiliary_metrix1(1) ...,
                    .*(Xh_Inner.^2./h2.^2+Yh_Inner.^2./h1.^2);
                for coe3 = 2:33
                    Ar(j1:jend,k1:kend,coe3) = ...，
                        Index(J_Index(j1:jend,k1:kend),K_Index(j1:jend,k1:kend));
                    Ac(j1:jend,k1:kend,coe3) = ...，
                        Index(J_Index(j1+Point_index(coe3,1):jend+Point_index(coe3,1),k1:kend), ...,
                        K_Index(j1:jend,k1+Point_index(coe3,2):kend+Point_index(coe3,2)));
                    
                    if Point_index(coe3,1)~=0 && Point_index(coe3,2)~=0
                        Anum(j1:jend,k1:kend,coe3) ...，
                            =Auxiliary_metrix1(coe3)./(h1.*h2).*Yh_Inner.*Xh_Inner; 

                    elseif abs(Point_index(coe3,1))>=1 && abs(Point_index(coe3,2))==0
                        Anum(j1:jend,k1:kend,coe3) ...,
                            =Auxiliary_metrix1(coe3)./(h1.^2).*Yh_Inner.^2 ...,
                            +Auxiliary_metrix2(coe3)./h1.*Xh_Inner;

                    elseif abs(Point_index(coe3,2))>=1 && abs(Point_index(coe3,1))==0
                        Anum(j1:jend,k1:kend,coe3) ...,
                            =Auxiliary_metrix1(coe3)./(h2.^2).*Xh_Inner.^2 ...,
                            +Auxiliary_metrix2(coe3)./h2.*Yh_Inner;
                    end
                end
            end
            Ar_new=reshape(Ar,number_of_x_points,[]);
            Ac_new=reshape(Ac,number_of_x_points,[]);
            Anum_new=reshape(Anum,number_of_x_points,[]);
            result=sparse(Ar_new,Ac_new,Anum_new,number_of_x_points.*number_of_y_points, ...,
                number_of_x_points.*number_of_y_points);

        elseif coe2 == 6 % tau.*lambda.*非线性项的系数矩阵
            Ar=Index(J_Index,K_Index);
            Ac=Index(J_Index,K_Index);
            Anum=ones(number_of_x_points,number_of_y_points);
            if strcmp(scheme_type,'LCFD') == 1          
                Anum(j1:jend,k1:kend)=tau.*lambda.*abs(psin(j1:jend,k1:kend)).^2; 
            elseif strcmp(scheme_type,'NCFD-I') == 1
                Anum(j1:jend,k1:kend)=0.5.*tau.*lambda.*(abs(psin_minus(j1:jend,k1:kend)).^2 ...,
                    +abs(uh_s(j1:jend,k1:kend,1)).^2);
            elseif strcmp(scheme_type,'NCFD-II') == 1
                Anum(j1:jend,k1:kend)=0.125.*tau.*lambda.* ...,
                    (abs(psin_minus(j1:jend,k1:kend)+psin(j1:jend,k1:kend)).^2 ...,
                    +abs(uh_s(j1:jend,k1:kend,1)+psin(j1:jend,k1:kend)).^2);
            elseif strcmp(scheme_type,'RCFD') == 1
                Anum(j1:jend,k1:kend)=0.5.*tau.*lambda.* ...,
                    (abs(psin(j1:jend,k1:kend)).^2 ...,
                    +abs(uh_s(j1:jend,k1:kend,1)).^2);
            end
            result=sparse(Ar,Ac,Anum,number_of_x_points.*number_of_y_points, ...,
                    number_of_x_points.*number_of_y_points);
        end
    end

%% 右端项R生成 =================================================================================
    function [result] = computer_R(method_type)
        if strcmp(method_type,'LCFD') == 1 
            result = 2./tau.*(2.*psin(j1:jend,k1:kend)-psin_minus(j1:jend,k1:kend)) ...,
                +tau.*Delta_h(psin_minus)-tau.*(1+lambda.*abs(psin(j1:jend,k1:kend)).^2) ...,
                .*psin_minus(j1:jend,k1:kend) ...,
                -2i.*gamma.*L_zh(psin_minus)+tau.*gamma.^2.*L_zh_2(psin_minus);
        elseif strcmp(method_type,'NCFD-I') == 1 
            result = 2./tau.*(2.*psin(j1:jend,k1:kend)-psin_minus(j1:jend,k1:kend)) ...,
                +tau.*Delta_h(psin_minus)-tau.*(1 ...,
                +0.5.*lambda.*(abs(psin_minus(j1:jend,k1:kend)).^2 ...,
                +abs(uh_s(j1:jend,k1:kend,1)).^2)).*psin_minus(j1:jend,k1:kend) ...,
                -2i.*gamma.*L_zh(psin_minus)+tau.*gamma.^2.*L_zh_2(psin_minus);
        elseif strcmp(method_type,'NCFD-II') == 1
            result = 4./tau.*(2.*psin(j1:jend,k1:kend)-psin_minus(j1:jend,k1:kend)) ...,
                +tau.*Delta_h(psin_minus+2.*psin)-tau.*(1 ...,
                +0.125.*lambda.*(abs(psin_minus(j1:jend,k1:kend)+psin(j1:jend,k1:kend)).^2 ...,
                +abs(uh_s(j1:jend,k1:kend,1)+psin(j1:jend,k1:kend)).^2)) ...,
                .*(psin_minus(j1:jend,k1:kend)+2.*psin(j1:jend,k1:kend)) ...,
                -4i.*gamma.*L_zh(psin_minus)+tau.*gamma.^2.*L_zh_2(psin_minus+2.*psin);
        elseif strcmp(method_type,'RCFD') == 1
            result = 4./tau.*(2.*psin(j1:jend,k1:kend)-psin_minus(j1:jend,k1:kend)) ...,
                +tau.*Delta_h(psin_minus+2.*psin) ...,
                -tau.*(psin_minus(j1:jend,k1:kend)+2.*psin(j1:jend,k1:kend)) ...,
                -0.5.*tau.*lambda.*( (abs(psin_minus(j1:jend,k1:kend)).^2 ...,
                +abs(psin(j1:jend,k1:kend)).^2).*(psin_minus(j1:jend,k1:kend) ...,
                +psin(j1:jend,k1:kend))+( abs(psin(j1:jend,k1:kend)).^2 ...,
                +abs(uh_s(j1:jend,k1:kend,1)).^2 ).*psin(j1:jend,k1:kend) )...,
                -4i.*gamma.*L_zh(psin_minus)+tau.*gamma.^2.*L_zh_2(psin_minus+2.*psin);
        end
    end

%% 计算守恒律辅助函数 =================================================================================
    function [result] = Inner_product(u1,u2)
        result = h1.*h2.*sum(sum(u1.*conj(u2)));
    end
    function [result] = H1_semi_norm(uh)
        result = (4./3).*(h2./h1) ...,
            .*sum(sum(abs(uh(j1:jend+1,k1:kend)-uh(j1-1:jend,k1:kend)).^2)) ...,
            +(4./3).*(h1./h2) ...,
            .*sum(sum(abs(uh(j1:jend,k1:kend+1)-uh(j1:jend,k1-1:kend)).^2)) ...,
            -(1./3).*h1.*h2.*sum(sum( abs(uh(j1:jend+2,k1-1:kend+1)-uh(j1-2:jend,k1-1:kend+1)).^2./(2.*h1).^2 ...,
            +abs(uh(j1-1:jend+1,k1:kend+2)-uh(j1-1:jend+1,k1-2:kend)).^2./(2.*h2).^2 ));
    end
    function [result] = L2_norm(uh)
        result = h1.*h2.*sum(sum(abs(uh).^2));
    end
    function [result] = Lz_h_energy(uh)
        result=Xh(j1-2:jend+2,k1-2:kend+2)./(12.*h2) ...,
            .*( uh(j1-2:jend+2,k1-4:kend) ...,
            -8.*uh(j1-2:jend+2,k1-3:kend+1)+8.*uh(j1-2:jend+2,k1-1:kend+3) ...,
            -uh(j1-2:jend+2,k1:kend+4) ) ...,
            ...,
            -Yh(j1-2:jend+2,k1-2:kend+2)./(12.*h1) ...,
            .*( uh(j1-4:jend,k1-2:kend+2) ...,
            -8.*uh(j1-3:jend+1,k1-2:kend+2)+8.*uh(j1-1:jend+3,k1-2:kend+2) ...,
            -uh(j1:jend+4,k1-2:kend+2) ) ;
        result=-1i.*result;
    end


%% 守恒律计算 =================================================================================
    function [Mn,En] = Comput_Mass_Energy(uh_minus,uh)
        if strcmp(scheme_type,'LCFD') == 1
            Mn = imag( Inner_product( uh(j1:jend,k1:kend)- ...,
                uh_minus(j1:jend,k1:kend),uh_minus(j1:jend,k1:kend) )./tau ) ...,
                ...,
                -0.5.*gamma.*real( Inner_product( L_zh(uh),uh(j1:jend,k1:kend) ) ...,
                +Inner_product( L_zh(uh_minus),uh_minus(j1:jend,k1:kend) ) );

            En = L2_norm(uh(j1:jend,k1:kend)-uh_minus(j1:jend,k1:kend))./(tau.^2) ...,
                ...,
                +0.5.*( H1_semi_norm(uh)+H1_semi_norm(uh_minus) ...,
                +L2_norm(uh(j1:jend,k1:kend))+L2_norm(uh_minus(j1:jend,k1:kend)) ) ...,
                ...,
                +0.5.*lambda.*L2_norm(uh_minus(j1:jend,k1:kend).*uh(j1:jend,k1:kend)) ...,
                ...,
                -0.5.*gamma.^2.*(L2_norm(Lz_h_energy(uh))+L2_norm(Lz_h_energy(uh_minus)));
        
        elseif strcmp(scheme_type,'NCFD-I') == 1
            Mn = imag( Inner_product( uh(j1:jend,k1:kend)- ...,
                uh_minus(j1:jend,k1:kend),uh_minus(j1:jend,k1:kend) )./tau ) ...,
                ...,
                -0.5.*gamma.*real( Inner_product( L_zh(uh),uh(j1:jend,k1:kend) ) ...,
                +Inner_product( L_zh(uh_minus),uh_minus(j1:jend,k1:kend) ) );

            En = L2_norm(uh(j1:jend,k1:kend)-uh_minus(j1:jend,k1:kend))./(tau.^2) ...,
                ...,
                +0.5.*( H1_semi_norm(uh)+H1_semi_norm(uh_minus) ...,
                +L2_norm(uh(j1:jend,k1:kend))+L2_norm(uh_minus(j1:jend,k1:kend)) ) ...,
                ...,
                +0.25.*lambda.*( L2_norm(uh_minus(j1:jend,k1:kend).^2) ...,
                +L2_norm(uh(j1:jend,k1:kend).^2) ) ..., 
                ...,
                -0.5.*gamma.^2.*(L2_norm(Lz_h_energy(uh))+L2_norm(Lz_h_energy(uh_minus)));

        elseif strcmp(scheme_type,'NCFD-II') == 1
            uh_star = 0.5.*(uh+uh_minus);

            Mn = imag( Inner_product( uh(j1:jend,k1:kend)- ...,
                uh_minus(j1:jend,k1:kend),uh_minus(j1:jend,k1:kend) )./tau ) ...,
                ...,
                -gamma.*real( Inner_product( L_zh(uh_star),uh_star(j1:jend,k1:kend) ) );

            En = L2_norm(uh(j1:jend,k1:kend)-uh_minus(j1:jend,k1:kend))./(tau.^2) ...,
                ...,
                +H1_semi_norm(uh_star)+L2_norm(uh_star(j1:jend,k1:kend)) ...,
                ...,
                +0.5.*lambda.*L2_norm(uh_star(j1:jend,k1:kend).^2) ..., 
                ...,
                -gamma.^2.*L2_norm(Lz_h_energy(uh_star));

        elseif strcmp(scheme_type,'RCFD') == 1
            Mn = imag( Inner_product( u(j1:jend,k1:kend),uh(j1:jend,k1:kend) )) ...,
                ...,
                -gamma.*real( Inner_product( L_zh(uh),uh(j1:jend,k1:kend) ) );

            En = L2_norm(u(j1:jend,k1:kend)) ...,
                ...,
                +H1_semi_norm(uh)+L2_norm(uh(j1:jend,k1:kend)) ...,
                ...,
                +0.5.*lambda.*L2_norm(uh(j1:jend,k1:kend).^2) ..., 
                ...,
                -gamma.^2.*L2_norm(Lz_h_energy(uh));
        end
    end
end

