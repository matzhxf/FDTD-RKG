function [psi_N] = I_scheme_sparse(phi_0,phi_1,gamma,lambda,h1,h2,tau,~, ...,
    x_left,x_right,y_left,y_right,t_begin,t_end,plot_mod,line_mod)
%LCFD_SCHEME_SPARSE 此处显示有关此函数的摘要
%   此处显示详细说明
tol = 1e-12;
iteration_max = 100;
[x_left,x_right,y_left,y_right,t_begin,t_end, ...,
    phi_0,phi_1,gamma,lambda,h1,h2,tau,J,K,N] = Variable_setting(num);


x_inner_point=x_left+h1:h1:x_left+(J-1)*h1;
y_inner_point=y_left+h2:h2:y_left+(K-1)*h2;
x_all_point=[x_left-3.*h1,x_left-2.*h1,x_left-h1,x_left,x_inner_point, ...,
    x_right,x_right+h1,x_right+2.*h1,x_right+3.*h1];
y_all_point=[y_left-3.*h2,y_left-2.*h2,y_left-h2,y_left,y_inner_point, ...,
    y_right,y_right+h2,y_right+2.*h2,y_right+3.*h2];
t_point=t_begin:tau:t_end;
[Y_Inner,X_Inner]=meshgrid(y_inner_point,x_inner_point);
[Y_All,X_All]=meshgrid(y_all_point,x_all_point);
% [Delta_x,Delta_y,Nable_x,Nable_y,~,~,~,~,~,~,~,~]= ...,
%     generate_Matrix(J,K,h1,h2,x_inner_point,y_inner_point);
r1=tau./h1.^2; r2=tau./h2.^2;
[jbegin_inner,jend_inner,kbegin_inner,kend_inner]=deal(5,J+3,5,K+3);
[K_matrix,J_matrix]=meshgrid(-3:1:K+3,-3:1:J+3);

psin_minus=zeros(J+7,K+7);
psin=zeros(J+7,K+7); 
psin_add=zeros(J+7,K+7);  % 分别预设第n-1、n、n+1层，n值待定
phi1=zeros(J+7,K+7);

% n=1时
psin_minus(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
    =feval(phi_0,X_Inner,Y_Inner); % 载入第0层
Infinity_norm(1,1)=max(max(abs(psin_minus)));

phi1(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
    =feval(phi_1,X_Inner,Y_Inner);

% 计算第1层(此时n=1)
psin(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
    =psin_minus(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
    +tau.*phi1(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
    +0.5.*tau.^2.*(Delta_h(psin_minus)-(1+lambda ...,
    .*abs(psin_minus(jbegin_inner:jend_inner,kbegin_inner:kend_inner)).^2) ...,
    .*psin_minus(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
    +2i.*gamma.*L_zh(phi1)+gamma.^2.*L_zh_2(psin_minus) );
Infinity_norm(1,2)=max(max(abs(psin)));

if plot_mod==1
    % 质量、能量
    Mh=zeros(1,N);
    Error_Mh=zeros(1,N);
    Eh=zeros(1,N);
    Error_Eh=zeros(1,N);
    [Mh(1,1),Eh(1,1)]=Comput_Mass_Energy(psin_minus,psin);
elseif plot_mod==2
    Dynal=moviein(41);   %=========== 画动态图
    split_N=ceil(N./40);
    split_n=1;
    surf(X_All,Y_All,abs(psin_minus));
    x_tick=x_left:2:x_right; y_tick=y_left:2:y_right;    
    axis([x_left x_right y_left y_right 0 1.5]);
    set(gca,'XTick',x_tick,'YTick',y_tick);
    set(gca,'Fontsize',15);
    xlabel('\itx','Fontname','Times New Roman','FontSize',20); 
    ylabel('\ity','Fontname','Times New Roman','FontSize',20);
    colorbar; shading interp;  view(0,90);
    title('{\itt} = 0, {\itL_z^{2,I}}','Fontname','Times New Roman','FontSize',20);
    Dynal(:,1)=getframe;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
Ar=ones(J+7,K+7,33);
Ac=ones(J+7,K+7,33);
Anum=zeros(J+7,K+7,33);
Point_index=[ 0,0;  -1,0; -1,-1;  0,-1;  1,-1;   1,0;   1,1;
    0,1;  -1,1;  -2,0; -2,-1; -2,-2; -1,-2;  0,-2;  1,-2;  2,-2;  2,-1;
    2,0;   2,1;   2,2;   1,2;   0,2;  -1,2;  -2,2;  -2,1;  -3,0;  0,-3;
    3,0;   0,3;  -4,0;  0,-4;   4,0;   0,4];
Coe1=[30*r1/12+30*r2/12, ...,
    -16*r1/12,  0,  -16*r2/12,  0,  -16*r1/12,  0,  -16*r2/12,  0, ...,
    r1/12,      0,          0,  0,      r2/12,  0,          0,  0, ...,
    r1/12,      0,          0,  0,      r2/12,  0,          0,  0, ...,
    0,          0,          0,  0,          0,  0,          0,  0];
Coe2=[0, ...,
    8/(12*h1),     0,  (-8)/(12*h2),  0,  (-8)/(12*h1),  0,  8/(12*h2),  0, ...,
    (-1)/(12*h1),  0,             0,  0,     1/(12*h2),  0,          0,  0, ...,
    1/(12*h1),     0,             0,  0,  (-1)/(12*h2),  0,          0,  0, ...,
    0,             0,             0,  0,             0,  0,          0,  0] ...,
    .*(-2).*gamma;
Coe3=[130,-16,64,-16,-64,-16,64,-16,-64,-64,-8,1,-8,-64,8,-1,8, ...,
    -64,-8,1,-8,-64,8,-1,8,16,16,16,16,-1,-1,-1,-1]./144;
Ar(:,:,1)=Index(J_matrix,K_matrix);
Ac(:,:,1)=Index(J_matrix,K_matrix);
Anum(:,:,1)=ones(J+7,K+7);
for num=2:33
    Ar(jbegin_inner:jend_inner,kbegin_inner:kend_inner,num) ...，
        =Index(Mat_tran(J_matrix,0,0),Mat_tran(K_matrix,0,0));
    Ac(jbegin_inner:jend_inner,kbegin_inner:kend_inner,num) ...，
        =Index(Mat_tran(J_matrix,Point_index(num,1),0), ...,
        Mat_tran(K_matrix,0,Point_index(num,2)));
    if Point_index(num,1)~=0 && Point_index(num,2)~=0
        Anum(jbegin_inner:jend_inner,kbegin_inner:kend_inner,num) ...，
            =-1.*tau.*gamma.^2.*Coe3(num)./(h1.*h2) ...,
            .*(Y_Inner.*Mat_tran(X_All,Point_index(num,1),0) ...,
            +X_Inner.*Mat_tran(Y_All,0,Point_index(num,2)));
    elseif abs(Point_index(num,1))>=3 && abs(Point_index(num,2))==0
        Anum(jbegin_inner:jend_inner,kbegin_inner:kend_inner,num) ...,
            =-1.*tau.*gamma.^2.*Coe3(num)./(h1.^2).*Y_Inner.^2;
    elseif abs(Point_index(num,1))==0 && abs(Point_index(num,2))>=3
        Anum(jbegin_inner:jend_inner,kbegin_inner:kend_inner,num) ...,
            =-1.*tau.*gamma.^2.*Coe3(num)./(h2.^2).*X_Inner.^2;
    elseif (abs(Point_index(num,1))==1 || abs(Point_index(num,1))==2) ...,
            && abs(Point_index(num,2))==0
        Anum(jbegin_inner:jend_inner,kbegin_inner:kend_inner,num) ...,
            =Coe1(num)+Coe2(num).*Y_Inner ...,
            -1.*tau.*gamma.^2.*Coe3(num)./(h1.^2).*Y_Inner.^2;
    elseif (abs(Point_index(num,2))==1 || abs(Point_index(num,2))==2) ...,
            && abs(Point_index(num,1))==0
        Anum(jbegin_inner:jend_inner,kbegin_inner:kend_inner,num) ...,
            =Coe1(num)+Coe2(num).*X_Inner ...,
            -1.*tau.*gamma.^2.*Coe3(num)./(h2.^2).*X_Inner.^2;
    end
end
Ar_new=reshape(Ar,J+7,[]);
Ac_new=reshape(Ac,J+7,[]);
Anum_new=reshape(Anum,J+7,[]);
A=sparse(Ar_new,Ac_new,Anum_new,(J+7)*(K+7),(J+7)*(K+7));

%% 计算第n+1层(n>=1)
R=zeros(J+7,K+7);
for n=1:N-1
    R(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
        =2./tau.*(2.*psin(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
        -psin_minus(jbegin_inner:jend_inner,kbegin_inner:kend_inner)) ...,
        +tau.*Delta_h(psin_minus)-tau.*(1+lambda ...,
        .*abs(psin(jbegin_inner:jend_inner,kbegin_inner:kend_inner)).^2) ...,
        .*psin_minus(jbegin_inner:jend_inner,kbegin_inner:kend_inner) ...,
        -2i.*gamma.*L_zh(psin_minus)+tau.*gamma.^2.*L_zh_2(psin_minus);
    R_vector=reshape(R,[],1);
    Anum(jbegin_inner:jend_inner,kbegin_inner:kend_inner,1) ...,
        =2./tau+tau.*(1+lambda ...,
        .*abs(psin(jbegin_inner:jend_inner,kbegin_inner:kend_inner)).^2) ...,
        +Coe1(1)-tau.*gamma.^2.*Coe3(1).*(X_Inner.^2./h2.^2 ...,
        +Y_Inner.^2./h1.^2 );
    for k=1:K-1
        for j=1:J-1
            A(Index(j,k),Index(j,k))=2./tau+tau.*(1+lambda ...,
                .*abs(psin(jbegin_inner+j-1,kbegin_inner+k-1)).^2) ...,
                +Coe1(1)-tau.*gamma.^2.*Coe3(1).*(x_inner_point(j).^2./h2.^2 ...,
                +y_inner_point(k).^2./h1.^2 );
        end
    end
%     tiaojianshu=cond(A,2)
    psin_add(:,:)=reshape(A\R_vector,J+7,K+7);
    Infinity_norm(1,n+2)=max(max(abs(psin_add)));

    psin_minus=psin;
    psin=psin_add;
    if plot_mod==1
        [Mh(1,n+1),Eh(1,n+1)]=Comput_Mass_Energy(psin_minus,psin);
        Error_Mh(1,n+1)=abs(Mh(1,n+1)-Mh(1,1))./Mh(1,1); 
        Error_Eh(1,n+1)=abs(Eh(1,n+1)-Eh(1,1))./Eh(1,1);
    elseif plot_mod==2 && mod(n+1,split_N)==0   %===========画动态图
        split_n=split_n+1;
        surf(X_All,Y_All,abs(psin));
        x_tick=x_left:2:x_right; y_tick=y_left:2:y_right;
        axis([x_left x_right y_left y_right 0 1.5]);
        set(gca,'XTick',x_tick,'YTick',y_tick);
        set(gca,'Fontsize',15);
        xlabel('\itx','Fontname','Times New Roman','FontSize',20);
        ylabel('\ity','Fontname','Times New Roman','FontSize',20);
        colorbar; shading interp;  view(0,90);
        title(['{\itt} = ',num2str((n+1).*tau),', {\itL_z^{2,I}}'], ...,
            'Fontname','Times New Roman','FontSize',20);
        Dynal(:,split_n)=getframe;
        stop=1;
    end
end
psi_N=psin( (jbegin_inner-1):(jend_inner+1),(kbegin_inner-1):(kend_inner+1));

%%  画图
if plot_mod==3
    plot(t_point,Infinity_norm,line_mod,'LineWidth',2);
    axis([0 t_end 0 5]);
    xlabel('\itt','Fontname','Times New Roman','FontSize',25); 
    ylabel('||{\it\psi^{ n}}||_{\itL^{\infty}}','Fontname','Times New Roman','FontSize',20);
    title('\itL_{z,I}^{2,h}','Fontname','Times New Roman','FontSize',20)
    % 坐标轴刻度设置
    y_tick=0:1:5;
    axis([t_begin t_end y_tick(1) y_tick(end)]);
    set(gca,'YTick',y_tick);
    set(gca,'Fontsize',15);

elseif plot_mod==1
    t_point(1)=[];
    % 质量图像
    plot(t_point,Mh,'Linewidth',1.2);
    xlabel('\itt','Fontname','Times New Roman','FontSize',25);
    ylabel('$Q^{n}$','Interpreter','Latex','FontSize',25);
    % 坐标轴刻度设置
    y_tick=0.7852:0.0001:0.7856;
    axis([t_begin t_end y_tick(1) y_tick(end)]);
    set(gca,'YTick',y_tick);
    set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.4f'));
    set(gca,'Fontsize',15);

    % 质量误差图像
    figure
    plot(t_point,Error_Mh);
    xlabel('\itt','Fontname','Times New Roman','FontSize',25); 
    ylabel('$|Q^{n}-Q^{0}|/Q^{0}$','Interpreter','Latex','FontSize',25);
    % 坐标轴刻度设置
%     y_tick=(-3:1:3).*1e-15;
%     axis([t_begin t_end y_tick(1) y_tick(end)]);
%     set(gca,'YTick',y_tick);
    set(gca,'Fontsize',15);

    % 能量图像
    figure
    plot(t_point,Eh,'Linewidth',1.2);
    xlabel('\itt','Fontname','Times New Roman','FontSize',25); 
    ylabel('$E^{n}$','Interpreter','Latex','FontSize',25);
    % 坐标轴刻度设置
    y_tick=5.4965:0.0001:5.4969;
    axis([t_begin t_end y_tick(1) y_tick(end)]);
    set(gca,'YTick',y_tick);
    set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.4f'));
    set(gca,'Fontsize',15);
    
    % 能量误差图像
    figure
    plot(t_point,Error_Eh);
    xlabel('\itt','Fontname','Times New Roman','FontSize',25); 
    ylabel('$|E^{n}-E^{0}|/E^{0}$','Interpreter','Latex','FontSize',25);
    % 坐标轴刻度设置
%     y_tick=(-3:1:3).*1e-14;
%     axis([t_begin t_end y_tick(1) y_tick(end)]);
%     set(gca,'YTick',y_tick);
    set(gca,'Fontsize',15);
end

%% 算子函数
    function [result] = Delta_h(uh)
        result=1./(12.*h1.^2).*( -1.*Mat_tran(uh,-2,0)+16.*Mat_tran(uh,-1,0) ...,
            -30.*Mat_tran(uh,0,0)+16.*Mat_tran(uh,1,0)-Mat_tran(uh,2,0) ) ...,
            +1./(12.*h2.^2).*( -1.*Mat_tran(uh,0,-2)+16.*Mat_tran(uh,0,-1) ...,
            -30.*Mat_tran(uh,0,0)+16.*Mat_tran(uh,0,1)-Mat_tran(uh,0,2) );
    end
    function [result] = L_zh(uh)
        result=X_Inner./(12.*h2).*( Mat_tran(uh,0,-2)-8.*Mat_tran(uh,0,-1) ...,
            +8.*Mat_tran(uh,0,1)-Mat_tran(uh,0,2) ) ...,
            -Y_Inner./(12.*h1).*( Mat_tran(uh,-2,0)-8.*Mat_tran(uh,-1,0) ...,
            +8.*Mat_tran(uh,1,0)-Mat_tran(uh,2,0) ) ;
        result=-1i.*result;
    end
    function [result] = L_zh_2(uh)
        result=(-1./(144.*h2.^2)).*X_Inner.^2.*( Mat_tran(uh,0,-4) ...,
            -16.*Mat_tran(uh,0,-3)+64.*Mat_tran(uh,0,-2)+16.*Mat_tran(uh,0,-1) ...,
            -130.*Mat_tran(uh,0,0)+16.*Mat_tran(uh,0,1)+64.*Mat_tran(uh,0,2) ...,
            -16.*Mat_tran(uh,0,3)+Mat_tran(uh,0,4) ) ...,
             ...,
            +Y_Inner./(144.*h1.*h2).*( Mat_tran(X_All,-2,0).*(Mat_tran(uh,-2,-2) ...,
            -8.*Mat_tran(uh,-2,-1)+8.*Mat_tran(uh,-2,1)-Mat_tran(uh,-2,2)) ...,
            -8.*Mat_tran(X_All,-1,0).*(Mat_tran(uh,-1,-2) ...,
            -8.*Mat_tran(uh,-1,-1)+8.*Mat_tran(uh,-1,1)-Mat_tran(uh,-1,2)) ...,
            +8.*Mat_tran(X_All,1,0).*(Mat_tran(uh,1,-2) ...,
            -8.*Mat_tran(uh,1,-1)+8.*Mat_tran(uh,1,1)-Mat_tran(uh,1,2)) ...,
            -Mat_tran(X_All,2,0).*(Mat_tran(uh,2,-2) ...,
            -8.*Mat_tran(uh,2,-1)+8.*Mat_tran(uh,2,1)-Mat_tran(uh,2,2)) ) ...,
             ...,
            +X_Inner./(144.*h1.*h2).*( Mat_tran(Y_All,0,-2).*(Mat_tran(uh,-2,-2) ...,
            -8.*Mat_tran(uh,-1,-2)+8.*Mat_tran(uh,1,-2)-Mat_tran(uh,2,-2)) ...,
            -8.*Mat_tran(Y_All,0,-1).*(Mat_tran(uh,-2,-1) ...,
            -8.*Mat_tran(uh,-1,-1)+8.*Mat_tran(uh,1,-1)-Mat_tran(uh,2,-1)) ...,
            +8.*Mat_tran(Y_All,0,1).*(Mat_tran(uh,-2,1) ...,
            -8.*Mat_tran(uh,-1,1)+8.*Mat_tran(uh,1,1)-Mat_tran(uh,2,1)) ...,
            -Mat_tran(Y_All,0,2).*(Mat_tran(uh,-2,2) ...,
            -8.*Mat_tran(uh,-1,2)+8.*Mat_tran(uh,1,2)-Mat_tran(uh,2,2)) ) ...,
             ...,
            -(1./(144.*h1.^2)).*Y_Inner.^2.*( Mat_tran(uh,-4,0) ...,
            -16.*Mat_tran(uh,-3,0)+64.*Mat_tran(uh,-2,0)+16.*Mat_tran(uh,-1,0) ...,
            -130.*Mat_tran(uh,0,0)+16.*Mat_tran(uh,1,0)+64.*Mat_tran(uh,2,0) ...,
            -16.*Mat_tran(uh,3,0)+Mat_tran(uh,4,0) ) ;
    end
    function [result] = Mat_tran(uh,j_var,k_var)
        result=uh((jbegin_inner+j_var):(jend_inner+j_var), ...,
            (kbegin_inner+k_var):(kend_inner+k_var));
    end
    function [result] = Index(j,k)
        result=(k+3).*(J+7)+j+4;
    end

%% 质量能量函数
    function [Mn,En] = Comput_Mass_Energy(uh_minus,uh)
        Mn=imag( h1.*h2./tau.*sum(sum( (uh(5:J+3,5:K+3)-uh_minus(5:J+3,5:K+3)) ...,
            .*conj(uh_minus(5:J+3,5:K+3)) )) ) ...,
            ...,
            -0.5.*gamma.*real( ...，
            h1.*h2.*sum(sum(-1i.*(X_Inner./(12.*h2) ...,
            .*(uh(5:J+3,3:K+1)-8.*uh(5:J+3,4:K+2)+8.*uh(5:J+3,6:K+4)-uh(5:J+3,7:K+5)) ...,
            -Y_Inner./(12.*h1) ...,
            .*(uh(3:J+1,5:K+3)-8.*uh(4:J+2,5:K+3)+8.*uh(6:J+4,5:K+3)-uh(7:J+5,5:K+3))) ...,
            .*conj(uh(5:J+3,5:K+3)) )) ...,
            ...,
            +h1.*h2.*sum(sum(-1i.*(X_Inner./(12.*h2) ...,
            .*(uh_minus(5:J+3,3:K+1)-8.*uh_minus(5:J+3,4:K+2) ...,
            +8.*uh_minus(5:J+3,6:K+4)-uh_minus(5:J+3,7:K+5)) ...,
            -Y_Inner./(2.*h1) ...,
            .*(uh_minus(3:J+1,5:K+3)-8.*uh_minus(4:J+2,5:K+3) ...,
            +8.*uh_minus(6:J+4,5:K+3)-uh_minus(7:J+5,5:K+3))) ...,
            .*conj(uh(5:J+3,5:K+3)) )) );

        En=h1.*h2./(tau.^2).*sum(sum(abs(uh(5:J+3,5:K+3)-uh_minus(5:J+3,5:K+3)).^2)) ...,
            ...,   
            +0.5.*( 4./3.*(h2./h1).*sum(sum(abs(uh(5:J+4,5:K+3)-uh(4:J+3,5:K+3)).^2)) ...,
            +4./3.*(h1./h2).*sum(sum(abs(uh(5:J+3,5:K+4)-uh(5:J+3,4:K+3)).^2))+ ...,
            h1.*h2.*sum(sum(abs(uh(5:J+3,5:K+3)).^2+abs(uh_minus(5:J+3,5:K+3)).^2)) ...,
            +4./3.*(h2./h1).*sum(sum(abs(uh_minus(5:J+4,5:K+3)-uh_minus(4:J+3,5:K+3)).^2)) ...,
            +4./3.*(h1./h2).*sum(sum(abs(uh_minus(5:J+3,5:K+4)-uh_minus(5:J+3,4:K+3)).^2))  ...,
            -h1.*h2./3.*sum(sum( abs(uh(5:J+5,4:K+4)-uh(3:J+3,4:K+4)).^2./(2.*h1).^2 ...,
            +abs(uh(4:J+4,5:K+5)-uh(4:J+4,3:K+3)).^2./(2.*h2).^2 ))  ...,
            -h1.*h2./3.*sum(sum( abs(uh_minus(5:J+5,4:K+4)-uh_minus(3:J+3,4:K+4)).^2./(2.*h1).^2 ...,
            +abs(uh_minus(4:J+4,5:K+5)-uh_minus(4:J+4,3:K+3)).^2./(2.*h2).^2 )) ) ...，
            ...,
            +0.5.*lambda.*h1.*h2.*sum(sum( ...,
            abs(uh_minus(5:J+3,5:K+3)).^2.*abs(uh(5:J+3,5:K+3)).^2 )) ...,
            ...,
            -0.5.*gamma.^2.*h1.*h2.*( sum(sum(abs(X_All(3:J+5,3:K+5)./(12.*h2) ...,
            .*(uh(3:J+5,1:K+3)-8.*uh(3:J+5,2:K+4)+8.*uh(3:J+5,4:K+6)-uh(3:J+5,5:K+7)) ...,
            -Y_All(3:J+5,3:K+5)./(12.*h2) ...,
            .*(uh(1:J+3,3:K+5)-8.*uh(2:J+4,3:K+5)+8.*uh(4:J+6,3:K+5)-uh(5:J+7,3:K+5)) ).^2)) ...,
            +sum(sum(abs(X_All(3:J+5,3:K+5)./(12.*h2) ...,
            .*(uh_minus(3:J+5,1:K+3)-8.*uh_minus(3:J+5,2:K+4) ...,
            +8.*uh_minus(3:J+5,4:K+6)-uh_minus(3:J+5,5:K+7)) ...,
            -Y_All(3:J+5,3:K+5)./(12.*h2) ...,
            .*(uh_minus(1:J+3,3:K+5)-8.*uh_minus(2:J+4,3:K+5) ...,
            +8.*uh_minus(4:J+6,3:K+5)-uh_minus(5:J+7,3:K+5)) ).^2)) );
end

end

