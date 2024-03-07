% This File computes the policy functions and distributions Psi(m,theta_1) and Psi(m,theta_2), 
% Using VFI.

clear all;
close all;
clc;

% Initialization

ns = 2;                                    %--number of shock states------@
QQ = [0.73,0.27;0.27,0.73];


theta = [0.74;1.36];                            %--shocks----------------@
bet   = 0.98;                              %--discount factor-------@
sig   = 2;                                 %--risk aversion---------@
y     = 1;                                 %--endowment-------------@
gama  = 0.02;                            %--inflation rate--------@
msup  = 1;                                 %--money supply----------@
tau   = 0.0234;                %--tau supply------------@

Np    = 1000;                                 %--number of points in the grid-----@

mlow  = y+tau;                             %--lower value of the grid----------@
mup   = 2.5;                                 %--upper value of the grid----------@
mgrid = linspace(mlow,mup,Np)';            %--grid for approximation-----------@

Nk   = 1000;                               %--number of iterations (?)
g0  = ones(Np,ns);                               %--initial guess for policy on future money, by index-----@
v0  = zeros(Np,ns);                        %--initial guess for value function

% Q2 a: Finding policy functions using VFI

%--initializing guess on policy function, by index--------@
g = g0;

%initialize error and value function
error = 1;
v = v0;
v_old = ones(Np,ns);
% while error is bigger than 1e-6, and there are more iterations left:
for iter = 1:1000;
    if error > 1e-6;
        %for each state value m
        for i=1:Np;
            m = mgrid(i);
            %possible c and u given m'
            max_c = m/(1+gama);
            c = max(min(y + tau + m/(1+gama) - mgrid, max_c),0);
            u = util(c,sig);
            v_rough = u * theta'+ bet .* (v * QQ');
            %Now update for each state s
            for j=1:ns
                max_v= max(v_rough(:,j));
                v(i,j) = max_v;
                g(i,j) = find(v_rough(:,j)==max_v,1,'last');
            end
        end;

        %update error in v estimation
        error = max(max(abs(v-v_old)));
        v_old = v;        
    else;
        break;
    end;
 end;

%Recovering policies*/

gopt = zeros(Np,ns);
copt = zeros(Np,ns);
for i=1:Np
    for j=1:ns
        gopt(i,j) = mgrid(g(i,j));
        copt(i,j) = mgrid(i)/(1+gama) - gopt(i,j)+y-tau;
    end
end;

%ploting the functions
figure(1);
title('Q2 a')

subplot(2,1,1);
plot(mgrid,[gopt,mgrid]);
grid;
% xticks(mgrid);
% xtickformat('%,.2f')
% yticks(mgrid);
% ytickformat('%,.1f')
xlabel('$m$','Interpreter','latex','Fontsize',10); 
ylabel('$m^{\prime}$','Interpreter','latex','Fontsize',10); 
title('Policy Functions g(m, $\theta$), using EGM', 'Interpreter','latex');
legend('$g(m,\theta_1)$','$g(m,\theta_2)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);

subplot(2,1,2);
plot(mgrid,[copt,mgrid]);
grid;
% xticks(mgrid);
% xtickformat('%,.2f')
% yticks(mgrid);
% ytickformat('%,.1f')
xlabel('$m$','Interpreter','latex','Fontsize',10); 
ylabel('$c$','Interpreter','latex','Fontsize',10); 
title('Policy Functions c(m, $\theta$), using EGM', 'Interpreter','latex');
legend('$g(m,\theta_1)$','$g(m,\theta_2)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);



function u = util(c,sigma)
    if sigma == 1;
        u = -log(c);
    else
        u = (c.^(1-sigma))./(1-sigma);
    end
    %return u;
end

% Recovering policies*/
% 
% gopt = zeros(Np,ns);
% copt = zeros(Np,ns);
% for i=1:ns;
%     gopt(:,i) = go(:,:,i);
%     copt(:,i) = mgrid./(1+gama) - gopt(:,i)+y-tau;
% end;
% 
% 
% figure(1);
% title('Q2 b')
% 
% subplot(2,1,1);
% plot(mgrid,[gopt,mgrid],'-o');
% grid;
% xticks(mgrid);
% xtickformat('%,.2f')
% yticks(mgrid);
% ytickformat('%,.1f')
% xlabel('$m$','Interpreter','latex','Fontsize',10); 
% ylabel('$m^{\prime}$','Interpreter','latex','Fontsize',10); 
% title('Policy Functions g(m, $\theta$), using EGM', 'Interpreter','latex');
% legend('$g(m,\theta_1)$','$g(m,\theta_2)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);
% 
% subplot(2,1,2);
% plot(mgrid,[copt,mgrid],'-o');
% grid;
% xticks(mgrid);
% xtickformat('%,.2f')
% yticks(mgrid);
% ytickformat('%,.1f')
% xlabel('$m$','Interpreter','latex','Fontsize',10); 
% ylabel('$c$','Interpreter','latex','Fontsize',10); 
% title('Policy Functions c(m, $\theta$), using EGM', 'Interpreter','latex');
% legend('$g(m,\theta_1)$','$g(m,\theta_2)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);
% 
% ---Computing Distributions---%
% 
% g1 = gopt(:,1);
% g2 = gopt(:,2);
% 
% 
% Dp     = 2*Np-1;                                          %--number of points in the grid of F--@                                      
% dgrid  = linspace(mlow,mup,Dp)';                         %--grid for distributions-------------@
% 
% Fini   = (dgrid-dgrid(1))   /(dgrid(Dp)-dgrid(1));         %--initial cummulative-----------@
% Fini = ones(Dp,1);
% 
% Prob      = QQ-eye(ns);
% Prob(:,ns) = ones(ns,1);
% 
% a = zeros(ns,1);
% a(ns) = 1;
% epr   = linsolve(Prob',a);
% 
% 
% 
% Dk=68;                                %--maximum number iterations--------------@
% 
% F1 = zeros(Dp,Dk);
% F1(:,1) = epr(1)*Fini;
% 
% F2 = zeros(Dp,Dk);
% F2(:,1) = epr(2)*Fini;
% 
% 
% for i=1:Dk-1;
%     for j=1:Dp;
% 
%         if dgrid(j) < g1(1);
%         F1(j,i+1) = QQ(1,1)*0+QQ(2,1)*LIP(dgrid,F2(:,i),LIP(g2,mgrid,dgrid(j)));
%         F2(j,i+1) = QQ(1,2)*0+QQ(2,2)*LIP(dgrid,F2(:,i),LIP(g2,mgrid,dgrid(j))); 
% 
%         elseif dgrid(j) >= g1(1) & dgrid(j)<g2(Np);
%         m1 = LIP(dgrid,F1(:,i),LIP(g1,mgrid,dgrid(j)));
%         m2 = LIP(dgrid,F2(:,i),LIP(g2,mgrid,dgrid(j)));
%         F1(j,i+1) = QQ(1,1)*m1+QQ(2,1)*m2;
%         F2(j,i+1) = QQ(1,2)*m1+QQ(2,2)*m2; 
% 
%         elseif dgrid(j)>=g2(Np) & dgrid(j)<g1(Np);
%         m2 = epr(2);
%         m1 = LIP(dgrid,F1(:,i),LIP(g1,mgrid,dgrid(j)));
%         F1(j,i+1) = QQ(1,1)*m1+QQ(2,1)*m2;
%         F2(j,i+1) = QQ(1,2)*m1+QQ(2,2)*m2;
% 
%         elseif dgrid(j)>=g1(Np);
%         m1 = epr(1);
%         m2 = epr(2);
%         F1(j,i+1) = QQ(1,1)*m1+QQ(2,1)*m2;
%         F2(j,i+1) = QQ(1,2)*m1+QQ(2,2)*m2;
% 
%         end;
%     end;
%     h=max(abs(F1(:,i+1)-F1(:,i)));
% 
%     if h<=0.00000000001;
%     F1N = F1(:,i+1);
%     F2N = F2(:,i+1);
%     display('distribution convergence achieved');
%     display(i);
%     break;
%     end;
% end;
% 
% 
% 
% figure(2);
% subplot(2,1,1);
% 
% plot(mgrid,[gopt,mgrid],'-o');
% grid;
% xticks(mgrid);
% xtickformat('%,.2f')
% yticks(mgrid);
% ytickformat('%,.1f')
% xlabel('$m$','Interpreter','latex','Fontsize',10); 
% ylabel('$m^{\prime}$','Interpreter','latex','Fontsize',10); 
% title('Policy Functions');
% legend('$g(m,\theta_1)$','$g(m,\theta_2)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);
% 
% subplot(2,1,2);
% plot(dgrid,F2,'-o');
% grid;
% xticks(dgrid);
% xtickformat('%,.2f')
% xlabel('$m$','Interpreter','latex','Fontsize',10); 
% ylabel('$\Psi(m,\theta_2)$','Interpreter','latex','Fontsize',10); 
% title('Finding $\Psi(m,\theta_2)$','Interpreter','latex');
% legend('$\Psi_0(m,\theta_2)$','$\Psi_1(m,\theta_2)$','$\Psi_2(m,\theta_2)$','$\Psi_3(m,\theta_2)$','Interpreter','latex', 'Location', 'Southeast','Fontsize',10);
% exportgraphics(gcf,['D:\Dropbox\Chicago 2024\winter 24\ecma 3114\classnotes\cn 10\fig7.pdf'])
