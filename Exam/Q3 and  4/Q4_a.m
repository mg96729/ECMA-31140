% This File computes the policy functions and distributions Psi(m,theta_1) and Psi(m,theta_2), 
% Using EGM.

clear all;
close all;
clc;

% Initialization

ns = 2;                                    %--number of shock states------@
QQ = [0.5,0.5;0.075,0.925];


ev = [1;0.1];                            %--endowment states----------------@
bet   = 0.99322;                              %--discount factor-------@
sig   = 1.5;                                 %--risk aversion---------@
q=1.0129;                        %--asset price


Np    = 1000;                                 %--number of points in the grid-----@

agrid = linspace(-2,2,Np)';            %--grid for approximation-----------@

Nk   = 1000;                               %--number of iterations (?)
g0  = agrid;                               %--initial guess for policy on future money-----@

%Finding policy functions using EGM

g = zeros(Np,Nk,ns);
for i=1:ns;                                %--initializing guess on policy function--------@
    g(:,1,i) = g0;
end;

for j=1:Nk-1;
    mc  = zeros(Np,ns);
    for i=1:ns;
        mc(:,i) = (agrid+ev(i)-q*g(:,j,i)).^(-sig);  %--this is u'(c')
    end;
    
    vf = (bet/q)*mc*QQ';
    
    ms  = zeros(Np,ns);
    for i=1:ns;    
        ms(:,i)  = vf(:,i).^(-1/sig)-ev(i)+q.*agrid;                        %--key step, finding m* the endogenous grid---@
    end;

    for i=1:Np;
        for s=1:ns;
            if agrid(i)<=ms(1,s);
            g(i,j+1,s) = agrid(1);
            elseif agrid(i)>=ms(Np,s);
            g(i,j+1,s) = agrid(Np);
            else;
            g(i,j+1,s) = LIP(ms(:,s),agrid,agrid(i));
            end;
        end;
    end;
    ac=abs(g(:,j+1,ns)-g(:,j,ns));
    if max(ac) <= 0.00000001;
    disp("convergence achieved");
    go = g(:,j+1,:);
    disp("number of iterations");
    disp(j);
    break;
    end;
end;

%Recovering policies*/

gopt = zeros(Np,ns);
%copt = zeros(Np,ns);
for i=1:ns;
    gopt(:,i) = go(:,:,i);
%    copt(:,i) = agrid + ev(i)-q.*gopt(:,i)
end;

%---Computing Distributions---%

g1 = gopt(:,1);
g2 = gopt(:,2);


Dp     = 2*Np-1;                                          %--number of points in the grid of F--@                                      
dgrid  = linspace(-2,2,Dp)';                         %--grid for distributions-------------@

Fini   = (dgrid-dgrid(1))   /(dgrid(Dp)-dgrid(1));         %--initial cummulative-----------@
%Fini = ones(Dp,1);

Prob      = QQ-eye(ns);
Prob(:,ns) = ones(ns,1);

a = zeros(ns,1);
a(ns) = 1;
epr   = linsolve(Prob',a);



Dk=68;                                %--maximum number iterations--------------@

F1 = zeros(Dp,Dk);
F1(:,1) = epr(1)*Fini;

F2 = zeros(Dp,Dk);
F2(:,1) = epr(2)*Fini;


for i=1:Dk-1;
    for j=1:Dp;

        if dgrid(j) < g1(1);
        F1(j,i+1) = QQ(1,1)*0+QQ(2,1)*LIP(dgrid,F2(:,i),LIP(g2,agrid,dgrid(j)));
        F2(j,i+1) = QQ(1,2)*0+QQ(2,2)*LIP(dgrid,F2(:,i),LIP(g2,agrid,dgrid(j))); 

        elseif dgrid(j) >= g1(1) & dgrid(j)<g2(Np);
        m1 = LIP(dgrid,F1(:,i),LIP(g1,agrid,dgrid(j)));
        m2 = LIP(dgrid,F2(:,i),LIP(g2,agrid,dgrid(j)));
        %Here is with Q update
        F1(j,i+1) = QQ(1,1)*m1+QQ(2,1)*m2;
        F2(j,i+1) = QQ(1,2)*m1+QQ(2,2)*m2; 

        elseif dgrid(j)>=g2(Np) & dgrid(j)<g1(Np);
        m2 = epr(2);
        m1 = LIP(dgrid,F1(:,i),LIP(g1,agrid,dgrid(j)));
        F1(j,i+1) = QQ(1,1)*m1+QQ(2,1)*m2;
        F2(j,i+1) = QQ(1,2)*m1+QQ(2,2)*m2;

        elseif dgrid(j)>=g1(Np);
        m1 = epr(1);
        m2 = epr(2);
        F1(j,i+1) = QQ(1,1)*m1+QQ(2,1)*m2;
        F2(j,i+1) = QQ(1,2)*m1+QQ(2,2)*m2;

        end;
    end;
    h=max(abs(F1(:,i+1)-F1(:,i)));

    if h<=0.00000000001;
    F1N = F1(:,i+1);
    F2N = F2(:,i+1);
    display('distribution convergence achieved');
    display(i);
    break;
    end;
end; 

% Plot
figure(2);
title('Distributions')

subplot(2,1,1);
plot(dgrid,F1());
grid;
% xticks(dgrid);
% xtickformat('%,.2f')
xlabel('$a$','Interpreter','latex','Fontsize',10); 
ylabel('$\Psi(a,e_1)$','Interpreter','latex','Fontsize',10); 
title('Finding $\Psi(a,e_1)$','Interpreter','latex');

subplot(2,1,2);
plot(dgrid,F2);
grid;
% xticks(dgrid);
% xtickformat('%,.2f')
xlabel('$m$','Interpreter','latex','Fontsize',10); 
ylabel('$\Psi(a,e_0)$','Interpreter','latex','Fontsize',10); 
title('Finding $\Psi(a,e_0)$','Interpreter','latex');
%legend('$\Psi_0(m,\theta_2)$','$\Psi_1(m,\theta_2)$','$\Psi_2(m,\theta_2)$','$\Psi_3(m,\theta_2)$','Interpreter','latex', 'Location', 'Southeast','Fontsize',10);
%exportgraphics(gcf,['D:\Dropbox\Chicago 2024\winter 24\ecma 3114\classnotes\cn 10\fig7.pdf'])

