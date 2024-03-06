% This File computes the policy functions and distributions Psi(m,theta_1) and Psi(m,theta_2).
% Sergio Salas

clear all;
close all;
clc;

ns = 2;                                    %--number of shocks------@
[z,QQ] = tauchen(ns,0,0.2,0.2042,3);


theta = exp(z);                            %--shocks----------------@
bet   = 0.95;                              %--discount factor-------@
sig   = 2;                                 %--risk aversion---------@
y     = 1;                                 %--endowment-------------@
gama  = 0.0526;                            %--inflation rate--------@
msup  = 1;                                 %--money supply----------@
tau   = msup*gama/(1+gama);                %--tau supply------------@

Np    = 8;                                 %--number of points in the grid-----@

mlow  = y+tau;                             %--lower value of the grid----------@
mup   = 2.5;                                 %--upper value of the grid----------@
mgrid = linspace(mlow,mup,Np)';            %--grid for approximation-----------@

Nk   = 1000;
g0  = mgrid;                               %--initial guess for policy on future money-----@

g = zeros(Np,Nk,ns);
for i=1:ns;                                %--initializing guess on policy function--------@
    g(:,1,i) = g0;
end;

for j=1:Nk-1;

    mc  = zeros(Np,ns);
    for i=1:ns;
        mc(:,i) = theta(i)*(mgrid/(1+gama)+y+tau-g(:,j,i)).^(-sig);  %--this is theta'u'(c')
    end;
    
    vf = bet/(1+gama)*mc*QQ';
    
    ms  = zeros(Np,ns);
    for i=1:ns;    
        ms(:,i)  = ((vf(:,i)/theta(i)).^(-1/sig)+mgrid-y-tau)*(1+gama);                        %--key step, finding m* the endogenous grid---@
    end;

    for i=1:Np;
        for s=1:ns;
            if mgrid(i)<=ms(1,s);
            g(i,j+1,s) = mgrid(1);
            elseif mgrid(i)>=ms(Np,s);
            g(i,j+1,s) = mgrid(Np);
            else;
            g(i,j+1,s) = LIP(ms(:,s),mgrid,mgrid(i));
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
for i=1:ns;
    gopt(:,i) = go(:,:,i);
end;



%---Computing Distributions---%

g1 = gopt(:,1);
g2 = gopt(:,2);


Dp     = 2*Np-1;                                          %--number of points in the grid of F--@                                      
dgrid  = linspace(mlow,mup,Dp)';                         %--grid for distributions-------------@

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
        F1(j,i+1) = QQ(1,1)*0+QQ(2,1)*LIP(dgrid,F2(:,i),LIP(g2,mgrid,dgrid(j)));
        F2(j,i+1) = QQ(1,2)*0+QQ(2,2)*LIP(dgrid,F2(:,i),LIP(g2,mgrid,dgrid(j))); 

        elseif dgrid(j) >= g1(1) & dgrid(j)<g2(Np);
        m1 = LIP(dgrid,F1(:,i),LIP(g1,mgrid,dgrid(j)));
        m2 = LIP(dgrid,F2(:,i),LIP(g2,mgrid,dgrid(j)));
        F1(j,i+1) = QQ(1,1)*m1+QQ(2,1)*m2;
        F2(j,i+1) = QQ(1,2)*m1+QQ(2,2)*m2; 

        elseif dgrid(j)>=g2(Np) & dgrid(j)<g1(Np);
        m2 = epr(2);
        m1 = LIP(dgrid,F1(:,i),LIP(g1,mgrid,dgrid(j)));
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


figure(1);
subplot(2,1,1);

plot(mgrid,[gopt,mgrid],'-o');
grid;
xticks(mgrid);
xtickformat('%,.2f')
yticks(mgrid);
ytickformat('%,.1f')
xlabel('$m$','Interpreter','latex','Fontsize',10); 
ylabel('$m^{\prime}$','Interpreter','latex','Fontsize',10); 
title('Policy Functions');
legend('$g(m,\theta_1)$','$g(m,\theta_2)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);

subplot(2,1,2);
plot(dgrid,F1,'-o');
grid;
xticks(dgrid);
xtickformat('%,.2f')
xlabel('$m$','Interpreter','latex','Fontsize',10); 
ylabel('$\Psi(m,\theta_1)$','Interpreter','latex','Fontsize',10); 
title('Finding $\Psi(m,\theta_1)$','Interpreter','latex');
%legend('$\Psi_0(m,\theta_1)$','$\Psi_1(m,\theta_1)$','$\Psi_2(m,\theta_1)$','$\Psi_3(m,\theta_1)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);
%exportgraphics(gcf,['D:\Dropbox\Chicago 2024\winter 24\ecma 3114\classnotes\cn 10\fig6.pdf'])


figure(2);
subplot(2,1,1);

plot(mgrid,[gopt,mgrid],'-o');
grid;
xticks(mgrid);
xtickformat('%,.2f')
yticks(mgrid);
ytickformat('%,.1f')
xlabel('$m$','Interpreter','latex','Fontsize',10); 
ylabel('$m^{\prime}$','Interpreter','latex','Fontsize',10); 
title('Policy Functions');
legend('$g(m,\theta_1)$','$g(m,\theta_2)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);

subplot(2,1,2);
plot(dgrid,F2,'-o');
grid;
xticks(dgrid);
xtickformat('%,.2f')
xlabel('$m$','Interpreter','latex','Fontsize',10); 
ylabel('$\Psi(m,\theta_2)$','Interpreter','latex','Fontsize',10); 
title('Finding $\Psi(m,\theta_2)$','Interpreter','latex');
%legend('$\Psi_0(m,\theta_2)$','$\Psi_1(m,\theta_2)$','$\Psi_2(m,\theta_2)$','$\Psi_3(m,\theta_2)$','Interpreter','latex', 'Location', 'Southeast','Fontsize',10);
%exportgraphics(gcf,['D:\Dropbox\Chicago 2024\winter 24\ecma 3114\classnotes\cn 10\fig7.pdf'])
