% This File computes the policy functions and distributions Psi(m,theta_1) and Psi(m,theta_2), 
% Using EGM.

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
g0  = mgrid;                               %--initial guess for policy on future money-----@

% Q2 b: Finding policy functions using EGM

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
copt = zeros(Np,ns);
for i=1:ns;
    gopt(:,i) = go(:,:,i);
    copt(:,i) = mgrid./(1+gama) - gopt(:,i)+y-tau;
end;

%ploting the functions
figure(1);
title('Q2 a')

subplot(2,1,1);
plot(mgrid,[gopt,mgrid],'-o');
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
plot(mgrid,[copt,mgrid],'-o');
grid;
% xticks(mgrid);
% xtickformat('%,.2f')
% yticks(mgrid);
% ytickformat('%,.1f')
xlabel('$m$','Interpreter','latex','Fontsize',10); 
ylabel('$c$','Interpreter','latex','Fontsize',10); 
title('Policy Functions c(m, $\theta$), using EGM', 'Interpreter','latex');
legend('$g(m,\theta_1)$','$g(m,\theta_2)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);

