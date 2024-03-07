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

%ploting the functions
figure(1);

subplot(2,1,1);
plot(agrid,[gopt,agrid]);
grid;
% xticks(mgrid);
% xtickformat('%,.2f')
% yticks(mgrid);
% ytickformat('%,.1f')
xlabel('$a$','Interpreter','latex','Fontsize',10); 
ylabel('$a^{\prime}$','Interpreter','latex','Fontsize',10); 
title('Policy Functions g(a, e), using EGM', 'Interpreter','latex');
legend('$g(a,e1)$','$g(a,e0)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);

