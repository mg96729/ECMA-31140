clear all;
close all;
clc;

% Initialization

ns = 2;                                    %--number of shock states------@
QQ = [0.5,0.5;0.075,0.925];


ev = [1;0.1];                            %--states of endowment----------------@
bet   = 0.99322;                              %--discount factor-------@
sig   = 1.5;                                 %--risk aversion---------@
q=1.0129;                        %--asset price

a_low = -2; %--borrowing constraint
Np    = 1000;                                 %--number of points in the grid-----@
agrid = linspace(a_low,2,Np)';            %--grid for approximation, over values of a-----------@

g  = ones(Np,ns);                         %--initial guess for policy on future money, by index-----@

%initialize error and value function
error = 1;
v = zeros(Np,ns);
v_old = zeros(Np,ns);
% while error is bigger than 1e-6, and there are more iterations left:
for iter = 1:1000;
    if error > 1e-6;
        %for each state values a and e
        for i=1:Np;
            for j=1:ns;
                a = agrid(i);
                e = ev(j);
                %possible c and u given m'
                max_c = e+a-q*a_low;
                c = max(min(a+e-q .* agrid, max_c),0);
                u = util(c,sig);
                v_rough = u + bet .* (v * QQ(j,:)');
                %Now update
                max_v= max(v_rough);
                v(i,j) = max_v;
                g(i,j) = find(v_rough==max_v,1,'first');
             end
        end

        %update error in v estimation
        error = norm(v(:,1)-v_old(:,1))*(1/Np);
        v_old = v;        
    else
        break;
    end
end

%Recovering policies*/

gopt = zeros(Np,ns);
for i=1:Np
    for j=1:ns
        gopt(i,j) = agrid(g(i,j));
    end
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
title('Policy Functions g(a, e), using VFI', 'Interpreter','latex');
legend('$g(a,e_1)$','$g(a,e_0)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);

subplot(2,1,2);
plot(agrid, v);
grid;
% xticks(mgrid);
% xtickformat('%,.2f')
% yticks(mgrid);
% ytickformat('%,.1f')
xlabel('$a$','Interpreter','latex','Fontsize',10); 
ylabel('$v$','Interpreter','latex','Fontsize',10); 
title('Value Functions v(a, e), using VFI', 'Interpreter','latex');
legend('$v(a,e_1)$','$v(a,e_0)$','Interpreter','latex', 'Location', 'Northwest','Fontsize',10);




function u = util(c,sigma)
    if sigma == 1;
        u = -log(c);
    else
        u = (c.^(1-sigma))./(1-sigma);
    end
    %return u;
end