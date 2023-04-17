clear;
close all;
data = readmatrix("ENGN1140_lab_1_data_2023/3_23_23data.txt")

t = data(:,1);
y = data(:,3);
delta = data(:,2);

istart = min(find(delta~=0));
iend = sum(t<725);

X = [t(istart:iend)-t(istart),delta(istart:iend)];
figure;
hold on;
plot(t,y,'o','MarkerFaceColor',[0  0.4470 0.7410],'markersize',2);


testf = @(b,x) (abs(x(:,1)-b(2))+x(:,1)-b(2))./(2*(x(:,1)-b(2))) * b(1).*x(1,2).*(1-exp(-(x(:,1)-b(2))./b(3)))+y(istart);
mdl = fitnlm(X,y(istart:iend),testf,[1,1,1,y(istart)]);


plot(t(istart:iend),testf(mdl.Coefficients.Estimate,X),'linewidth',2);

param = mdl.Coefficients.Estimate;
k_fit = param(1)
theta_fit = param(2)
tau_fit = param(3)
y0_fit = param(4)

[ti,yi] = rk2(@(it,iy) -(iy)/tau_fit+k_fit/tau_fit*interp1(t,delta,max(it-theta_fit,0)),[0,max(t)],0,.1);
[tt,yt] = rk2(@(it,iy) -(iy)/120+0.5/120*interp1(t,delta,max(it-10,0)),[0,max(t)],0,.1);

plot(ti,yi+y0_fit,'linewidth',2)
plot(t,delta)


%% Plot regular lab data
figure;
hold on;
plot(data(:,1),data(:,4),'-','linewidth',2)
plot(data(:,1),data(:,5),'-','linewidth',2)
plot(ti,yi+y0_fit,'linewidth',2)
%plot(tt,yt+y0_fit,'linewidth',2)

plot(data(:,1),data(:,3),'k-','linewidth',2)

%plot(x,y_fit,'linewidth',2)

format("t (s)","T ($^\circ$ C)")
legend(["First Principles","FOPDT (assumed)","FOPD (fit)","Measured"],'box','off');

function format(xl,yl)
    xlabel(xl,'Interpreter','Latex','FontSize',20)
    ylabel(yl,'Interpreter','Latex','FontSize',20)
    
    box on;
    set(gcf,'color','w')
    set(gca,'FontName','Times','FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    xa = gca;
    xa.TickLength = [0.025,0.025];
    xa.LineWidth = 1.5;
end

function [t,y] = rk2(odefun, tspan, y0, h)
% Custom Implementation of rk2

    t = tspan(1):h:tspan(end);
    y = zeros(numel(y0),numel(t));
    y(:,1) = y0;
    
    alpha = 0.5;
    beta = alpha;
    gamma1 = 1-1/(2*alpha);
    gamma2 = 1/(2*alpha);
    for i = 1:numel(t)-1
        k1 = h * odefun(t(i),y(:,i));
        k2 = h * odefun(t(i) + alpha * h,y(:,i) + beta * k1);
        y(:,i+1) = y(:,i) + gamma1 * k1 + gamma2 * k2;
    end
end
