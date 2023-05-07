clear;
close all;
%% Define Default Parameters

P_B = 325e6; % Rated power, VA
V_r = 20e3; % Rated voltage, V
V_B = V_r * sqrt(2)/sqrt(3);
I_B = P_B * (2/3) / V_B;
Z_B = V_B/I_B;
pf_B = 0.85; % Rated power factor
P = 64; % Number of poles
rpm = 112.5; % Rated mechanical speed, rpm
J = 35.1e6; % Inertia constant, J*s^2 = kg*m
r_s = .0019; % Stator winding resistance
X_ls = .1200; % Stator leakage reactance at rated freq
X_q = .4800; % Stator q-axis reactance at rated freq
X_d = .8500; % Stator d-axis reactance at rated freq
r_fd_p = .00041; % Field winding resistance (ref to stator side)
X_lfd_p = .2049; % Field winding leakage reactance at rated freq (ref to stator side)
r_kq1_p = Inf; % q-axis Damper winding 1 resistance at rated freq (ref to stator side)
r_kq2_p = .0136; % q-axis Damper winding 2 resistance at rated freq (ref to stator side)
r_kd_p = .0131; % d-axis Damper winding resistance at rated freq (ref to stator side)
X_lkq1_p = Inf; % q-axis Damper winding 1 leakage reactance at rated freq (ref to stator side)
X_lkq2_p = .1029; % q-axis Damper winding 2 leakage reactance at rated freq (ref to stator side)
X_lkd_p = .1600; % d-axis Damper winding leakage reactance at rated freq (ref to stator side)

v_in = 1;
e_in = 1;
T_I = -1;

X_mq = X_q - X_ls; % q-axis magnetizing reactance 
X_md = X_d - X_ls; % d-axis magnetizing reactance 

omega_b = rpm * (P/2) * 2*pi / 60; % electrical rad / s
T_B = P_B / (2/P * omega_b);

H = 0.5 * (2/P)^2 * J * omega_b ^2 / P_B;

X_aq = 1/(1/X_mq + 1/X_ls + 1/X_lkq1_p + 1/X_lkq2_p);
X_ad = 1/(1/X_md + 1/X_ls + 1/X_lfd_p + 1/X_lkd_p);

param = [
    omega_b
    omega_b
    r_s
    X_ls
    X_q
    X_d
    r_fd_p
    X_lfd_p
    r_kq1_p
    r_kq2_p
    r_kd_p
    X_lkq1_p
    X_lkq2_p
    X_lkd_p
    v_in
    e_in
    T_I
    H
];

%% Solve for Steady-State conditions
% Steady-state
delta = atan(X_q * 1 * pf_B / (3 * V_r/V_B) / (V_r/V_B));
v_qs_init = cos(delta);
v_ds_init = sin(delta);
vfd_init = e_in * r_fd_p / X_md;

A = [
    r_s X_d X_md
    -X_q r_s 0
    0 0 r_fd_p
    ];

i_init = A^-1 * [v_qs_init; v_ds_init; vfd_init];

i_qs_init = i_init(1);
i_ds_init = i_init(2);
i_fd_init = i_init(3);

psi_qs_init = X_ls * i_qs_init + X_mq * i_qs_init;
psi_ds_init = X_ls * i_ds_init + X_md * (i_ds_init+i_fd_init);
psi_0s_init = 0;

psi_fd_init = X_lfd_p * i_fd_init + X_md * (i_ds_init + i_fd_init);

psi_kq1_init = 0;
psi_kq2_init = X_mq * i_qs_init;
psi_kd_init = X_md * (i_ds_init + i_fd_init);

y0 = [
    psi_qs_init
    psi_ds_init
    psi_0s_init
    psi_kq1_init
    psi_kq2_init
    psi_fd_init
    psi_kd_init
    omega_b
    0
];

%% Equilibrate to steady state with given parameters
[t,y] = rk4(@(t,y) simulate(t,y,param,17,{@(t,y) 0 .* (t>4)}), [0,20], y0, 1e-4);
y0 = y(:,end);

figure;
plot(t,y(8,:)/omega_b,'linewidth',2);
format_plot('t','$\omega_r/\omega_b$');
%xline(4,'x--','linewidth',1.5)

%% Verification Simulations
datav = readmatrix("data_omega.csv");
[t,y] = rk4(@(t,y) simulate(t,y,param,17,{@(t,y) T_I * (t>=1)}), [0,10], y0, 1e-4);
%[tr,yr] = rk4(@(t,y) simulate(t,y,param,17,{@(t,y) (t-1)*T_I*(t<2&&t>=1) + T_I * (t>=2)}), [0,10], y0, 1e-4);

psi_qs_r = y(1,:);
psi_ds_r = y(2,:);
psi_0s = y(3,:);
psi_kq1_pr = y(4,:);
psi_kq2_pr = y(5,:);
psi_fd_pr = y(6,:);
psi_kd_pr = y(7,:);
omega_r = y(8,:);
theta_r = y(9,:);

psi_mq_r = X_aq * (psi_qs_r/X_ls + psi_kq1_pr/X_lkq1_p + psi_kq2_pr/X_lkq2_p);
psi_md_r = X_ad * (psi_ds_r/X_ls + psi_fd_pr/X_lfd_p + psi_kd_pr/X_lkd_p);

i_qs_r = (psi_qs_r - psi_mq_r) / X_ls;
i_ds_r = (psi_ds_r - psi_md_r) / X_ls;

i_fd_pr = (psi_fd_pr - psi_md_r) / X_lfd_p;

T_e = psi_ds_r .* i_qs_r - psi_qs_r .* i_ds_r;
delta = theta_r-omega_b.*t-theta_r(1);

figure;
hold on;
plot(t,omega_r./omega_b,'LineWidth',2)
plot(datav(:,1)-datav(min(find(datav(2:end,2)-datav(1:end-1,2) > .001)),1)+1,datav(:,2)-datav(1,2)+1,'LineWidth',2)
%plot(t,omega_r_r./omega_b,'LineWidth',2)
%legend(["Step Change","Ramp"],'box','off')

axis([0,10,0.99,1.01])
format_plot('t','$\frac{\omega_r}{\omega_b}$')

figure;
data2d = readmatrix("t.csv");
hold on;
plot(delta.*180/pi,-T_e*T_B/1e6,'linewidth',2);
plot(data2d(:,1),data2d(:,2),'o','markersize',3,'markerfacecolor',"#D95319",'linewidth',2)
axis([0,60,0,55.2])
format_plot("$\delta$ (electrical degrees)","$T_e$ ($10^6$ N$\cdot$m)")
legend(["Simulated","Reference"],'box','off','interpreter','latex')
yticks([0,round(T_B/1e6,1),round(2*T_B/1e6,1)])

figure;
hold on;
datat = readmatrix("data_Te.csv")
plot(t,-T_e*T_B/10^6,'linewidth',2)
plot(datat(:,1)-datat(min(find(datat(2:end,2)-datat(1:end-1,2) > .001))+1,1)+1,datat(:,2)-datat(1,2),'linewidth',2);
legend(["Simulated","Reference"],'box','off','interpreter','latex')
format_plot("t (s)","$T_e$ ($10^6$ N$\cdot$m)")
yticks([0,round(T_B/1e6,1),round(2*T_B/1e6,1)])

axis([0,6,0,55.2])

%% New problems
delta = .3;
deltar = 1*2*pi;

dt = 1;

ap = @(t) omega_b*(delta/2*(cos((t-1)/dt*pi)-1).*(t<1+2*dt&t>=1) + 1);
an = @(t) omega_b*(1 + sin(2*pi/dt*10*t)*delta.*(t>=1 & t<1+2*dt));
ap2 = @(t) omega_b*(delta/2*(sin((t-1)/dt*pi)).*(t<1+2*dt&t>=1) + 1);
ar = @(t) omega_b - deltar * (t-1).*(t>=1) + 2*deltar*(t-(dt+1)).*(t>=(dt+1))-deltar*(t-(2*dt+1)).*(t>=(2*dt+1));
we = ap2;
endt = dt*50;

[tp,yp] = rk4(@(t,y) simulate(t,y,param,[2],{@(t,y) we(t)}), [0,endt], y0, 1e-4);
[tp2,yp2] = rk4(@(t,y) simulate(t,y,param,[2,17],{@(t,y) we(t),@(t,y) param(17) + control(y(8)/we(t),1,50)}), [0,endt], y0, 1e-4);
[tp3,yp3] = rk4(@(t,y) simulate(t,y,param,[2,17],{@(t,y) we(t),@(t,y) max(min(param(17) + control(y(8)/we(t),1,50),0),-2)}), [0,endt], y0, 1e-4);
%%
ts = 0:.01:endt;
figure;
hold on;
plot(ts,we(ts)/omega_b,'k','linewidth',4)
plot(tp,yp(8,:)/omega_b,'color','#D95319','linewidth',2)
%plot(tp2,yp2(8,:)/omega_b,'color','#0072BD','linewidth',2)
%plot(tp2,yp3(8,:)/omega_b,'color','#77AC30','linewidth',2)
axis([0,endt,0,1.5])
format_plot("t","$\omega/\omega_b$")
xline(1,'k--','linewidth',1)
xline(1+2*dt,'k--','linewidth',1)
%legend(["$\omega_e$","$\omega_r$ (no control)","$\omega_r$ (P control)","$\omega_r$ (clamped P control)"],'box','on','interpreter','latex','fontsize',10)
legend(["$\omega_e$","$\omega_r$ (no control)"],'box','on','interpreter','latex','fontsize',10)
title(delta * 100 + "\% Fluctuation in $\omega_e$",'interpreter','latex')
%%
figure;
hold on;
plot(tp,ones(size(tp)),'color','#D95319','linewidth',2)
plot(tp2,-(param(17) + control(yp2(8,:)./we(tp2),1,50)),'color','#0072BD','linewidth',2)
cclamp = param(17) + control(yp3(8,:)./we(tp3),1,50);
cclamp(cclamp>0) = 0;
cclamp(cclamp<-2) = -2;
plot(tp3,-cclamp,'color','#77AC30','linewidth',2)
format_plot("t","$T_I/T_B$")
xline(1,'k--','linewidth',1)
xline(1+dt*2,'k--','linewidth',1)
axis([0,endt,-inf,inf])
legend(["No control","P control","Clamped P control"],'box','on','interpreter','latex','fontsize',10)
title(delta * 100 + "\% Fluctuation in $\omega_e$",'interpreter','latex')

%% Power Calculations

v_abcs = v_in .* [cos(tp .* we(tp)); cos(tp .* we(tp) - (2*pi/3)); cos(tp.*we(tp)+(2*pi/3))];

psi_qs_r = yp(1,:);
psi_ds_r = yp(2,:);
psi_0s = yp(3,:);
psi_kq1_pr = yp(4,:);
psi_kq2_pr = yp(5,:);
psi_fd_pr = yp(6,:);
psi_kd_pr = yp(7,:);
omega_r = yp(8,:);
theta_r = yp(9,:);

psi_mq_r = X_aq * (psi_qs_r/X_ls + psi_kq1_pr/X_lkq1_p + psi_kq2_pr/X_lkq2_p);
psi_md_r = X_ad * (psi_ds_r/X_ls + psi_fd_pr/X_lfd_p + psi_kd_pr/X_lkd_p);

K_s_r_i = zeros(numel(theta_r),3,3);

K_s_r_i(:,:,1) = [cos(theta_r);cos(theta_r-2*pi/3);cos(theta_r+2*pi/3)]';
K_s_r_i(:,:,2) = [sin(theta_r);sin(theta_r-2*pi/3);sin(theta_r+2*pi/3)]';
K_s_r_i(:,:,3) = [ones(size(theta_r));ones(size(theta_r));ones(size(theta_r))]';

i_qd0s_r = zeros(numel(theta_r),1,3);

i_qs_r = (psi_qs_r - psi_mq_r) / X_ls;
i_ds_r = (psi_ds_r - psi_md_r) / X_ls;
i_qd0s_r(:,1,:)  = [i_qs_r', i_ds_r', zeros(size(i_qs_r))'];

K_s_r_i = permute(K_s_r_i,[2,3,1]);
i_qd0s_r = permute(i_qd0s_r,[3,2,1]);

i_abcs = (squeeze(pagemtimes(K_s_r_i, i_qd0s_r)));
ptot = sum(i_abcs.*v_abcs,1);
figure;
hold on;
plot(tp,(i_abcs(1,:)*I_B).^2*r_s*Z_B/1e6,'linewidth',2)
plot(tp-10,(i_abcs(1,:)*I_B).^2 * r_s*Z_B/1e6,'linewidth',2)
axis([0,0.05,0,12])
format_plot('t','$i_{a}^2 \cdot r_s$ (MW)')
legend(["Steady-State","Stalled"],'interpreter','latex')

%% Function Definisions

% derivative function
% use ps, ts, and vs to specify parameter changes (ps) with values (vs) at

function dydt = simulate(t,y,param,ps,vs)
    arguments
        t double
        y double
        param double
        ps (:,1) double = []
        vs cell = {@(t,y) 0}
    end
    
    for i = 1:numel(ps)
        param(ps(i)) = vs{i}(t,y);    
    end
    dydt = om_der(t,y,param);
end
function dydt = om_der(t,y,param)
    % Def param :
    omega_b = param(1);
    omega_e = param(2);
    r_s = param(3); % Stator winding resistance
    X_ls = param(4); % Stator leakage reactance at rated freq
    X_q = param(5); % Stator q-axis reactance at rated freq
    X_d = param(6); % Stator d-axis reactance at rated freq
    r_fd_p = param(7); % Field winding resistance (ref to stator side)
    X_lfd_p = param(8); % Field winding leakage reactance at rated freq (ref to stator side)
    r_kq1_p = param(9); % q-axis Damper winding 1 resistance at rated freq (ref to stator side)
    r_kq2_p = param(10); % q-axis Damper winding 2 resistance at rated freq (ref to stator side)
    r_kd_p = param(11); % d-axis Damper winding resistance at rated freq (ref to stator side)
    X_lkq1_p = param(12); % q-axis Damper winding 1 leakage reactance at rated freq (ref to stator side)
    X_lkq2_p = param(13); % q-axis Damper winding 2 leakage reactance at rated freq (ref to stator side)
    X_lkd_p = param(14); % d-axis Damper winding leakage reactance at rated freq (ref to stator side)
    
    v_in = param(15); % input stator voltages
    e_xfd_pr = param(16); % excitation
    
    T_I = param(17);
    
    H = param(18);
        
    X_mq = X_q - X_ls; % q-axis magnetizing reactance 
    X_md = X_d - X_ls; % d-axis magnetizing reactance 
    
    % Def y :
    psi_qs_r = y(1);
    psi_ds_r = y(2);
    psi_0s = y(3);
    psi_kq1_pr = y(4);
    psi_kq2_pr = y(5);
    psi_fd_pr = y(6);
    psi_kd_pr = y(7);
    omega_r = y(8);
    theta_r = y(9);
    
    
    % Transform to rotor reference frame
    K_s_r = (2/3) * [
        cos(theta_r) cos(theta_r - (2*pi/3)) cos(theta_r + (2*pi/3))
        sin(theta_r) sin(theta_r - (2*pi/3)) sin(theta_r + (2*pi/3))
        0.5 0.5 0.5
    ];
    
    v_abcs = v_in * [cos(t * omega_e); cos(t * omega_e - (2*pi/3)); cos(t*omega_e+(2*pi/3))];
    v_qd0s_r = (K_s_r * v_abcs);
    
    v_qs_r = v_qd0s_r(1);
    v_ds_r = v_qd0s_r(2);
    v_0s = v_qd0s_r(3);
    
    X_aq = 1/(1/X_mq + 1/X_ls + 1/X_lkq1_p + 1/X_lkq2_p);
    X_ad = 1/(1/X_md + 1/X_ls + 1/X_lfd_p + 1/X_lkd_p);

    psi_mq_r = X_aq * (psi_qs_r/X_ls + psi_kq1_pr/X_lkq1_p + psi_kq2_pr/X_lkq2_p);
    psi_md_r = X_ad * (psi_ds_r/X_ls + psi_fd_pr/X_lfd_p + psi_kd_pr/X_lkd_p);
    
    i_qs_r = (psi_qs_r - psi_mq_r) / X_ls;
    i_ds_r = (psi_ds_r - psi_md_r) / X_ls;
    
    i_kq1_pr = (psi_kq1_pr - psi_mq_r) / X_lkq1_p;
    i_kq2_pr = (psi_kq2_pr - psi_mq_r) / X_lkq1_p;
    i_kd_pr = (psi_kd_pr - psi_md_r) / X_lkd_p;
    i_fd_pr = (psi_fd_pr - psi_md_r) / X_lfd_p;

    
    T_e = psi_ds_r * i_qs_r - psi_qs_r * i_ds_r;
    
    v_kq1_pr = 0;
    v_kq2_pr = 0;
    v_kd_pr = 0;

    dydt = [
        omega_b * (v_qs_r - omega_r/omega_b * psi_ds_r + r_s / X_ls * (psi_mq_r - psi_qs_r))
        omega_b * (v_ds_r + omega_r/omega_b * psi_qs_r + r_s / X_ls * (psi_md_r - psi_ds_r))
        omega_b * (v_0s - r_s / X_ls * psi_0s)
        0 %omega_b * (v_kq1_pr + r_kq1_p / X_lkq1_p * (psi_mq_r - psi_kq1_r))
        omega_b * (v_kq2_pr + r_kq2_p / X_lkq2_p * (psi_mq_r - psi_kq2_pr))
        omega_b * (r_fd_p/X_md * e_xfd_pr + r_fd_p / X_lfd_p * (psi_md_r - psi_fd_pr))
        omega_b * (v_kd_pr + r_kd_p / X_lkd_p * (psi_md_r - psi_kd_pr))
        omega_b / (2*H) * (T_e - T_I)
        omega_r
    ];

end

% RK4 solver
function [t,y] = rk4(odefun, tspan, y0, h)
    t = tspan(1):h:tspan(end);
    y = zeros(numel(y0),numel(t));
    y(:,1) = y0;
    
    for i = 1:numel(t)-1
        k1 = h*odefun(t(i),y(:,i));
        k2 = h*odefun(t(i) + h/2,y(:,i)+k1/2);
        k3 = h*odefun(t(i) + h/2,y(:,i)+k2/2);
        k4 = h*odefun(t(i) + h,y(:,i)+k3);
        y(:,i+1) = y(:,i) + k1/6 + 1/3 * (k2+k3) + k4/6;
    end
end

function a = control(val,sp,k)
    err = val-sp;
    a = k * err;
end

% Plot formatter
function format_plot(xl,yl)
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