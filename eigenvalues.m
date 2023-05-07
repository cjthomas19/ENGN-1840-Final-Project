
syms omega_b v_ds_r v_qs_r v_0s v_kq2_pr e_xfd_pr v_kd_pr H

syms r_s X_ls X_q X_d r_fd_p r_lfd_p r_kq1_p r_kq2_p r_kd_p X_lfd_p X_lkq1_p X_lkq2_p X_lkd_p T_I v_in e_in

syms psi_qs_r psi_ds_r psi_0s psi_kq1_pr psi_kq2_pr psi_fd_pr psi_kd_pr omega_r theta_r

ys = [psi_qs_r psi_ds_r psi_0s psi_kq1_pr psi_kq2_pr psi_fd_pr psi_kd_pr omega_r theta_r];

param = [
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

X_mq = X_q - X_ls;
X_md = X_d - X_ls;

X_aq = 1/(1/X_mq + 1/X_ls + 1/X_lkq1_p + 1/X_lkq2_p);
X_ad = 1/(1/X_md + 1/X_ls + 1/X_lfd_p + 1/X_lkd_p);

psi_mq_r = X_aq * (psi_qs_r/X_ls + psi_kq1_pr/X_lkq1_p + psi_kq2_pr/X_lkq2_p);
psi_md_r = X_ad * (psi_ds_r/X_ls + psi_fd_pr/X_lfd_p + psi_kd_pr/X_lkd_p);

i_qs_r = (psi_qs_r - psi_mq_r) / X_ls;
i_ds_r = (psi_ds_r - psi_md_r) / X_ls;

T_e = psi_ds_r * i_qs_r - psi_qs_r * i_ds_r;



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
rpm = 112.5; % Rated mechanical speed, rpm
J = 35.1e6; % Inertia constant, J*s^2 = kg*m
r_sa = .0019; % Stator winding resistance
X_lsa = .1200; % Stator leakage reactance at rated freq
X_qa = .4800; % Stator q-axis reactance at rated freq
X_da = .8500; % Stator d-axis reactance at rated freq
r_fd_pa = .00041; % Field winding resistance (ref to stator side)
X_lfd_pa = .2049; % Field winding leakage reactance at rated freq (ref to stator side)
r_kq1_pa = Inf; % q-axis Damper winding 1 resistance at rated freq (ref to stator side)
r_kq2_pa = .0136; % q-axis Damper winding 2 resistance at rated freq (ref to stator side)
r_kd_pa = .0131; % d-axis Damper winding resistance at rated freq (ref to stator side)
X_lkq1_pa = Inf; % q-axis Damper winding 1 leakage reactance at rated freq (ref to stator side)
X_lkq2_pa = .1029; % q-axis Damper winding 2 leakage reactance at rated freq (ref to stator side)
X_lkd_pa = .1600; % d-axis Damper winding leakage reactance at rated freq (ref to stator side)
P = 64; % Number of poles
P_B = 325e6; % Rated power, VA

v_ina = 1;
e_ina = 1;
T_Ia = -1;

omega_ba = rpm * (P/2) * 2*pi / 60; % electrical rad / s
Ha = 0.5 * (2/P)^2 * J * omega_ba ^2 / P_B;

parama = [
    omega_ba
    r_sa
    X_lsa
    X_qa
    X_da
    r_fd_pa
    X_lfd_pa
    r_kq1_pa
    r_kq2_pa
    r_kd_pa
    X_lkq1_pa
    X_lkq2_pa
    X_lkd_pa
    v_ina
    e_ina
    T_Ia
    Ha
];
%%
dydt = subs(dydt,param,parama);

j = jacobian(dydt,ys);

%%
me = zeros(1,floor(numel(y(1,:))/100));
e = eig(j);
for i = 1:numel(me)
    se = double(subs(e,ys,y(:,i*100)'));
    me(i) = se(min(find(abs(se) == max(abs(se)))));
    disp(i);
end

%%
figure;
hold on;
%%
plot(real(me(end)).*[1e-4 1e-3 1e-2],imag(me(end)).*[1e-4 1e-3 1e-2],'k--','linewidth',2)
plot(real(me(end))*1e-4,imag(me(end))*1e-4,'.','color',"#0072BD",'markersize',30,'linewidth',2)
plot(real(me(end))*1e-3,imag(me(end))*1e-3,'.','color',"#D95319",'markersize',30,'linewidth',2)
plot(real(me(end))*1e-2,imag(me(end))*1e-2,'.','color',"#EDB120",'markersize',30,'linewidth',2)

%%
figure;
plot(real(me),imag(me),'linewidth',2);
format_plot("$\lambda_R$","$\lambda_I$")
%% 
draw_stability([1/24 1/6 1/2 1 1])


function draw_stability(coeff)
%Generate stability diagrams. Takes in coefficients of sigma (amplification
%factor)
%Solves e^i theta = sigma for different values of theta
    th = linspace(0,2*pi,1000);
    n = numel(coeff);
    pts = zeros(numel(th)*(n-1),1);
    for ii = 1:numel(th)
        t = th(ii);
        a = (n-1)*(ii-1) + 1;
        pts(a:a+(n-2)) = roots([coeff(1:numel(coeff)-1), coeff(numel(coeff))+exp(1i*t)]);
    end
    hold on;
    %Plot points and format
    %plot([-4,4],[0,0],'k','LineWidth',.5)
    %plot([0,0],[-4,4],'k','LineWidth',.5)
    plot(real(pts), imag(pts),'k.','MarkerSize',4)
    axis square;
    axis([-5,5,-5,5])

    xlabel('$\lambda_R h$','Interpreter','Latex','FontSize',20)
    ylabel('$\lambda_I h$','Interpreter','Latex','FontSize',20)
    
    box on;
    set(gcf,'color','w')
    set(gca,'FontName','Times','FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    xa = gca;
    xa.TickLength = [0.025,0.025];
    xa.LineWidth = 1.5;
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