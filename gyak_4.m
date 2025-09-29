% JR-2, 4. gyakorlat (09.26.)

%% Skálázási tétel

% Azonos energiájú jelek!
% fesz/áram * idő területek egyenlőek: A_{i}*tau_{i} = const.

clc;
clear all;

A = [2,3,4,10,100];
T = 1;
w_0 = 2*pi/T;
tau = T./A;
w = -12*pi/T : 0.1*pi/T : 12*pi/T;
F = zeros([length(w),length(tau)]);
G = zeros([length(w),length(tau)]);
for d = 1 : length(tau)
    for k = 1 : length(w)
        if w(k) ~= 0
            F(k,d) = A(d)*tau(d)*sin( w(k)*tau(d)/2 )/( w(k)*tau(d)/2 );
            G(k,d) = (F(k,d))^2;
        else
            F(k,d) = A(d)*tau(d);
            G(k,d) = (F(k,d))^2;
        end
    end
end

size = 12;
figure(1)
colororder({'r','g','b','m','k'})
title('Skálázási tétel');
set(gca,'fontsize',size)
xlabel('\omega','FontSize',size);
ylabel('F(j\omega)','FontSize',size);
grid on;
legend('Location','northeast','NumColumns',1,'Interpreter','latex','FontSize',size);
hold on
for d = 1 : length(tau)
    plot(w,F(:,d),...
        'DisplayName',['$\Delta$: $\tau$ = T/',num2str(A(d))],...
        'LineStyle','-',...
        'Marker','none',...
        'LineWidth',2);
end
for d = 1 : length(tau)
    plot(w,G(:,d),...
        'DisplayName',['$\sqcap$: $\tau$ = T/',num2str(A(d))],...
        'LineStyle','-.',...
        'Marker','none',...
        'LineWidth',2);
end
hold off

t = -1*T : 0.005*T : 1*T;
f = zeros(length(tau),length(t));
g = zeros(length(tau),length(t));
for d = 1 : length(tau)
    f(d,:) = A(d)*tau(d);
    g(d,:) = (f(d,:)).^2;
    for p = 1 : 1000
        F_p = A(d)*tau(d)*sin( p*w_0*tau(d)/2 )/( p*w_0*tau(d)/2 );
        G_p = (F_p)^2;
        f(d,:) = f(d,:) + F_p*exp(1i*p*w_0*t) + F_p*exp(1i*(-p)*w_0*t);
        g(d,:) = g(d,:) + G_p*exp(1i*p*w_0*t) + G_p*exp(1i*(-p)*w_0*t);
    end
end

size = 12;
figure(2)
colororder({'r','g','b','m','k'})
title('Skálázási tétel');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('f(t)','FontSize',size);
xlim([-0.5*T,0.5*T]);
grid on;
legend('Location','northeast','NumColumns',1,'Interpreter','latex','FontSize',size);
hold on
for d = 1 : length(tau)
    plot(t,f(d,:),...
        'DisplayName',['$\Delta$: $\tau$ = T/',num2str(A(d))],...
        'LineStyle','-',...
        'Marker','none',...
        'LineWidth',2);
end
for d = 1 : length(tau)
    plot(t,g(d,:),...
        'DisplayName',['$\sqcap$: $\tau$ = T/',num2str(A(d))],...
        'LineStyle','-.',...
        'Marker','none',...
        'LineWidth',2);
end
hold off


%% Két eltolt négyszögimpulzus, háromszögimpulzus

clc;
clear all;

A = 1;
tau = 1;
N = 20;
w = -N*2*pi/tau : 0.01*pi/tau : N*2*pi/tau;
% Szimmetrikus, tau széles impulzus: F0 "alapjel"
F0 = zeros([length(w),1]);
% Két eltolt alapimpulzus: F 
F = zeros([length(w),1]);
for k = 1 : length(w)
    if w(k) ~= 0
        F0(k) = tau*sin(w(k)*tau/2)/(w(k)*tau/2);
        % F = F0*exp(-jw tau/2) - F0*exp(-jw 3tau/2)
        F(k) = A * F0(k) .* ( exp(-1i*w(k)*tau/2) - exp(-1i*w(k)*3*tau/2) );
    else
        F0(k) = tau;
        F(k) = A*F0(k);
    end
end


% Spektrumból időfv implementálása:
% f(t) = 1/2pi int_{-inf}^{+inf} F(jw) exp(jwt) dw
t = -tau : 0.01*tau : 3*tau;
f0 = 0;
f = 0;
g = 0;
for k = 1 : length(w)-1
    f0 = f0 + (1/(2*pi)) * F0(k)*exp(1i*w(k)*t)*( w(k+1)-w(k) );
    f = f + (1/(2*pi)) * F(k)*exp(1i*w(k)*t)*( w(k+1)-w(k) );
end
% 

size = 12;
figure(1)
title('Két eltolt négyszögimpulzus');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('x(t)','FontSize',size);
xlim([min(t),max(t)]);
xticks(-tau : 0.5*tau : 3*tau);
xticklabels({'-\tau','-\tau/2','0','\tau/2','\tau','3\tau/2','2\tau','5\tau/2','3\tau'})
grid on;
% 'Interpreter','latex' --> "szép" dolgok írásához!!!
legend('Location','northeast','NumColumns',1,'Interpreter','latex','FontSize',size);
hold on
plot(t,real(f0),...
    'DisplayName','$\sqcap(t)$',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(t,real(f),...
    'DisplayName','$\sqcap(t - \tau/2) + \sqcup(t + 3\tau/2)$',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
hold off

% --- Háromszög imp. a négyszögből ---
% G: időbeli konvolúció = spektrális szorzás
G = zeros([length(w),1]);
% H: deriválással
H = zeros([length(w),1]);
for k = 1 : length(w)
    if w(k) ~= 0
        % időbeli konv. önmagával = spektrális négyzetre emelés
        G(k) = A * F0(k).*F0(k) .* exp(-1i*w(k)*tau);
        % időbeli deriválás = spektrumon jw-val szorzás
        % --> DE! most "fordított" az irány, jw-val osztunk
        H(k) = F(k)/(1i*w(k));
    else
        G(k) = A*F0(k)*F0(k);
        H(k) = 0;
    end
end

% Spektrumból időfv implementálása:
g = 0;
h = 0;
for k = 1 : length(w)-1
    g = g + (1/(2*pi)) * G(k)*exp(1i*w(k)*t)*( w(k+1)-w(k) );
    h = h + (1/(2*pi)) * H(k)*exp(1i*w(k)*t)*( w(k+1)-w(k) );
end

size = 12;
figure(2)
title('Háromszögimpulzus');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('x(t)','FontSize',size);
%xlim([-T,T]);
%xticks(-T : 0.5*T : T);
%xticklabels({'-T','-T/2','0','T/2','T'})
%ylim([-1.5,1.5]);
grid on;
legend('Location','northeast','NumColumns',1,'Interpreter','latex','FontSize',size);
hold on
plot(t,real(g),...
    'DisplayName','convolution',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
plot(t,real(h),...
    'DisplayName','$\frac{\partial}{\partial t} \rightarrow j\omega$',...
    'LineStyle','--',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[1.0,0.0,0.0]);
hold off



%% Belépő, lecsengő jel

clc;
clear all;

A = 1;
% három különböző "sajátérték"
lam = [0.25, 1, 4];
t = 0 : 0.01 : 3*(1/min(lam));
f = zeros([length(t),length(lam)]);
for k = 1 : length(lam)
    f(:,k) = A*exp(-lam(k)*t);
end


size = 12;
figure(1)
title('Belépő, lecsengő jel (időben)');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('x(t)','FontSize',size);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(t,f(:,1),...
    'DisplayName',['Ae^{-\alpha t}, \alpha = ',num2str(lam(1))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(t,f(:,2),...
    'DisplayName',['Ae^{-\alpha t}, \alpha = ',num2str(lam(2))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
plot(t,f(:,3),...
    'DisplayName',['Ae^{-\alpha t}, \alpha = ',num2str(lam(3))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
hold off

w = -10 : 0.1 : 10;
F = zeros([length(w),length(lam)]);
for k = 1 : length(lam)
    F(:,k) = 1./(lam(k) + 1i*w);
end

size = 12;
figure(2)
title('Belépő, lecsengő jel (frekvenciában)');
set(gca,'fontsize',size)
xlabel('\omega','FontSize',size);
ylabel('X(j\omega)','FontSize',size);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(w,abs(F(:,1)),...
    'DisplayName',['1/(\alpha + j\omega), \alpha = ',num2str(lam(1))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(w,abs(F(:,2)),...
    'DisplayName',['1/(\alpha + j\omega), \alpha = ',num2str(lam(2))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
plot(w,abs(F(:,3)),...
    'DisplayName',['1/(\alpha + j\omega), \alpha = ',num2str(lam(3))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
hold off

%% Sávszélesség

clc;
clear all;

sigma = 0.1;
A = 8;
lam = 0.2;

w = -15 : 0.1 : 15;
U = A./sqrt(lam^2 + w.^2);

w_max = lam*sqrt(1/sigma^2 - 1);

size = 12;
figure(2)
title('Belépő, lecsengő jel - sávkorlát');
set(gca,'fontsize',size)
xlabel('\omega [MHz]','FontSize',size);
ylabel('U(j\omega) [Vs]','FontSize',size);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(w,abs(U),...
    'DisplayName','|U(j\omega)|',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
yline(sigma*max(abs(U)), ...
    'DisplayName',['\sigma |U(j\omega)|_{max} = ',num2str(sigma*max(abs(U)))], ...
    'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0.5,0.5,0.5])
xline(-w_max, ...
    'DisplayName',['-\omega_{max} = ',num2str(-w_max,'%.3f'), '[MHz]'], ...
    'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0.5,0.5,0.5])
xline(w_max, ...
    'DisplayName',['\omega_{max} = ',num2str(w_max,'%.3f'),' [MHz]'], ...
    'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0.5,0.5,0.5])
hold off

%% Modulált exponenciális jel + korlátok

clc;
clear all;

sigma = 0.1;
A = 8;
lam = 0.2;
tau = 1/lam;
w_0 = 5;
t = 0 : 0.01 : 5*tau;
u = A*exp(-lam*t).*cos(w_0*t);

size = 12;
figure(1)
title('Sávszélesség');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('x(t)','FontSize',size);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(t,A*exp(-lam*t),...
    'DisplayName','A e^{-\alpha t} \epsilon (t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,0.0]);
plot(t,-A*exp(-lam*t),...
    'DisplayName','-A e^{-\alpha t} \epsilon (t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,0.0]);
plot(t,u,...
    'DisplayName',['1/(\alpha + j\omega), \alpha = ',num2str(lam(1))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
hold off



w = -15 : 0.1 : 15;
U_1 = 0.5*A./sqrt(lam^2 + (w+w_0).^2);
U_2 = 0.5*A./sqrt(lam^2 + (w-w_0).^2);
U = U_1 + U_2;

w_max = lam*sqrt(1/sigma^2 - 1);

size = 12;
figure(2)
title('Modulált jel - sávkorlát');
set(gca,'fontsize',size)
xlabel('\omega [MHz]','FontSize',size);
ylabel('U(j\omega) [Vs]','FontSize',size);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(w,abs(U),...
    'DisplayName','|U(j(\omega \pm \omega_{0}))|',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,0.0]);
plot(w,abs(U_1),...
    'DisplayName','|U(j(\omega+\omega_{0}))|',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(w,abs(U_2),...
    'DisplayName','|U(j(\omega-\omega_{0}))|',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
yline(sigma*max(abs(U)), ...
    'DisplayName',['\sigma |U(j\omega)|_{max} = ',num2str(sigma*max(abs(U)))], ...
    'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0.5,0.5,0.5])
xline(-w_max-w_0, ...
    'DisplayName',['-\omega_{0}-\omega_{max} = ',num2str(-w_0-w_max,'%.3f'), '[MHz]'], ...
    'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0.5,0.5,0.5])
xline(w_max-w_0, ...
    'DisplayName',['-\omega_{0}+\omega_{max} = ',num2str(-w_0+w_max,'%.3f'),' [MHz]'], ...
    'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0.5,0.5,0.5])
xline(-w_max+w_0, ...
    'DisplayName',['\omega_{0}-\omega_{max} = ',num2str(w_0-w_max,'%.3f'), '[MHz]'], ...
    'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0.5,0.5,0.5])
xline(w_max+w_0, ...
    'DisplayName',['\omega_{0}+\omega_{max} = ',num2str(w_0+w_max,'%.3f'),' [MHz]'], ...
    'LineStyle','-.',...
    'LineWidth',2,...
    'Color',[0.5,0.5,0.5])
hold off

%% Sok koszinusz

clc;
clear all;

sigma = 0.1;
A = 8;
lam = 0.2;
tau = 1/lam;
w_0 = 5; % 20;

t = 0 : 0.01 : 5*tau;
u = A*exp(-lam*t).*(0 + ...
    4/(1*pi)*cos(w_0*t-pi/2) + 4/(3*pi)*cos(3*w_0*t-pi/2) + 4/(5*pi)*cos(5*w_0*t-pi/2) + ...
    4/(7*pi)*cos(7*w_0*t-pi/2) + 4/(9*pi)*cos(9*w_0*t-pi/2) + 4/(11*pi)*cos(11*w_0*t-pi/2)+ ...
    4/(13*pi)*cos(13*w_0*t-pi/2) + 4/(15*pi)*cos(15*w_0*t-pi/2) + 4/(17*pi-pi/2)*cos(17*w_0*t-pi/2) );


size = 12;
figure(1)
title('Belépő, lecsengő jel szorozva sok koszinusszal');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('u(t)','FontSize',size);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(t,u,...
    'DisplayName','Ae^{-\alpha t} * négyszögjel',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
hold off


w = -100 : 0.1 : 100;
U = 0;
for k = 1 : length(t)-1
    U = U + u(k)*exp(-1i*w*t(k))*(t(k+1)-t(k));
end

size = 12;
figure(2)
set(gca,'fontsize',size)
xlabel('\omega [MHz]','FontSize',size);
ylabel('U(j\omega) [Vs]','FontSize',size);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(w,abs(U),...
    'DisplayName','\Sigma_{p}|U(j(\omega - p\omega_{0}))|',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
hold off

%% Exponenciálisan fel- és lefutó négyszögimpulzus
clc;
clear all;

A = 1;
lam = 1;
tau = 1 / lam;
% három időintervallum
t_0 = -2.5*tau : 0.01 : 0;
t_1 = 0 : 0.01 : 5*tau;
t_2 = 5*tau : 0.01 : 10*tau;
t = [t_0,t_1,t_2];

% lecsengő fv.
f_0 = zeros([1,length(t_0)]);
f_1 = A*(1 - exp(-lam*t_1));
f_2 = f_1(length(t_1))*exp(-lam*t_1);
f = [f_0,f_1,f_2];

% gyorsabb lecsengés
lam_g = 2;
g_0 = f_0;
g_1 = A*(1 - exp(-lam_g*t_1));
g_2 = g_1(length(t_1))*exp(-lam_g*t_1);
g = [g_0,g_1,g_2];

%tau = 1/lam;
w = -10 : 0.1 : 10;
F = 0;
G = 0;
for k = 1 : length(t)-1
    F = F + f(k)*exp(-1i*w*t(k))*(t(k+1)-t(k));
    G = G + g(k)*exp(-1i*w*t(k))*(t(k+1)-t(k));
end


%w0 = 0.01*pi/tau : 0.01*pi/tau : 4*pi/tau;
%F0 = A*tau*sin(w0*tau/2)./(w0*tau/2);
%w0 = [-flip(w0),0,w0];
%F0 = [flip(F0),A*tau,F0];

w = -10 : 0.1 : 10;
F0 = zeros([1,length(w)]);
tau_5 = 5*tau;
for k = 1 : length(w)
    F0(k) = A*tau_5*sin(w(k)*tau_5/2)./(w(k)*tau_5/2)*exp(-1i*w(k)*5*tau_5/2);
end


%t_0 = -2.5*tau : 0.01 : 0;
f_0 = zeros([1,length(t_0)]);

%t_1 = 0 : 0.01 : 5*tau;
f_1 = A*ones([1,length(t_1)]);

%t_2 = 5*tau : 0.01 : 10*tau;
f_2 = zeros([1,length(t_2)]);

%t0 = [t_0,t_1,t_2];
f0 = [f_0,f_1,f_2];


size = 12;
figure(1)
title('Fel- és lefutó négyszögimpulzus (frekiben)');
set(gca,'fontsize',size)
xlabel('\omega','FontSize',size);
ylabel('F(j\omega)','FontSize',size);
grid on;
%legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(w,abs(F0),...
    'DisplayName','asd',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(w,abs(F),...
    'DisplayName','asd',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
plot(w,abs(G),...
    'DisplayName','asd',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
hold off


size = 12;
figure(2)
title('Fel- és lefutó négyszögimpulzus (időben)');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('f(t)','FontSize',size);
xlim([min(t),max(t)]);
grid on;
legend('Location','northeast','NumColumns',1,'Interpreter','latex','FontSize',size);
hold on
plot(t,abs(f0),...
    'DisplayName','$A(\varepsilon(t)-\varepsilon(t-\tau))$',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(t,abs(f),...
    'DisplayName',['$A(e^{-\lambda_{1} t}\varepsilon(t)-' ...
    'e^{-\lambda_{1} t}\varepsilon(t-\tau))$'],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
plot(t,abs(g),...
    'DisplayName',['$A((e^{-\lambda_{2} t}\varepsilon(t)-' ...
    'e^{-\lambda_{2} t}\varepsilon(t-\tau))$'],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
hold off

%% Négyszögjel szorozva szinusz
% időben szorzás: frekiben konvolválás
% spec eset: cos(w_c*t)-vel szorzás --> moduláció

clc;
clear all;

% "alapsávi" f(t) = A*cos(w_0*) vagy f(t) = A*( eps(t+tau/2)-eps(t-tau/2) )
% "vivőjel" w_0 << w_c carrier frequency --> g(t) = B*cos(w_c*t)
% h(t) = f(t)*g(t) = f(t)*B*cos(w_c*t)


% időbeli szorzás: frekiben konvolúció:
% H(jw) = int[ h(t) exp(-j*w*t) ]dt = 
% int[ f(t)*g(t) exp(-j*w*t) ]dt =
% int[ f(t)*cos(w_c*t) exp(-j*w*t) ]dt

% kifejezzük g(t)-t a Fourier-sorával
% g(t) = G_0 + \sum_q[ G_q exp(j*q*w_0*t) ] = 
% G_-1*exp(j*-1*w_c*t) + G_+1*exp(j*+1*w_c*t) = 
% B*exp(j*-1*w_c*t) + B*exp(j*+1*w_c*t)

% H(jw) = int[ f(t)*cos(w_c*t) exp(-j*w*t) ]dt
% \int[ f(t) \sum_q[ G_q exp(j*q*w_0t) ] exp(-j*w*t) ]dt = 
% \sum_q[ G_q * \int[ f(t) * exp(j*q*w_0t) * exp(-j*w*t) ]dt ] = 
% \sum_q[ G_q * F(j( w - p*w_0 ))]


A = 1;
B = 1;
T = 1;
tau = T;
w_0 = 2*pi/T;
w_c = 2*w_0;
w = -100 : 0.1 : 100;
F = zeros([length(w),1]);
H = zeros([length(w),1]);
K = zeros([length(w),1]);
for k = 1 : length(w)
    if w(k) ~= 0
        F(k) = A*tau*sin( w(k)*tau/2 )/( w(k)*tau/2 );
    else
        F(k) = A*tau;
    end
    if w(k) ~= -1*w_c || w(k) ~= 1*w_c
        H(k) = 0.5*B*A*tau*sin( (w(k)-1*w_c)*tau/2 )/( (w(k)-1*w_c)*tau/2 ) + ...
        0.5*B*A*tau*sin( (w(k)+1*w_c)*tau/2 )/( (w(k)+1*w_c)*tau/2 );
    else
        H(k) = 0.5*B*A*tau*1 + 0.5*B*A*tau*1;
    end
    if w(k) ~= -2*w_c || w(k) ~= 2*w_c
        K(k) = 0.5*B*A*tau*sin( (w(k)-2*w_c)*tau/2 )/( (w(k)-2*w_c)*tau/2 ) + ...
        0.5*B*A*tau*sin( (w(k)+2*w_c)*tau/2 )/( (w(k)+2*w_c)*tau/2 );
    else
        K(k) = 0.5*B*A*tau*1 + 0.5*B*A*tau*1;
    end
end



size = 12;
figure(1)
title('Moduláció frekvenciában');
set(gca,'fontsize',size)
xlabel('\omega','FontSize',size);
ylabel('|X(j\omega)|','FontSize',size);
xlim([-3*w_c,3*w_c]);
xticks(-3*w_c : w_c : 3*w_c);
xticklabels({'-3\omega_{c}','-2\omega_{c}','-\omega_{c}','0','\omega_{c}','2\omega_{c}','3\omega_{c}'})
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(w,abs(F),...
    'DisplayName','|F(j\omega)|',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(w,abs(H),...
    'DisplayName','|F(j(\omega \pm \omega_{c}))|',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
plot(w,abs(K),...
    'DisplayName','|F(j(\omega \pm 2\omega_{c}))|',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
hold off


A = 1;
B = 1;
T = 1;
tau = T/2;
w_0 = 2*pi/T;
w_c = 25;
% Időfv. visszaállítás
t = -1*T : 0.005*T : 1*T; % <---- Figyelem, tetszőleges periódust tudunk csinálni!
%f = zeros([length(t),1]);
%h = zeros([length(t),1]);
g = 0.5;
h = 0.0;
d = 0.0;
% H_p = \sum_q[ G_q * F_p-q ]
q = [-1,1];

for p = 1 : 1000
    F_p = A*(tau/T)*sin( (p*w_0)*tau/2 )/( (p*w_0)*tau/2 );
    H_p = 0;
    D_p = 0;
    % Innentől a konvolúció szummája
    for k = 1 : length(q)
        H_p = H_p + 0.5*B*A*(tau/T)*sin( (p*w_0+q(k)*w_c)*tau/2 )/( (p*w_0+q(k)*w_c)*tau/2 );
        D_p = D_p + 0.5*B*A*(tau/T)*sin( (p*w_0+2*q(k)*w_c)*tau/2 )/( (p*w_0+2*q(k)*w_c)*tau/2 );
    end
    g = g + F_p*exp(1i*p*w_0*t) + F_p*exp(1i*(-p)*w_0*t);
    h = h + H_p*exp(1i*p*w_0*t) + H_p*exp(1i*(-p)*w_0*t);
    d = d + D_p*exp(1i*p*w_0*t) + D_p*exp(1i*(-p)*w_0*t);
end


size = 12;
figure(2)
title('Szorzás időben');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('x(t)','FontSize',size);
xlim([-0.5*T,0.5*T]);
xticks(-T : 0.5*T : T);
xticklabels({'-T','-T/2','0','T/2','T'})
ylim([-1.5,1.5]);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(t,d,...
    'DisplayName','f(t)*cos(2\omega_{c}t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
plot(t,h,...
    'DisplayName','f(t)*cos(\omega_{c}t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
plot(t,g,...
    'DisplayName','f(t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
hold off






%% Emelt cos-os jel

clc;
clear all;

A = 1;
T = 1;
tau = T/2;
alfa = [tau/4, tau/8, tau/16];
T_cos = 4*alfa;

t = -T/2 : T/1000 : T/2;
f_square = zeros([1,length(t)]);
f = zeros([length(alfa),length(t)]);
f_cos = 0.5*A*cos(2*pi*t/T) + 0.5*A;

for k = 1 : length(t)
    if -tau/2 < t(k) && t(k) < tau/2
        f_square(k) = A;
    end
end

for p = 1 : length(alfa)
    for k = 1 : length(t)
        if -tau/2 - alfa(p) < t(k) && t(k) < - tau/2 + alfa(p)
            delay = -2*pi*(-tau/2 - alfa(p))/T_cos(p);
            f(p,k) = 0.5*A*cos( 2*pi/(T_cos(p)) * t(k) - delay) + 0.5*A;
        end
        if -tau/2 + alfa(p) <= t(k) && t(k) <= tau/2 - alfa(p)
            f(p,k) = A;
            %g(k) = A;
            %h(k) = A;
        end
        if tau/2 - alfa(p) < t(k) && t(k) < tau/2 + alfa(p)
            delay = 2*pi*(tau/2 - alfa(p))/T_cos(p);
            f(p,k) = 0.5*A*cos( 2*pi/(T_cos(p)) * t(k) - delay) + 0.5*A;
        end
    end
end

size = 12;
figure(1)
title('Emelt cos-os impulzus (időfv)');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('f(t)','FontSize',size);
xlim([-T/2,T/2]);
xticks(-T/2 : T/4 : T/2);
xticklabels({'-T/2','-T/4','0','T/4','T/2'});
%ylim([-1.5,1.5]);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(t,f_square,...
    'DisplayName','négyzögjel',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,0.0]);
plot(t,f(3,:),...
    'DisplayName',['\alpha = \tau/',num2str(0.5/alfa(3))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(t,f(2,:),...
    'DisplayName',['\alpha = \tau/',num2str(0.5/alfa(2))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
plot(t,f(1,:),...
    'DisplayName',['\alpha = \tau/',num2str(0.5/alfa(1))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
plot(t,f_cos,...
    'DisplayName','emelt cos',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,1.0]);
hold off

w_0 = 2*pi/tau;
w = -10*w_0 : 0.01*w_0 : 10*w_0;
F_square = zeros([1,length(w)]);
F = zeros([length(alfa),length(w)]);
F_cos = zeros([1,length(w)]);

for k = 1 : length(t)-1
    F_square = F_square + f_square(k) * exp(-1i*w*t(k)) * (t(k+1)-t(k));
    F_cos = F_cos + f_cos(k) * exp(-1i*w*t(k)) * (t(k+1)-t(k));
end

for p = 1 : length(alfa)
    for k = 1 : length(t)-1
        F(p,:) = F(p,:) + f(p,k) * exp(-1i*w*t(k)) * (t(k+1)-t(k));
    end
end
%}


size = 12;
figure(2)
title('Emelt cos-os impulzus (spektrum)');
set(gca,'fontsize',size)
xlabel('\omega','FontSize',size);
ylabel('F(j\omega)','FontSize',size);
xlim([-7*w_0,7*w_0]);
xticks(-7*w_0 : 1*w_0 : 7*w_0);
%xticklabels({'-T','-T/2','0','T/2','T'})
%ylim([-1.5,1.5]);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
plot(w,real(F_square),...
    'DisplayName','négyszögjel',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,0.0]);
plot(w,real(F(3,:)),...
    'DisplayName',['\alpha = \tau/',num2str(0.5/alfa(3))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
plot(w,real(F(2,:)),...
    'DisplayName',['\alpha = \tau/',num2str(0.5/alfa(2))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[0.2,0.8,0.2]);
plot(w,real(F(1,:)),...
    'DisplayName',['\alpha = \tau/',num2str(0.5/alfa(1))],...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,0.0]);
plot(w,real(F_cos),...
    'DisplayName','emelt cos',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',2,...
    'Color',[1.0,0.0,1.0]);
hold off


