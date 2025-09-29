% JR-2, 3. gyakorlat (06.19.)
% Fourier-sorfejtés vége, Fourier-transzformáció

%% Háromszögjel négyszögjelből - deriválással

% ---------------------------------------------------------
% Periodikus négyszögjel: g(t) = B, t \in [0;T/2], különben g(t) = -B

clc;
clear all;

A = 1;
T = 1;
B = 4*A/T;
w_0 = 2*pi/T;
p = 1 : 1000;
w = w_0*p;

Gp = zeros([length(p),1]);
Gn = zeros([length(p),1]);
G0 = 0;

for k = 1 : length(p)
    Gp(k) = 1i*B/(k*pi)*( cos(k*pi) - 1 );
    Gn(k) = 1i*B/((-k)*pi)*( cos((-k)*pi) - 1 );
end

Fp = zeros([length(p),1]);
Fn = zeros([length(p),1]);
F0 = 0;

for k = 1 : length(p)
    Fp(k) = Gp(k)/(1i*w(k));
    Fn(k) = Gn(k)/(-1i*w(k));
end

size = 12;
figure(1)
colororder({'k','k'});
title('Háromszögjel deriválással');
set(gca,'fontsize',size)
xlabel('p*\omega_{0}','FontSize',size);
ylabel('|G_{p}^{c}|,|F_{p}^{c}|','FontSize',size);
xlim([0,w(20)]);
%ylim([-1.5,2.5]);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
text(25,1.5,'Gyorsabb (négyzetes) lecsengés, páros p\omega_{0} továbbra is 0!');
plot([0,w],[G0;abs(Gp)],...
    'DisplayName','|G_{p}^{c}|',...
    'LineStyle','none',...
    'Marker','o',...
    'LineWidth',5,...
    'Color',[1.0,0.0,0.0]);
plot([0,w],[F0;abs(Fp)],...
    'DisplayName','|F_{p}^{c}| = |G_{p}^{c}|/(j*p*\omega_{0})',...
    'LineStyle','none',...
    'Marker','+',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
hold off

t = -1.0*T : 0.01*T : 1.0*T;

f = F0;
g = G0;

for k = 1 : length(w)
    f = f + Fp(k)*exp(1i*w(k)*t) + Fn(k)*exp(1i*(-w(k))*t);
    g = g + Gp(k)*exp(1i*w(k)*t) + Gn(k)*exp(1i*(-w(k))*t);
end

size = 12;
figure(2)
colororder({'k','k'});
title('Háromszögjel deriválással');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('y(t)','FontSize',size);
%xlim([-T,T]);
%ylim([-1.5,2.5]);
grid on;
legend('Location','northwest','NumColumns',1,'FontSize',size);
hold on
plot(t,f,...
    'DisplayName','f(t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
plot(t,g,...
    'DisplayName','g(t) = d( f(t) )/dt',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[1.0,0.0,0.0]);

hold off


%% Két négyszögjel időtartománybeli konvolúciója

A = 1;
T = 1;
tau = T/2;
w_0 = 2*pi/T;
w = w_0*(1 : 1000);

Fp = zeros([length(w),1]);
Fn = zeros([length(w),1]);
F0 = A*(tau/T);

% Négyszögjel
for p = 1 : length(w)
    Fp(p) = A*(tau/T)*exp(-1i*w(p)*tau/2)*sin(w(p)*tau/2)/(w(p)*tau/2);
    Fn(p) = A*(tau/T)*exp(-1i*(-w(p))*tau/2)*sin(-w(p)*tau/2)/(-w(p)*tau/2);
end

Gp = zeros([length(w),1]);
Gn = zeros([length(w),1]);
G0 = F0*F0;

% Háromszögjel: négyszögjelek szorzata!!!
for p = 1 : length(w)
    %Gp(p) = Fp(p)*Fp(p);
    %Gn(p) = Fn(p)*Fn(p);
    Gp(p) = A^2*(tau/T)^2*exp(-1i*w(p)*tau/2*2)*(sin(w(p)*tau/2)/(w(p)*tau/2))^2;
    Gn(p) = A^2*(tau/T)^2*exp(-1i*(-w(p))*tau/2*2)*(sin(-w(p)*tau/2)/(-w(p)*tau/2))^2;
end


% ---------------------------------------------------------
t = -1.0*T : 0.01*T : 1.0*T;

f = F0;
g = G0;
for p = 1 : length(w)
    f = f + Fp(p)*exp(1i*w(p)*t) + Fn(p)*exp(1i*(-w(p))*t);
    g = g + Gp(p)*exp(1i*w(p)*t) + Gn(p)*exp(1i*(-w(p))*t);
end

size = 12;
figure(2)
colororder({'k','k'});
title('Konvolúció');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('y(t)','FontSize',size);
%xlim([-T,T]);
%ylim([-1.5,2.5]);
grid on;
legend('Location','northwest','NumColumns',1,'FontSize',size);
hold on
plot(t,f,...
    'DisplayName','f(t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
plot(t,g,...
    'DisplayName','g(t) = conv(f(t),f(t))',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[1.0,0.0,0.0]);
hold off


%% Vadabb konvolválások

Fp = zeros([length(w),1]);
Fn = zeros([length(w),1]);
F0 = 0; %F0 = A*(tau/T);
Gp = zeros([length(w),1]);
Gn = zeros([length(w),1]);
G0 = 0; %G0 = F0*F0;
Hp = zeros([length(w),1]);
Hn = zeros([length(w),1]);
H0 = 0; %H0 = G0*G0;
Ip = zeros([length(w),1]);
In = zeros([length(w),1]);
I0 = 0; %I0 = H0*H0;
Jp = zeros([length(w),1]);
Jn = zeros([length(w),1]);
J0 = 0; %J0 = I0*I0;


for p = 1 : length(w)
    Fp(p) = A*(tau/T)*exp(-1i*w(p)*tau/2)*sin(w(p)*tau/2)/(w(p)*tau/2);
    Fn(p) = A*(tau/T)*exp(-1i*(-w(p))*tau/2)*sin(-w(p)*tau/2)/(-w(p)*tau/2);
    % ---------------------------------
    Gp(p) = Fp(p)*Fp(p);
    Gn(p) = Fn(p)*Fn(p);
    % ---------------------------------
    %Hp(p) = (Fp(p)*Fp(p))^2;
    %Hn(p) = (Fn(p)*Fn(p))^2;
    Hp(p) = Gp(p)*Gp(p);
    Hn(p) = Gn(p)*Gn(p);
    % ---------------------------------
    Ip(p) = Hp(p)*Hp(p);
    In(p) = Hn(p)*Hn(p);
    % ---------------------------------
    Jp(p) = Ip(p)*Ip(p);
    Jn(p) = In(p)*In(p);
end

f = F0;
h = H0;
g = G0;
i = I0;
j = J0;
for p = 1 : length(w)
    f = f + Fp(p)*exp(1i*w(p)*t) + Fn(p)*exp(1i*(-w(p))*t);
    g = g + Gp(p)*exp(1i*w(p)*t) + Gn(p)*exp(1i*(-w(p))*t);
    h = h + Hp(p)*exp(1i*w(p)*t) + Hn(p)*exp(1i*(-w(p))*t);
    i = i + Ip(p)*exp(1i*w(p)*t) + In(p)*exp(1i*(-w(p))*t);
    j = j + Jp(p)*exp(1i*w(p)*t) + Jn(p)*exp(1i*(-w(p))*t);
end

% folyton csökken az amplitúdó --> érdemes a max-hoz normalizálni
%
f = f./max(f);
g = g./max(g);
h = h./max(h);
i = i./max(i);
j = j./max(j);
%}


size = 12;
figure(2)
colororder({'k','k'});
title('Konvolúció');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('y(t)','FontSize',size);
ylim([-1.25,1.75]);
grid on;
legend('Location','northwest','NumColumns',1,'FontSize',size);
hold on
plot(t,f,...
    'DisplayName','f(t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
plot(t,g,...
    'DisplayName','g(t) = conv(f(t),f(t))',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[1.0,0.0,0.0]);
plot(t,h,...
    'DisplayName','h(t) = conv(g(t),g(t))',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[0.2,0.8,0.2]);
% i és j már nem "érdekes", ami "érdekes": a MATLAB engedi a
% felüldefiniálást, mert amúgy az i és j a komplex egység jele!!!
%{
plot(t,i,...
    'DisplayName','i(t) = conv(h(t),h(t))',...
    'LineStyle',':',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[1.0,0.0,1.0]);
plot(t,j,...
    'DisplayName','j(t) = conv(i(t),i(t))',...
    'LineStyle','none',...
    'Marker','+',...
    'LineWidth',2,...
    'Color',[0.0,0.0,0.0]);
%}
hold off

%% Innentől Fourier-transzformáció -------------------------------------

%% Szimmetrikus négszögimpulzus

clc;
clear all;

% Szimmetrikusan, t \in [-tau/2; +tau/2] között f(t) = A, különben 0
A = 1;
tau = 1;

% pozitív féloldal:
w = 0.1*pi/tau : 0.1*pi/tau : 10*pi/tau;
F = A*tau*sin(w*tau/2)./(w*tau/2);

% teljes, kétoldalas jel + nullabeli érték
%w = -10*pi/tau : 0.1*pi/tau : 10*pi/tau;
w = [-flip(w),0,w];
F = [flip(F),A*tau,F];


size = 12;
figure(1)
title('Szimmetrikus négyszögimpulzus');
set(gca,'fontsize',size)
xlabel('\omega','FontSize',size);
ylabel('F(j\omega)','FontSize',size);
xlim([-10*pi/tau,10*pi/tau]);
xticks(-10*pi/tau : 2*pi/tau : 10*pi/tau);
xticklabels({'-10\pi/\tau','-8\pi/\tau','-6\pi/\tau','-4\pi/\tau','-2\pi/\tau', ...
             '0','2\pi/\tau','4\pi/\tau','6\pi/\tau','8\pi/\tau','10\pi/\tau',})
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
rectangle('Position',[-9.5*pi/tau, 0.55 5*pi/tau 0.4],'FaceColor',[0.9,0.9,0.9])
text(-9*pi/tau,0.9,'Nullhelyek:');
text(-9*pi/tau,0.8,'sin(\omega \tau/2) = 0');
text(-9*pi/tau,0.7,'\omega \tau/2 = k \pi');
text(-9*pi/tau,0.6,'\omega = k 2\pi/\tau');
plot(w,F,...
    'DisplayName','F(j\omega)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
plot(w,abs(F),...
    'DisplayName','|F(j\omega)|',...
    'LineStyle','-.',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[1.0,0.0,0.0]);
hold off

%% Teljesítmény, Perseval-tétele
A = 2;
% -----------------------------------------------------
% frekvencia tartomány

% Pozitív oldalon numerikusan integrálunk
w = 0.01*pi/tau : 0.01*pi/tau : 25*pi/tau;
F = A*tau*sin(w*tau/2)./(w*tau/2);
E_w = 0;
for k = 1 : length(w)-1
    E_w = E_w + ( abs(F(k))^2 )*(w(k+1)-w(k));
end
% Negatív w-k oldala (szimmetria okán csak x2)
E_w = 2*E_w;
% w=0 esete
E_w = E_w + A*tau*(2*w(1));
% 1/2pi tényező:
E_w = E_w/(2*pi);

% -----------------------------------------------------
% idő tartomány

E_t = abs(A)^2 *tau;


% -----------------------------------------------------
% Összehasonlítás:
clc;
disp('Energia:');
disp(['E_w = ',num2str(E_w,'%.3f')]);
disp(['E_t = ',num2str(E_t,'%.3f')]);
disp('(koherens egységrendszerben)');

%% Szögfüggvény

clc;
clear all;

A = 1;
T = 1;
w_0 = 2*pi/T;
p = [-1,1];
w = w_0*p;
F = [0.5*A*(2*pi),0.5*A*(2*pi)];

size = 12;
figure(1)
title('Koszinusz spektruma');
set(gca,'fontsize',size)
xlabel('\omega','FontSize',size);
ylabel('|X(j\omega)|','FontSize',size);
xlim([-2*w_0,2*w_0]);
xticks(-2*w_0 : w_0 : 2*w_0);
xticklabels({'-2\omega_{0}','-1\omega_{0}','0','1\omega_{0}','2\omega_{0}'})
ylim([0,1.25*pi]);
yticks(0 : 0.5*pi : 1.5*pi);
yticklabels({'0','0.5\pi','\pi','1.5\pi'});
grid on;
%legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
yline(0.5*A*(2*pi),'LineStyle','--','Color',[0.5,0.5,0.5]);
stem(w,abs(F),...
    'DisplayName','|F(j\omega)|',...
    'LineStyle','-',...
    'Marker','*',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
hold off


%% Sok szögfüggvény

clc;
clear all;

A = 1;
T = 1;
w_0 = 2*pi/T;
N = 7;
F = zeros([1,N]);
F_0 = A/2;
for p = 1 : N
    if mod(p,2) == 1
        F(p) = 4*A/(p*pi);
    end
end
w = (-N : 1 : N)*w_0;
F = [flip(F),F_0,F];

size = 12;
figure(1)
title('Szimmetrikus négyszögJEL spektruma');
set(gca,'fontsize',size)
xlabel('\omega','FontSize',size);
ylabel('|X(j\omega)|','FontSize',size);
xlim([-7*w_0,7*w_0]);
xticks([-7*w_0 : 2*w_0 : -1*w_0,0,1*w_0 : 2*w_0 : 7*w_0]);
xticklabels({'-7\omega_{0}','-5\omega_{0}','-3\omega_{0}','-1\omega_{0}','0' ...
    '1\omega_{0}','3\omega_{0}','5\omega_{0}','7\omega_{0}',})
ylim([0,0.6*pi]);
yticks(0 : 0.125*pi : 0.5*pi);
yticklabels({'0','\pi/8','\pi/4','3\pi/8','\pi/2'});
grid on;
%legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
yline(4*A/(1*pi),'LineStyle','--','Color',[0.5,0.5,0.5]);
yline(4*A/(3*pi),'LineStyle','--','Color',[0.5,0.5,0.5]);
yline(4*A/(5*pi),'LineStyle','--','Color',[0.5,0.5,0.5]);
yline(4*A/(7*pi),'LineStyle','--','Color',[0.5,0.5,0.5]);
text(-6*w_0,pi/2+0.1,'Periodikus jel spektruma vonalas!')
text(-6*w_0,pi/2,['f(t) = 1 + 4/\pi cos(\omega_{0}t) + ' ...
    '4/3\pi cos(3\omega_{0}t) + 4/5\pi cos(5\omega_{0}t) + ...'])
stem(w,F,...
    'DisplayName','F(j\omega)',...
    'LineStyle','-',...
    'Marker','*',...
    'LineWidth',2,...
    'Color',[0.0,0.0,1.0]);
hold off











