% JR-2, 2. gyak
% Fourier sorfejtés


%% Szimm négyszögjel

% clear command window
clc;
% clear workspace undsoweiter
clear all;

% ---------------------------------------------------------
% négyszögjel amplitúdója
A = 1;
% periódusidő
T = 1;
% periódus felében f(t) = A, különben f(t) = 0
tau = T/2;
% alapkörfrekvencia
w_0 = 2*pi/T;
% harmonikus frekik --> w = w_0*(1 : 20)-ra még jól látszódik a "hullámzás"
w = w_0*(1 : 1000);

% ---------------------------------------------------------
% Komplex Fourier-együtthatók
% pozitív p indexek
Fp = zeros([length(w),1]);
% negatív p indexek
Fn = zeros([length(w),1]);
% p=0 (DC komponens=átlag)
F0 = A*1*(tau/T);

% F = A/T * sin(p*w_0*tau/2)/(p*w_0*tau/2) * tau
for p = 1 : length(w)
    Fp(p) = A*(tau/T)*sin(w(p)*tau/2)/(w(p)*tau/2);
    Fn(p) = A*(tau/T)*sin(-w(p)*tau/2)/(-w(p)*tau/2);
end


% ---------------------------------------------------------
% Frekvenciatománybeli leírás

size = 12;
figure(1)
colororder({'k','k'});
title('Szimmetrikus négyszögjel');
set(gca,'fontsize',size)
xlabel('p*\omega_{0}','FontSize',size);
ylabel('|G_{p}^{c}|,|F_{p}^{c}|','FontSize',size);
xlim([0,w(20)]);
%ylim([-1.5,2.5]);
grid on;
legend('Location','northwest','NumColumns',1,'FontSize',size);
hold on
plot([0,w],[F0;abs(Fp)],...
    'DisplayName','|F_{p}^{c}|',...
    'LineStyle','none',...
    'Marker','+',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
hold off

% ---------------------------------------------------------
% Időtartománybeli jel

% ábrázolt időtartomány: t \in [-T;T]
t = -1.0*T : 0.01*T : 1.0*T;

% Fourier-sor: f(t) = F_0 + Summa_{P=-N}^{+N} F_p*exp(j*p*w_0*t)
f = F0;
for p = 1 : length(w)
    f = f + Fp(p)*exp(1i*w(p)*t) + Fn(p)*exp(1i*(-w(p))*t);
end

%{
% csak pozitív Fourier együtthatók --> teljesítmény fele "elveszik"
f = F0/2;
for p = 1 : length(w)
    f = f + Fp(p)*exp(1i*w(p)*t);
end
%}

size = 12;
figure(2)
colororder({'k','k'});
title('Szimmetrikus négyszögjel');
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
hold off

%% eltolt négyszögjel

% négyszögjel amplitúdója
A = 1;
% periódusidő
T = 1;
% periódus felében g(t) = A, különben g(t) = 0
tau = T/2;
% alapkörfrekvencia
w_0 = 2*pi/T;
% felharmonikus frekik
w = w_0*(1 : 1000);

% ---------------------------------------------------------
% Komplex Fourier-együtthatók
% pozitív indexek
Gp = zeros([length(w),1]);
% negatív indexek
Gn = zeros([length(w),1]);
% p=0 (DC komponens=átlag)
G0 = A*1*1*(tau/T);

% F = A/T * [ exp(-j*p*w_0*tau/2) ] * sin(p*w_0*tau/2)/(p*w_0*tau/2) * tau
for p = 1 : length(w)
    Gp(p) = A*(tau/T)*exp(-1i*w(p)*tau/2)*sin(w(p)*tau/2)/(w(p)*tau/2);
    Gn(p) = A*(tau/T)*exp(-1i*(-w(p))*tau/2)*sin(-w(p)*tau/2)/(-w(p)*tau/2);
end

% ---------------------------------------------------------
% Frekvenciatománybeli leírás

size = 12;
figure(1)
colororder({'k','k'});
title('Eltolt négyszögjel');
set(gca,'fontsize',size)
xlabel('p*\omega_{0}','FontSize',size);
ylabel('|G_{p}^{c}|,|F_{p}^{c}|','FontSize',size);
xlim([0,w(20)]);
%ylim([-1.5,2.5]);
grid on;
legend('Location','northwest','NumColumns',1,'FontSize',size);
hold on
text(50,0.25,'Időbeli eltolás = Frekiben fázisforgatás');
text(50,0.225,'Nem látszódik meg F_{p}^{c}-k absz. értékeiben!');
plot([0,w],[G0;abs(Gp)],...
    'DisplayName','|G_{p}^{c}| = |F_{p}^{c}|*|e^{-j\omega \tau/2}|',...
    'LineStyle','none',...
    'Marker','o',...
    'LineWidth',5,...
    'Color',[1.0,0.0,0.0]);
plot([0,w],[F0;abs(Fp)],...
    'DisplayName','|F_{p}^{c}|',...
    'LineStyle','none',...
    'Marker','+',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
hold off


% ---------------------------------------------------------
% Időtartománybeli jel 

% ábrázolt időtartomány: t \in [-T;T]
t = -1.0*T : 0.01*T : 1.0*T;

% Fourier-sor: g(t) = G_0 + Summa_{P=-N}^{+N} G_p*exp(j*p*w_0*t)
g = G0;
for p = 1 : length(w)
    g = g + Gp(p)*exp(1i*w(p)*t) + Gn(p)*exp(1i*(-w(p))*t);
end

size = 12;
figure(2)
colororder({'k','k'});
title('Eltolt négyszögjel');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('y(t)','FontSize',size);
%xlim([-T,T]);
%ylim([-1.5,2.5]);
grid on;
legend('Location','northwest','NumColumns',1,'FontSize',size);
hold on
plot(t,g,...
    'DisplayName','g(t) = f(t - \tau/2)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[1.0,0.0,0.0]);
plot(t,f,...
    'DisplayName','f(t)',...
    'LineStyle','-',...
    'Marker','none',...
    'LineWidth',3,...
    'Color',[0.0,0.0,1.0]);
hold off

%% Kitöltési tényező

% négyszögjel amplitúdója
A = 1;
% periódusidő
T = 1;
% periódus felében/harmadában/negyedében/... f(t) = A, különben f(t) = 0
tau = [ T/2, T/3, T/4, T/5, T/10 ];
% alapkörfrekvencia
w_0 = 2*pi/T;
% felharmonikus frekik
w = w_0*(1 : 1000);

Fp = zeros([length(w),length(tau)]);
Fn = zeros([length(w),length(tau)]);
F0 = A*1*1*(tau/T);

for k = 1 : length(tau)
    for p = 1 : length(w)
        Fp(p,k) = A*(tau(k)/T)*exp(-1i*w(p)*tau(k)/2)*sin(w(p)*tau(k)/2)/(w(p)*tau(k)/2);
        Fn(p,k) = A*(tau(k)/T)*exp(-1i*(-w(p))*tau(k)/2)*sin(-w(p)*tau(k)/2)/(-w(p)*tau(k)/2);
    end
end

% ---------------------------------------------------------
% Frekvenciatománybeli leírás

size = 12;
figure(1)
colororder({'k','r','g','b','m'});
title('Kitöltési tényező');
set(gca,'fontsize',size)
xlabel('p*\omega_{0}','FontSize',size);
ylabel('|F_{p}^{c}|','FontSize',size);
xlim([0,w(20)]);
%ylim([-1.5,2.5]);
grid on;
legend('Location','northeast','NumColumns',2,'FontSize',size);
hold on
text(30,0.300,'\tau = 0.50*T: T/\tau = 2, minden második p*\omega_{0}-ra nulla');
text(30,0.275,'\tau = 0.33*T: T/\tau = 3, minden harmadik p*\omega_{0}-ra nulla');
text(30,0.250,'\tau = 0.25*T: T/\tau = 4, minden negyedik p*\omega_{0}-ra nulla');
for k = 1 : length(tau)
    plot([0,w],[F0(k);abs(Fp(:,k))],...
        'DisplayName',['|F_{p}^{c}|, \tau = ',num2str(tau(k)/T,'%.2f'),'*T'],...
        'LineStyle','-',...
        'Marker','none',...
        'LineWidth',3);
end
hold off

% ---------------------------------------------------------
% Időtartománybeli jel 

% ábrázolt időtartomány: t \in [-T;T]
t = -1.0*T : 0.01*T : 1.0*T;

% Fourier-sor: g(t) = G_0 + Summa_{P=-N}^{+N} G_p*exp(j*p*w_0*t)
f = zeros([length(t),length(tau)]);
for k = 1 : length(tau)
    f(:,k) = F0(k);
    for p = 1 : length(w)
        f(:,k) = f(:,k) + Fp(p,k)*exp(1i*w(p)*t') + Fn(p,k)*exp(1i*(-w(p))*t');
    end
end

size = 12;
figure(2)
colororder({'k','r','g','b','m'});
title('Kitöltési tényező');
set(gca,'fontsize',size)
xlabel('t','FontSize',size);
ylabel('y(t)','FontSize',size);
%xlim([-T,T]);
%ylim([-1.5,2.5]);
grid on;
legend('Location','northeast','NumColumns',1,'FontSize',size);
hold on
for k = 1 : length(tau)
    plot(t,f(:,k),...
        'DisplayName',['f(t), \tau = ',num2str(tau(k)/T,'%.2f'),'*T'],...
        'LineStyle','-',...
        'Marker','none',...
        'LineWidth',3);
end
hold off

%% Háromszögjel négyszögjelből

% ---------------------------------------------------------
% Négyszögjel: g(t) = B, t \in [0;T/2], különben g(t) = -B

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
title('Eltolt négyszögjel');
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
title('Eltolt négyszögjel');
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

% ---------------------------------------------------------------------------------------------------
% 3. gyakorlat anyaga

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
