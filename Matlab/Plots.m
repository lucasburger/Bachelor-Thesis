% Erstellt die Plots der Bachelorarbeit
clear all; close all; clc

%% Plotten der Auszahlungsfunktionen der Optionen
clear all; close all; clc
K = 30;
S = 0:30:60;
Call = max(S-K,0);
Put = max(K-S,0);
payoff = figure(1);

%%%%%%%%%%%% CALL-OPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(S,Call,'r','LineWidth',5);

%%%%%%%%%%%%%%%%%%%% PUT-OPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on 
plot(S,Put,'b','LineWidth',5);
axis([0 60 0 40]);
axis off;
titlestring = sprintf('Auszahlungsfunktion der Optionen \n Strike Price = K');
title(titlestring)
xlabel('Aktie');
ylabel('Auszahlung');
text(28.5,-2,'\fontsize{20} K');
text(-3.5,29.7,'\fontsize{20} K');
text(61,-1.5,'\fontsize{20} S');
text(-7,42.5,'\fontsize{12} Auszahlung');
%annotation('textarrow',[.8 .5],[.3 .5],'String','Textpfeil','HeadStyle','plain');

% Position des aktuellen Axes-Objekts ermitteln
p = get(gca, 'Position');
 
% Pfeile an die X- und Y-Achsen setzen
a = annotation('arrow');
set(a, 'X', [p(1) p(1)], 'Y', [p(2) p(2)+p(4)+0.025],'LineWidth',2);
a = annotation('arrow');
set(a, 'X', [p(1) p(1)+p(3)+0.025], 'Y', [p(2) p(2)],'LineWidth',2);

legend({'Call','Put'},'Location','north','FontSize',15);

if (1==0)
%%%%%%% SPEICHERN DES PLOTS %%%%%%%%%%%%%%%%
print(payoff,'-djpeg','-r250',['/Users/lucasburger/Dropbox/Studium/Bachelorarbeit'...
     '/Bilder/PayoffOptionen.jpg'])
end

%% Plot der geometrischen Brownschen Bewegung
clear all; clc; close all;

T = 1; N = 14; dt = T/N; 
t = 0:dt:T; 
dW = sqrt(dt)*[1.1 0.3 0.6 -0.4 0.8 -0.3 0.4 -0.7 0.05 0.2 -0.4 0.45 -0.3 0.5]; 
mu = 0.2; sigma = 0.3; 
BBW = figure(2);
y = exp(mu*t + sigma*[0,cumsum(dW)]);

plot(t,y,'LineWidth',7,'Color',[0.64 0.71 0.8]);
axis([min(t) max(t)+0.2 min(y) max(y)+0.1]);
%pb = get(gca,'Position');
pb = [0.130000000000000; 0.110000000000000; 0.775000000000000; 0.815000000000000];
axis off;
a = annotation('arrow');
set(a, 'X', [pb(1) pb(1)], 'Y', [pb(2) pb(2)+pb(4)+0.025],'LineWidth',6,'Color',[0.64 0.71 0.8]);
a = annotation('arrow');
set(a, 'X', [pb(1) pb(1)+pb(3)-0.01], 'Y', [pb(2) pb(2)],'LineWidth',6,'Color',[0.64 0.71 0.8]);
%a = annotation('arrow');
%set(a,'X', [725/1000 786/1000],'Y', [600/1000 770/1000],'LineWidth',7,'Color',[0.64 0.71 0.8]);

if (1==0)
    print(BBW,'-djpeg','-r350',['/Users/lucasburger/Dropbox/Studium/Bachelorarbeit/Bilder/BBW.jpg']);
end


%% Plots der amerikanischen und Europäischen Put Option mit diskretisierung
clear all; clc; close all;

K = 100;
r = 0.1;
sigma = 0.35;
T=1;
D0 = 0.05;

% Optionswerte ausrechnen
resultA = AmericanPut(K,r,D0,sigma,T);
resultE = EuropeanPut(K,r,D0,sigma,T);


Putoption = figure(3);
% Payoff plotten
plot(resultA(1,:),max(0,K-resultA(1,:)),'--k');
hold on
% Optionswerte Plotten
plot(resultA(1,:),resultA(2,:),'k');
plot(resultE(1,:),resultE(2,:),'-.k');

i = 0:0.5:120;
[~,X] = blsprice(i,K,r,T,sigma,D0);
hold on
%plot(i,X,'-.b');


% Plot verschönern
xlim([0 2*K]);
legend('Payoff','amerikanisch','europäisch','Location','northeast');
xlabel('Aktienkurs S');
ylabel('Optionspreis P');

if (1==0)
    print(Putoption,'-djpeg','-r350',['/Users/lucasburger/Dropbox/Studium/Bachelorarbeit/Bilder/PlotPutOption.jpg']);
end


%% Plots der amerikanischen und Europäischen Call Option mit diskretisierung
clear all; clc; close all;

K = 100;
r = 0.1;
sigma = 0.35;
T=1;
D0 = 0.08;


% Optionswerte ausrechnen
resultA = AmericanCall(K,r,D0,sigma,T);
resultE = EuropeanCall(K,r,D0,sigma,T);


Calloption = figure(4);
% Payoff plotten
plot(resultA(1,:),max(0,resultA(1,:)-K),'--k');
hold on
% Optionswerte Plotten
plot(resultA(1,:),resultA(2,:),'k');
plot(resultE(1,:),resultE(2,:),'-.k');

% Plot verschönern
xlim([80 2.5*K]);
legend('Payoff','amerikanisch','europäisch','Location','northwest');
xlabel('Aktienkurs S');
ylabel('Optionspreis P');

if (1==0)
    print(Calloption,'-djpeg','-r350',['/Users/lucasburger/Dropbox/Studium/Bachelorarbeit/Bilder/PlotCallOption.jpg']);
end

%% Erstellen der Konvergenzplots der Call-Option
clear all; clc; close all;
S0 = 20;
K = 18;
r = 0.1;
sigma = 0.35;
T=1;

KonvergenzCall = figure(5);
Von = 100;
Schritt = 3;
Bis = 600;

N = Von:Schritt:Bis;
s = size(N);
s = s(2);

[CallBS,PutBS] = blsprice(S0,K,r,T,sigma);
CallPut = nan(2,s);
for i = 1:s
    [Call,Put] = BinbaumEuro(S0,K,r,sigma,T,N(i));
    CallPut(1,i) = Call;
    %CallPut(2,i) = Put;
end
plot([Von Bis],[CallBS CallBS],'g','Linewidth',2);
hold on
for i = 1:s
    if (mod(N(i),2) == 0)
        plot(N(i),CallPut(1,i),'*b');
    elseif (mod(N(i),2) == 1)
        plot(N(i),CallPut(1,i),'*r');
    end
end

xlabel('Schrittanzahl');
ylabel('Optionspreis');
legend('Black-Scholes','Gerade N','Ungerade N');
titleString = sprintf('K = %d',K) ;
title(titleString) ;


% Fehlerplot
fehler = figure(5);
plot(N,CallPut(1,:)-CallBS,'*','Color',[0.09 0.45 0.8]);
hold on
plot(N,1./N,'Color',[0.9 0 0])
plot(N,CallPut(1,:)-CallBS,'-k');

plot(N,-1./N,'Color',[0.9 0 0])

xlabel('Schrittanzahl');
ylabel('Fehler');
axis([100 600 -0.006 0.006]);
legend('Fehler','Schranke')
title(titleString);

if (1==0)
    zielfehler = sprintf('/Users/lucasburger/Dropbox/Studium/Bachelorarbeit/Bilder/FehlerCall%d.jpg',K);
    zielkonvergenz = sprintf('/Users/lucasburger/Dropbox/Studium/Bachelorarbeit/Bilder/KonvergenzCall%d.jpg',K);
    print(fehler,'-djpeg','-r350',zielfehler);
    print(KonvergenzCall,'-djpeg','-r350',zielkonvergenz);
    close all;
end

%% Verfüngung Binomialbaum
clear all; close all; clc
S0 = 20;
K = 18;
r = 0.1;
sigma = 0.35;
T=1;

Konvergenz = figure(6);
Strike18 = [18 18];
Strike20 = [20 20];
hold on
plot([0 50],Strike18,'-.g','Linewidth',1.5);
plot([0 50],Strike20,'-.b','Linewidth',1.5);

for N = 1:50
    dt = T/N;
    u = exp(sigma*sqrt(dt));
    d = 1/u;
    S = S0*(u.^(N:-1:0)'.*d.^(0:N)');
    n = N*ones(1,length(S));
    hold on
    plot(n,S,'Color',[0.74,0.82,0.93]);
    plot(n,S,'*k');
    axis([0 50 10 30])
end
xlabel('Schrittanzahl');
ylabel('Aktienkurse in T');
legend('Strike = 18','Strike = 20');

if (1==0)
    print(Konvergenz,'-djpeg','-r350',['/Users/lucasburger/Dropbox/Studium/Bachelorarbeit/Bilder/Verjuengung.jpg']);
end

%% Plot der amerikanischen Put-Option im zeitverlauf

clear all; clc; close all;

K = 100;
r = 0.1;
sigma = 0.35;
T=1;
D0 = 0.05;
dp = 0:0.2:1;
for i = 1:length(dp)
    X = AmericanPut(K,r,D0,sigma,dp(i));
    if i == 1
        V = X;
    else
        V = [V;X(2,:)];
    end
end

AmericanTime = figure(7);

for i = 2:length(V(:,1))
    hold on
    plot(V(1,:),V(i,:));
end


axis([65 150 0 40]);
legend('T = 0','T = 0.2','T = 0.4','T = 0.6','T = 0.8','T = 1');
xlabel('Aktienkurs S');
ylabel('Optionspreis P');


if (1==0)
    print(AmericanTime,'-djpeg','-r350',['/Users/lucasburger/Dropbox/Studium/Bachelorarbeit/Bilder/AmerikanischeZeitverlauf.jpg']);
end





