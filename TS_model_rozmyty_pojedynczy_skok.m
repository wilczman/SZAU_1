clear all;
%% Sta³e
load('parameters.mat');
timespan = 10000;
liczba_regulatorow=2;
% I
Fd = 14.2; %dop³yw zak³ócaj¹cy

F1ster(1:2000) = 98.5; %dop³yw wody do zbiornika
delta_u = 15;
F1ster(2000:timespan) = 98.5 + delta_u;
t=tau+2;
warpoczv1=C1*h1lin^2; %warunki pocz¹tkowe
warpoczv2=C2*h2lin^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;
h2roz = nthroot(V2/C2, 3);

h2_pocz=15;
h2_koniec=45; %[h2_pocz, h2_koniec] - zakres w jakim dziala model rozmyty
centra=linspace(h2_pocz,h2_koniec,liczba_regulatorow);
%% Model rozmyty
while t<timespan %wykonuj na przedziale [0,15]
    [V1(t), V2(t)] = objectLin(t-1,h,V1,V2,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    h2roz(t) = (V2lin(t)/C2-h2lin^3)/(3*h2lin^2)+h2lin; 
    for nr = 1:liczba_regulatorow
        w(nr) = gaussmf(h2roz(t), [gausy(liczba_regulatorow, h2_pocz, h2_koniec) centra(nr)]);
    end
    h2roz_suma=0;
    for nr=1:liczba_regulatorow 
        h2roz_suma = h2roz_suma+h2roz(t) * w(nr); 
    end
    h2roz(t) = h2roz_suma / sum(w);
    V2(t)=C2*h2roz(t)^3;
    %h1=((V2(t)+alfa2*(h2roz(t))^(1/2))/alfa1)^2;
    %V1(t)=C1*h1^2;
    t=t+h;
end

%% Model zwykly
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;
h2 = nthroot(V2/C2, 3);
t=tau+2;
while t<=timespan %wykonuj na przedziale [0,15]
    
    [V1(t), V2(t)] = object(t-1,h,V1,V2,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2);
    h2(t) = nthroot(V2(t)/C2, 3);
    t=t+h;
end
%% Pobranie odpowiedzi skokowej z obiektu liniowego

%% Wykresy


figure;
plot(h2);
% plot( (V2/C2-h2zero^3)/(3*h2zero^2)+h2zero);
hold on
plot(h2roz,'--');
legend('h2', 'h2_r_o_z');
ylabel('wysokoœæ s³upa wody [cm]');
xlabel('czas symulacji [s]');
%title('Model normalny');

figure;
plot(h2roz);
% plot( (V2/C2-h2zero^3)/(3*h2zero^2)+h2zero);
hold on
legend('h2_r_o_z');
ylabel('wysokoœæ s³upa wody [cm]');
xlabel('czas symulacji [s]');
title('Model rozmyty');