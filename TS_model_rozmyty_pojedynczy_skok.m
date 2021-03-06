clear all; close all;
%% Sta�e
load('parameters.mat');
timespan = 10000;
liczba_regulatorow=5;
% I
Fd = 14.2; %dop�yw zak��caj�cy

F1ster(1:2000) = 98.5; %dop�yw wody do zbiornika
delta_u = 15;
F1ster(2000:timespan) = 98.5 + delta_u;
t=tau+2;
warpoczv1=C1*h1lin^2; %warunki pocz�tkowe
warpoczv2=C2*h2lin^3;
V1roz(1:timespan)=warpoczv1;
V2roz(1:timespan)=warpoczv2;
h2roz(1:timespan) = h2lin;

h2_pocz=15;
h2_koniec=45; %[h2_pocz, h2_koniec] - zakres w jakim dziala model rozmyty
centra=linspace(h2_pocz,h2_koniec,liczba_regulatorow);
%% Model rozmyty
while t<timespan %wykonuj na przedziale [0,15]
    
    for nr = 1:liczba_regulatorow
        h2lin_mod_lok = centra(nr);
        h1lin_mod_lok = h2lin_mod_lok;
        [V1roz(t), V2roz(t)] = objectLin(t-1,h,V1roz,V2roz,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2,h1lin_mod_lok,h2lin_mod_lok);
        h2roz(t) = (V2roz(t)/C2-h2lin_mod_lok^3)/(3*h2lin_mod_lok^2)+h2lin_mod_lok;
        %h2_mod_lok(nr)=h2roz(t);
        V2_mod_lok(nr)=V2roz(t);
        %V1_mod_lok(nr)=V1roz(t);
        w(nr) = gaussmf(h2roz(t), [gausy(liczba_regulatorow, h2_pocz, h2_koniec) centra(nr)]);
    end

%     h2roz_suma=0;
%     for nr=1:liczba_regulatorow 
%         h2roz_suma = h2roz_suma + h2roz(t)*w(nr); 
%     end
%     h2roz(t) = h2roz_suma / sum(w);
    
%     h2roz_suma=0;
%     for nr=1:liczba_regulatorow 
%         h2roz_suma = h2roz_suma + h2_mod_lok(nr)*w(nr); 
%     end
%     h2roz(t) = h2roz_suma / sum(w);
    
    V2roz_suma=0;
    %V1roz_suma=0;
    for nr=1:liczba_regulatorow 
        V2roz_suma = V2roz_suma + V2_mod_lok(nr)*w(nr); 
       % V1roz_suma = V1roz_suma + V1_mod_lok(nr)*w(nr); 
    end
    V2roz(t)=V2roz_suma/sum(w);
   % V1roz(t)=V1roz_suma/sum(w);
    h2roz(t)=nthroot(V2roz(t)/C2, 3);%(V2roz(t)/C2-h2lin^3)/(3*h2lin^2)+h2lin;
    %h2roz(t) = h2roz_suma / sum(w);
    %V2(t)=C2*h2roz(t)^3;
    %h1=((V2roz(t)+alfa2*(h2roz(t))^(1/2))/alfa1)^2;
    %V1roz(t)=C1*h1^2;
    
    %h1_t_minus_1=((V2roz(t)-V2roz(t-1)+alfa2*(h2roz(t-1)))/alfa1)^(1/2);
    %V1roz(t)=V1roz(t-1)+F1ster(t-1-tau)+Fd-alfa1*(h1_t_minus_1)^(1/2);
    V1roz(t)=V1roz(t-1);
    [V1roz(t), xx] = object(t,h,V1roz,V2roz,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2);
    %V1    
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
ylabel('wysoko�� s�upa wody [cm]');
xlabel('czas symulacji [s]');
%title('Model normalny');

figure;
plot(V2);
% plot( (V2/C2-h2zero^3)/(3*h2zero^2)+h2zero);
hold on
plot(V2roz,'--');
legend('V2', 'V2_r_o_z');
ylabel('objetosc [cm]');
xlabel('czas symulacji [s]');

figure;
plot(V1);
% plot( (V2/C2-h2zero^3)/(3*h2zero^2)+h2zero);
hold on
plot(V1roz,'--');
legend('V1', 'V1_r_o_z');
ylabel('objetosc [cm]');
xlabel('czas symulacji [s]');

figure;
plot(h2roz);
% plot( (V2/C2-h2zero^3)/(3*h2zero^2)+h2zero);
hold on
legend('h2_r_o_z');
ylabel('wysoko�� s�upa wody [cm]');
xlabel('czas symulacji [s]');
title('Model rozmyty');