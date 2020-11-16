clear all;
%% Sta³e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
h1zero = 28.63;
h2zero = 28.63;
h1lin = 28.63;
h2lin = 28.63;
timespan = 15000;
tau = 50;
% I
Fd = 14.2; %dop³yw zak³ócaj¹cy
F1in(1:timespan) = 0; %dop³yw wody do zbiornika - wielkoœæ steruj¹ca
F1ster(1:5998) = 98.5; %dop³yw wody do zbiornika
delta_u = 10;
F1ster(5999:timespan) = 98.5 + delta_u;
% F1 = 125; %dop³yw wody do zbiornika 
% Fd = 11; %dop³yw zak³ócaj¹cy
% F2 = 0; %dop³yw do drugiego zbiornika = wyp³yw z pierwszego

h=1;  %ustawienie kroku
% starttime=1; %pocz¹tek przedzia³u%
t=tau+1;
warpoczv1=C1*h1zero^2; %warunki pocz¹tkowe
warpoczv2=C2*h2zero^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;

V1lin(1:timespan)=warpoczv1;
V2lin(1:timespan)=warpoczv2;

%% Liczenie modelu
while t<timespan %wykonuj na przedziale [0,15]
    F1 = F1ster(t-tau);
    [V1(1+t), V2(1+t)] = object(t,h,V1,V2,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2);
    [V1lin(1+t), V2lin(1+t)] = objectLin(t,h,V1lin,V2lin,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    t=t+h;
end
h2=nthroot(V2lin/C2, 3);

%% Pobranie odpowiedzi skokowej z obiektu liniowego
Ypp = h2(5998);
s=(h2(6000:timespan) - Ypp)/delta_u;


%% Wykresy


figure
plot(h2);
title('Wysokosc h_2 podczas calej symulacji')
ylabel('wysokoœæ s³upa wody [cm]');
xlabel('czas symulacji [s]');


figure
plot(s);
title('OdpowiedŸ skokowa')

save('rw_odp_skok','s')