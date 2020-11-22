function s = fun_odp_skok(center,delta_u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load('parameters.mat')
%Dlugosc symulacji
kk=10000;


timespan = 10000;
tau = 50;
% I
Fd = 14.2; %dop³yw zak³ócaj¹cy
F1in(1:timespan) = 0; %dop³yw wody do zbiornika - wielkoœæ steruj¹ca
F1ster(1:1998) = center; %dop³yw wody do zbiornika
F1ster(1999:timespan) = center + delta_u;
t=tau+1;
warpoczv1=C1*h1zero^2; %warunki pocz¹tkowe
warpoczv2=C2*h2zero^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;
V1lin(1:timespan)=warpoczv1;
V2lin(1:timespan)=warpoczv2;

%% Liczenie modelu
while t<timespan %wykonuj na przedziale [0,15]
    [V1(1+t), V2(1+t)] = object(t,h,V1,V2,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2);
    t=t+h;
end
h2=nthroot(V2lin/C2, 3);
Ypp = h2(1998);
s=(h2(2000:timespan) - Ypp)/delta_u;
