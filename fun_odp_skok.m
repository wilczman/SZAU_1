function s = fun_odp_skok(h2_0,delta_u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load('parameters.mat')
%Dlugosc symulacji
kk=10000;


timespan = 10000;
tau = 50;
% I
Fd = 14.2; %dop�yw zak��caj�cy

h1_0=h2_0;
F2=alfa1*h1_0^(1/2);
F1ster(1:1998) = F2-Fd; %dop�yw wody do zbiornika
F1ster(1999:timespan) = F2-Fd + delta_u;
t=tau+1;
warpoczv1=C1*h1_0^2; %warunki pocz�tkowe
warpoczv2=C2*h2_0^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;


%% Liczenie modelu
while t<timespan %wykonuj na przedziale [0,15]
    [V1(1+t), V2(1+t)] = object(t,h,V1,V2,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2);
    t=t+h;
end
h2=nthroot(V2/C2, 3);
Ypp = h2(1998);
s=(h2(2000:timespan) - Ypp)/delta_u;
