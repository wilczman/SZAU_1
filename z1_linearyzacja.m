clear;
%% Sta³e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
h1zero = 28.63;
h2zero = 28.63;
h1lin = 28.63;
h2lin = 28.63;
timespan = 7000;
tau = 50;

%% I
F1in(1:timespan) = 0; %dop³yw wody do zbiornika - wielkoœæ steruj¹ca
F1ster(1:2000) = 98.5; %dop³yw wody do zbiornika 
F1ster(2000:timespan) = 98.5+10;
Fd = 14.2;%14.2; %dop³yw zak³ócaj¹cy
h=1;  %ustawienie kroku
t=tau+1;
% starttime=1; %pocz¹tek przedzia³u%
warpoczv1=C1*h1zero^2; %warunki pocz¹tkowe
warpoczv2=C2*h2zero^3;
V1lin(1:timespan)=warpoczv1;
V2lin(1:timespan)=warpoczv2;

%% symulacja
while t<timespan %wykonuj na przedziale [0,15]
    F1 = F1ster(t-tau);
    k11=lin_dv1dt(V1lin(t),V2lin(t),F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin); %obliczanie wspó³czynników k dla obu zmiennych 
    k12=lin_dv2dt(V1lin(t),V2lin(t),F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k21=lin_dv1dt(V1lin(t)+0.5*h*k11,V2lin(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k22=lin_dv2dt(V1lin(t)+0.5*h*k11,V2lin(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k31=lin_dv1dt(V1lin(t)+0.5*h*k21,V2lin(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k32=lin_dv2dt(V1lin(t)+0.5*h*k21,V2lin(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k41=lin_dv1dt(V1lin(t)+h*k31,V2lin(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    k42=lin_dv2dt(V1lin(t)+h*k31,V2lin(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
    V1lin(t+1)=V1lin(t)+1/6*h*(k11+2*k21+2*k31+k41);
    if V1lin(t+1)<0
        V1lin(t+1) = 0;
    end
    V2lin(t+1)=V2lin(t)+1/6*h*(k12+2*k22+2*k32+k42);   
    if V2lin(t+1)<0
        V2lin(t+1) = 0;
    end
    t=t+h;
end


%% rysowanie wykresu
plot(V1lin);
title('wykres V1lin,V2lin')
hold on
plot(V2lin);
legend('V1lin','V2lin');

figure;
title('wykres h1,h2')
plot(nthroot(V1lin/C1, 2));
hold on
plot(nthroot(V2lin/C2, 3));
legend('h1','h2');

% figure
% odp_skok = nthroot(V2lin/C2, 3);
% odp_skok = odp_skok(2000:timespan-1)
% odp_skok_scaled = (odp_skok-28.77)/10
% plot(odp_skok_scaled)
% save odp_skokH28.mat