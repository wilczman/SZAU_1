clear;
% Sta³e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
h1zero = 28.63;
h2zero = 28.63;
timespan = 30000;
tau = 50;
% C1 = 0.5;
% C2 = 0.4;
% alfa1 = 18;
% alfa2 = 24;
% I
F1in(1:timespan) = 0; %dop³yw wody do zbiornika - wielkoœæ steruj¹ca
F1 = 98.5; %dop³yw wody do zbiornika 
Fd = 14.2; %dop³yw zak³ócaj¹cy
% F1 = 125; %dop³yw wody do zbiornika 
% Fd = 11; %dop³yw zak³ócaj¹cy
% F2 = 0; %dop³yw do drugiego zbiornika = wyp³yw z pierwszego

% dV1dt = F1 + Fd - alfa1*(V1/C1)^(1/4);
% dV2dt = alfa1*(V1/C1)^(1/4) - alfa2*(V2/C2)^(1/6);

h=1;  %ustawienie kroku
% starttime=1; %pocz¹tek przedzia³u%
t=tau;
warpoczv1=C1*28.8^2; %warunki pocz¹tkowe
warpoczv2=C2*26.3^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;
while t<timespan %wykonuj na przedziale [0,15]
    k11=lin_dv1dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero); %obliczanie wspó³czynników k dla obu zmiennych 
    k12=lin_dv2dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero);
    k21=lin_dv1dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero);
    k22=lin_dv2dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero);
    k31=lin_dv1dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero);
    k32=lin_dv2dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero);
    k41=lin_dv1dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero);
    k42=lin_dv2dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2,h1zero,h2zero);
    V1(t+1)=V1(t)+1/6*h*(k11+2*k21+2*k31+k41); % wyznaczanie kolejnych wartoœci x1,y1
    V2(t+1)=V2(t)+1/6*h*(k12+2*k22+2*k32+k42);    
    t=t+h;
end
%rysowanie wykresu i liczenie œrednich b³êdów
plot(V1);
title('wykres V1,V2')
hold on
plot(V2);
legend('V1','V2');

figure;
title('wykres h1,h2')
plot(nthroot(V1/C1, 2));
hold on
plot(nthroot(V2/C2, 3));
legend('h2','h1');