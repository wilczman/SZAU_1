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
timespan = 3000;
tau = 50;
% I
Fd = 0; %dop³yw zak³ócaj¹cy
F1in(1:timespan) = 0; %dop³yw wody do zbiornika - wielkoœæ steruj¹ca
F1ster(1:2000) = 0; %dop³yw wody do zbiornika 
F1ster(2000:timespan) = 98.5;
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

%% Liczenie modelu
while t<timespan %wykonuj na przedziale [0,15]
    F1 = F1ster(t-tau);
    k11=dv1dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2); %obliczanie wspó³czynników k dla obu zmiennych 
    k12=dv2dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2);
    k21=dv1dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2);
    k22=dv2dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2);
    k31=dv1dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2);
    k32=dv2dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2);
    k41=dv1dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2);
    k42=dv2dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2);
    V1(t+1)=V1(t)+1/6*h*(k11+2*k21+2*k31+k41);
    if V1(t+1)<0
        V1(t+1) = 0;
    end
    V2(t+1)=V2(t)+1/6*h*(k12+2*k22+2*k32+k42);   
    if V2(t+1)<0
        V2(t+1) = 0;
    end    
    t=t+h;
end
t=tau+1;
V1lin(1:timespan)=warpoczv1;
V2lin(1:timespan)=warpoczv2;
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

%% Wykresy
% %rysowanie wykresu i liczenie œrednich b³êdów
% plot(V1);
% hold on
% plot(V2);
% legend('V1','V2');

figure;
plot(nthroot(V2lin/C2, 3));
% plot( (V2/C2-h2zero^3)/(3*h2zero^2)+h2zero);
hold on
plot(nthroot(V2/C2, 3));
legend('h2lin','h2');
ylabel('wysokoœæ s³upa wody [cm]');
xlabel('czas symulacji [s]');
title(sprintf('F1 = %.2f, Fd = %.2f, h1start = %.2f, h2start = %.2f',F1,Fd,h1zero,h2zero));

figure;
plot(nthroot(V1lin/C1, 2))
% plot( (V1/C1-h1zero^2)/(2*h1zero)+h1zero );
hold on
plot(nthroot(V1/C1, 2));
legend('h1lin','h1');
ylabel('wysokoœæ s³upa wody [cm]');
xlabel('czas symulacji [s]');
title(sprintf('F1 = %.2f, Fd = %.2f, h1start = %.2f, h2start = %.2f',F1,Fd,h1zero,h2zero));

