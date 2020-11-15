clear;
% Sta³e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
timespan = 3050;
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
warpoczv1=0; %warunki pocz¹tkowe
warpoczv2=0;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;
while t<timespan %wykonuj na przedziale [0,15]
    k11=dv1dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2); %obliczanie wspó³czynników k dla obu zmiennych 
    k12=dv2dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2);
    k21=dv1dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2);
    k22=dv2dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2);
    k31=dv1dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2);
    k32=dv2dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2);
    k41=dv1dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2);
    k42=dv2dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2);
    V1(t+1)=V1(t)+1/6*h*(k11+2*k21+2*k31+k41); % wyznaczanie kolejnych wartoœci x1,y1
    V2(t+1)=V2(t)+1/6*h*(k12+2*k22+2*k32+k42);    
    t=t+h;
end
%rysowanie wykresu i liczenie œrednich b³êdów
plot(V1);
hold on
plot(V2);
legend('V1','V2');

figure;
plot(nthroot(V2/C2, 3));
hold on
plot(nthroot(V1/C1, 2));
legend('h2','h1');
ylabel('wysokoœæ s³upa wody [cm]');
xlabel('czas symulacji [s]');