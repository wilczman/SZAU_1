clear;
% Sta�e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
timespan = 13000;
tau = 50;
% C1 = 0.5;
% C2 = 0.4;
% alfa1 = 18;
% alfa2 = 24;
% I
F1in(1:timespan) = 0; %dop�yw wody do zbiornika - wielko�� steruj�ca
F1 = 98.5; %dop�yw wody do zbiornika 
Fd = 14.2; %dop�yw zak��caj�cy
% F1 = 125; %dop�yw wody do zbiornika 
% Fd = 11; %dop�yw zak��caj�cy
% F2 = 0; %dop�yw do drugiego zbiornika = wyp�yw z pierwszego

% dV1dt = F1 + Fd - alfa1*(V1/C1)^(1/4);
% dV2dt = alfa1*(V1/C1)^(1/4) - alfa2*(V2/C2)^(1/6);

h=1;  %ustawienie kroku
% starttime=1; %pocz�tek przedzia�u%
t=tau;
warpoczv1=0; %warunki pocz�tkowe
warpoczv2=0;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;
while t<timespan %wykonuj na przedziale [0,15]
    if t == 4001
        Fd = 30;
    end
    k11=dv1dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2); %obliczanie wsp�czynnik�w k dla obu zmiennych 
    k12=dv2dt(V1(t),V2(t),F1,Fd,alfa1,alfa2,C1,C2);
    k21=dv1dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2);
    k22=dv2dt(V1(t)+0.5*h*k11,V2(t)+0.5*h*k12,F1,Fd,alfa1,alfa2,C1,C2);
    k31=dv1dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2);
    k32=dv2dt(V1(t)+0.5*h*k21,V2(t)+0.5*h*k22,F1,Fd,alfa1,alfa2,C1,C2);
    k41=dv1dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2);
    k42=dv2dt(V1(t)+h*k31,V2(t)+h*k32,F1,Fd,alfa1,alfa2,C1,C2);
    V1(t+1)=V1(t)+1/6*h*(k11+2*k21+2*k31+k41); % wyznaczanie kolejnych warto�ci x1,y1
    V2(t+1)=V2(t)+1/6*h*(k12+2*k22+2*k32+k42);    
    t=t+h;
end
%rysowanie wykresu i liczenie �rednich b��d�w
plot(V1);
hold on
plot(V2);
legend('V1','V2');

figure;
plot(nthroot(V2/C2, 3));
hold on
plot(nthroot(V1/C1, 2));
legend('h2','h1');
ylabel('wysoko�� s�upa wody [cm]');
xlabel('czas symulacji [s]');
h2 = nthroot(V2/C2, 3);
odp_skok_fd = (h2(4001:13000)-h2(4000) )/15.8;
figure;
plot(odp_skok_fd)