clear all;
%% Sta³e
load('parameters.mat')
timespan = 8500;
Fd = 0; %dop³yw zak³ócaj¹cy

%% Liczenie charakterystyki statycznej dla F1ster = [2*1,2*kk]
kk=600;
Ustat=zeros(kk,1);
Ystat=zeros(kk,1);
Ustat_lin=zeros(kk,1);
Ystat_lin=zeros(kk,1);
for i=1:kk
    V1(1:timespan)=0;
    V2(1:timespan)=0;

    V1lin(1:timespan)=0;
    V2lin(1:timespan)=0;
    
    F1ster(1:timespan) = i;
    t=tau+1;
    while t<timespan %wykonuj na przedziale [0,15]
        F1 = F1ster(t-tau);
        [V1(1+t), V2(1+t)] = object(t,h,V1,V2,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2);
        [V1lin(1+t), V2lin(1+t)] = objectLin(t,h,V1lin,V2lin,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
        t=t+h;
    end
    h2_out = nthroot(V2(timespan)/C2, 3);
    h2lin_out = (V2lin(timespan)/C2-h2lin^3)/(3*h2lin^2)+h2lin;
    
    Ustat(i)=F1ster(timespan);
    Ystat(i)=h2_out;
    
    Ustat_lin(i)=F1ster(timespan);
    Ystat_lin(i)=h2lin_out;
    
end

%WYKRES CHARAKTERYSTYKI STATYCZNEJ
figure
plot(Ustat,Ystat);
xlabel('U');
ylabel('Y');
title('Charakterystyka statyczna Y(U)');

figure
plot(Ustat_lin,Ystat_lin);
xlabel('U');
ylabel('Y');
title('Charakterystyka statyczna liniowa Y(U)');


