function s = fun_odp_skok(h2_0,delta_u,objectType,D)
%objectType - string 'liniowy', 'nieliniowy'
load('parameters.mat')
tau = 50;
Fd = 14.2; %dop³yw zak³ócaj¹cy
%Dlugosc symulacji
step_time = 2000;
timespan = step_time+D;


h1_0=h2_0;
F2=alfa1*h1_0^(1/2);


F1ster(1:step_time-2) = F2-Fd; %dop³yw wody do zbiornika
F1ster(step_time-1:timespan) = F2-Fd + delta_u;

t=tau+1;
warpoczv1=C1*h1_0^2; %warunki pocz¹tkowe
warpoczv2=C2*h2_0^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;


%% Liczenie modelu
while t<timespan %wykonuj na przedziale [0,15]
    
    if strcmp(objectType,'nieliniowy')
        [V1(1+t), V2(1+t)] = object(t,h,V1,V2,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2);
    elseif strcmp(objectType,'liniowy')
        [V1(1+t), V2(t+1)] = objectLin(t,h,V1,V2,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2,h1_0,h2_0);
    end
    t=t+h;
end
h2=nthroot(V2/C2, 3);
Ypp = h2(step_time-2);
s=(h2(step_time:timespan) - Ypp)/delta_u;
