clear all;
%% Sta³e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
h1zero = 28.8;
h2zero = 28.8;
h1lin = 28.63;
h2lin = 28.63;
F1zero = 98.5;
timespan = 10000;
tau = 50;
% I
Fd = 14.2; %dop³yw zak³ócaj¹cy
h=1;  %ustawienie kroku
% starttime=1; %pocz¹tek przedzia³u%
t=tau+1;

F1ster = zeros(1,timespan)+F1zero;
warpoczv1=C1*h1zero^2; %warunki pocz¹tkowe
warpoczv2=C2*h2zero^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;

V1lin(1:timespan)=warpoczv1;
V2lin(1:timespan)=warpoczv2;

skok = 9;
delta_u = [skok*3, skok*2, skok*1, 0, -skok*1, -skok*2, -skok*3];
h2_out = zeros(timespan,length(delta_u));
h2lin_out = zeros(timespan,length(delta_u));

%% Skoki sterowania
figure
hold on
title('Sterowanie F1_i_n')
xlabel('t [s]')
ylabel('F1_i_n [cm^3/s]')
for i=1:length(delta_u)
    F1ster(200:timespan)=F1zero + delta_u(i);
    
    V1(1:timespan)=warpoczv1;
    V2(1:timespan)=warpoczv2;
    
    while t<timespan %wykonuj na przedziale [0,15]
        F1 = F1ster(t-tau);
        [V1(1+t), V2(1+t)] = object(t,h,V1,V2,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2);
        
        t=t+h;
    end
    
    plot(F1ster,'DisplayName',sprintf('\\Delta u=%d', delta_u(i)));
    h2_out(:,i) = nthroot(V2/C2, 3); 
    t=tau+1;
end
hold off
legend


for i=1:length(delta_u)
    F1ster(200:timespan)=F1zero + delta_u(i);
   
    V1lin(1:timespan)=warpoczv1;
    V2lin(1:timespan)=warpoczv2;
    
    while t<timespan %wykonuj na przedziale [0,15]
        F1 = F1ster(t-tau);
        [V1lin(1+t), V2lin(1+t)] = objectLin(t,h,V1lin,V2lin,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
        t=t+h;
    end
    
    h2lin_out(:,i) = (V2lin/C2-h2lin^3)/(3*h2lin^2)+h2lin;%nthroot(V2lin/C2, 3);
    t=tau+1;
end


figure
title('Wyjœcie modelu liniowego')
xlabel('t [s]')
ylabel('Wysokoœæ s³upa cieczy h_2 [cm^3/s]')
hold on
for i=1:length(delta_u)
    plot(h2lin_out(:,i),'DisplayName',sprintf('\\Delta u=%d', delta_u(i)));
end
hold off
legend


figure
title('Wyjœcie modelu nieliniowego')
xlabel('t [s]')
ylabel('Wysokoœæ s³upa cieczy h_2 [cm^3/s]')
hold on
for i=1:length(delta_u)
    plot(h2_out(:,i),'DisplayName',sprintf('\\Delta u=%d', delta_u(i)));
end
hold off
legend