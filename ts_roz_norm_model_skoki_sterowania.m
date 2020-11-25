clear all;
%% Sta�e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
h1zero = 28.8;
h2zero = 28.8;
h1lin = 28.63;
h2lin = 28.63;
F1zero = 98.5;
timespan = 5000;
tau = 50;
% I
Fd = 14.2; %dop�yw zak��caj�cy
h=1;  %ustawienie kroku
% starttime=1; %pocz�tek przedzia�u%
t=tau+1;

F1ster = zeros(1,timespan)+F1zero;
warpoczv1=C1*h1zero^2; %warunki pocz�tkowe
warpoczv2=C2*h2zero^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;

V1lin(1:timespan)=warpoczv1;
V2lin(1:timespan)=warpoczv2;

skok = 10;
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
    
    V1lin(1:timespan)=warpoczv1;
    V2lin(1:timespan)=warpoczv2;
    
    while t<timespan %wykonuj na przedziale [0,15]
        F1 = F1ster(t-tau);
        [V1(1+t), V2(1+t)] = object(t,h,V1,V2,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2);
        [V1lin(1+t), V2lin(1+t)] = objectLin(t,h,V1lin,V2lin,F1ster(t-tau),Fd,alfa1,alfa2,C1,C2,h1lin,h2lin);
        
        t=t+h;
    end
    
    plot(F1ster,'DisplayName',sprintf('\\Delta u=%d', delta_u(i)));
    h2_out(:,i) = nthroot(V2/C2, 3);
    h2lin_out(:,i) = (V2lin/C2-h2lin^3)/(3*h2lin^2)+h2lin;%nthroot(V2lin/C2, 3);
    
    t=tau+1;
end
hold off
legend

%%%%%---Model rozmyty---%%%%%
liczba_regulatorow=5;
warpoczv1=C1*h1lin^2; %warunki pocz�tkowe
warpoczv2=C2*h2lin^3;
V1roz(1:timespan)=warpoczv1;
V2roz(1:timespan)=warpoczv2;
h2roz(1:timespan) = h2lin;

h2_pocz=15;
h2_koniec=45; %[h2_pocz, h2_koniec] - zakres w jakim dziala model rozmyty
centra=linspace(h2_pocz,h2_koniec,liczba_regulatorow);
t=tau+2;
for i=1:length(delta_u)
    F1ster(200:timespan)=F1zero + delta_u(i);
    d_u=delta_u(i);
    V1lin(1:timespan)=warpoczv1;
    V2lin(1:timespan)=warpoczv2;
    while t<timespan %wykonuj na przedziale [0,15]
    
    for nr = 1:liczba_regulatorow
        h2lin_mod_lok = centra(nr);
        h1lin_mod_lok = h2lin_mod_lok;
        [V1roz(t), V2roz(t)] = objectLin(t-1,h,V1roz,V2roz,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2,h1lin_mod_lok,h2lin_mod_lok);
        h2roz(t) = (V2roz(t)/C2-h2lin_mod_lok^3)/(3*h2lin_mod_lok^2)+h2lin_mod_lok;
        %h2_mod_lok(nr)=h2roz(t);
        V2_mod_lok(nr)=V2roz(t);
        %V1_mod_lok(nr)=V1roz(t);
        w(nr) = gaussmf(h2roz(t), [gausy(liczba_regulatorow, h2_pocz, h2_koniec) centra(nr)]);
    end    
    V2roz_suma=0;
    for nr=1:liczba_regulatorow 
        V2roz_suma = V2roz_suma + V2_mod_lok(nr)*w(nr);  
    end
    V2roz(t)=V2roz_suma/sum(w);
    h2roz(t)=nthroot(V2roz(t)/C2, 3);%(V2roz(t)/C2-h2lin^3)/(3*h2lin^2)+h2lin;
    
    V1roz(t)=V1roz(t-1);
    [V1roz(t), xx] = object(t,h,V1roz,V2roz,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2);
    %V1    
    t=t+h;
    end
    h2roz_out(:,i) =nthroot(V2roz/C2, 3); %(V2roz/C2-h2lin^3)/(3*h2lin^2)+h2lin;%nthroot(V2lin/C2, 3);
    t=tau+2;
end
figure
title('Wyj�cie modelu liniowego')
xlabel('t [s]')
ylabel('Wysoko�� s�upa cieczy h_2 [cm^3/s]')
hold on
for i=1:length(delta_u)
    plot(h2lin_out(:,i),'DisplayName',sprintf('\\Delta u=%d', delta_u(i)));
end
hold off
legend


figure
title('Wyj�cie modelu nieliniowego')
xlabel('t [s]')
ylabel('Wysoko�� s�upa cieczy h_2 [cm^3/s]')
hold on
for i=1:length(delta_u)
    plot(h2_out(:,i),'DisplayName',sprintf('\\Delta u=%d', delta_u(i)));
end
hold off
legend

figure
title('Wyj�cie modelu rozmytego')
xlabel('t [s]')
ylabel('Wysoko�� s�upa cieczy h_2 [cm^3/s]')
hold on
for i=1:length(delta_u)
    plot(h2roz_out(:,i),'DisplayName',sprintf('\\Delta u=%d', delta_u(i)));
end
hold off
legend