clear all

%%%%%%%    inicjalizacja    %%%%%%
load('parameters.mat')
liczba_regulatorow = 3; %ilo�� rozmytych
timespan=3000;
F1ster(1:timespan) = 98.5;
Fd = 14.2;
%[h2_pocz, h2_koniec] - zakres w jakim dziala model rozmyty
h2_pocz=15;
h2_koniec=45; 
centra = linspace(h2_pocz,h2_koniec,liczba_regulatorow);

%warunki pocz�tkowe
warpoczv1=C1*h1lin^2; 
warpoczv2=C2*h2lin^3;
V1roz(1:timespan)=warpoczv1;
V2roz(1:timespan)=warpoczv2;

%parametry DMC
D=5000; N=200; Nu=50; lambda = 1;
params=[D N Nu];
params_DMC{1} = params;
params_DMC{2} = [ params;  params];
params_DMC{3} = [ params; params; params ];
params_DMC{4} = [ params; params; params; params ] ;
params_DMC{5} = [ params; params ;params;params; params ];

%%%%%%%%Punkt Pracy%%%%%%%%
Upp=98.5;
Ypp=0%h2lin;

%%%%%%%%%Ograniczenia%%%%%%%
%u_min=-1-Upp;
%u_max=1-Upp;
delta_u_max=10;
delta_u_min=-delta_u_max;

%%%%%%%%%-----Parametry DMC------%%%%%%%
 lambda=[5,5,5,5,5]; %parametr lambda np. 1
% D=[50,50,50,50,50]; %horyzont dynamiki (D)
% N=[50,50,50,50,50];%horyzont predykcji (N)
% Nu=[30,30,30,30,30]; %horyzont sterowania (Nu)(ilosc przyszlych przyrostow wartosci sterowania)

%load('optymalne_parametry_DMC.mat');
%parametry dobrane przez funkcj� fmincon()
%D=150; N=nastawy_DMC_fmincon(1); Nu=nastawy_DMC_fmincon(2); lambda=nastawy_DMC_fmincon(3);
%D=150; N=35; Nu=5; lambda=1 %NAJLEPSZE EKSPERYMENTALNE PARAMETRY SWIATA
%D=200; N=50; Nu=10; % robocze parametry


%%%%%%deklaracja wektor�w sygna��w oraz b��d�w%%%%%%
for nr=1:liczba_regulatorow
    
    D(nr)=params_DMC{liczba_regulatorow}(nr,1); 
    N(nr)=params_DMC{liczba_regulatorow}(nr,2);
    Nu(nr)=params_DMC{liczba_regulatorow}(nr,3);  %parametry dobrane funkcj� ga

end

Y = zeros(1,timespan) + h2lin;
y = zeros(1,timespan)+ h2lin;

%Y(1:51)=Ypp;
U=zeros(1, timespan);  
U(:)=Upp;
e=zeros(1, timespan);
u_final=zeros(1,timespan);


yzad=zeros(1, timespan)+28.63;
yzad(round(1*timespan/6):round(1*timespan/3))=30;
yzad(round(1*timespan/3):round(2*timespan/3))=26;
yzad(round(2*timespan/3):round(3*timespan/3))=27;
yzad=yzad-Ypp;

% yzad=zeros(1, timespan)+28.63;
% 
% yzad(round(1*timespan/7):round(2*timespan/7))=26;
% yzad(round(2*timespan/7):round(3*timespan/7))=29;
% yzad(round(3*timespan/7):round(4*timespan/7))=30;
% yzad(round(4*timespan/7):round(5*timespan/7))=27;
% yzad(round(5*timespan/7):round(7*timespan/7))=25;
% yzad=yzad-Ypp;

%%%%%%%   tycie odpowiedzi skokowe  %%%%%%% 
for nr=1:liczba_regulatorow
    u{nr}=zeros(1,timespan);
    %odp_skok{nr} = gotowa_odp_skokowa;
    odp_skok{nr} = fun_odp_skok(centra(nr),2,'nieliniowy',max(D)+max(N));
end

%%%%%%% Macierze M, K, Mp %%%%%%%
for nr=1:liczba_regulatorow
    Mp{nr}=zeros(N(nr),D(nr)-1);        %macierz ma wymiary Nx(D-1)
    for i=1:D(nr)-1 %wypelnianie macierzy Mp
       Mp{nr}(1:N(nr), i)=odp_skok{nr}(i+1:N(nr)+i)-odp_skok{nr}(i);
    end
    %macierz wsp�czynnik�w odpowiedzi skokowej wymiary(NxNu)
    M=zeros(N(nr), Nu(nr));  
    i=0;
    for j=1:Nu(nr)  %wypelnianie macierzy trojkatnej dolnej M
       M(j:N(nr),j)=odp_skok{nr}(1:N(nr)-i).';  
       i=i+1;
    end

    I=eye(Nu(nr));  %tworzenie macierzy jednostkowej o wymiarach NuxNu
    K{nr}=inv(M.'*M+lambda(nr)*I)*M.';   %macierz K

    deltaUP{nr}(1:D(nr)-1,1)=0;
    deltaU{nr}=0;
end

%%%%%%%%% Algorytm DMC %%%%%%%%%
for t=tau+2:timespan-max(N) 
    %symulacja obiektu i regulatora
    [V1roz(t), V2roz(t)] = object(t-1,h,V1roz,V2roz,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2);
    Y(t) = nthroot(V2roz(t)/C2, 3);%(V2roz(t)/C2-h2lin^3)/(3*h2lin^2)+h2lin;
    y(t)=Y(t)-Ypp;
    e(t)=yzad(t)-y(t);
    
%     for nr = 1:liczba_regulatorow
%         h2lin_mod_lok = centra(nr);
%         h1lin_mod_lok = h2lin_mod_lok;
%         [V1roz(t), V2roz(t)] = objectLin(t-1,h,V1roz,V2roz,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2,h1lin_mod_lok,h2lin_mod_lok);
%         h2roz(t) = (V2roz(t)/C2-h2lin_mod_lok^3)/(3*h2lin_mod_lok^2)+h2lin_mod_lok;
%         %h2_mod_lok(nr)=h2roz(t);
%         V2_mod_lok(nr)=V2roz(t);
%         %V1_mod_lok(nr)=V1roz(t);
%         w(nr) = gaussmf(h2roz(t), [gausy(liczba_regulatorow, h2_pocz, h2_koniec) centra(nr)]);
%     end
    
    for nr=1:liczba_regulatorow
        deltaUP{nr}(2:D(nr)-1)=deltaUP{nr}(1:D(nr)-2);
        deltaUP{nr}(1) = u_final(t-1)-u_final(t-2);
        Y0{nr}=Mp{nr}*deltaUP{nr}+y(t);
        Yzad=yzad(t+1:t+N(nr));
        deltaU{nr}=K{nr}*(Yzad-Y0{nr});	
        delta_u{nr}=deltaU{nr}(1);

        u{nr}(t)=u_final(t-1)+delta_u{nr};
    end
    
    %%%%%%%------ROZMYCIE------%%%%%%%%
    for nr = 1:liczba_regulatorow
        w(nr) = gaussmf(y(t), [gausy(liczba_regulatorow, h2_pocz, h2_koniec) centra(nr)]);
    end
    
    %�rednia wa�ona z f.przynale�no�ci
    for nr=1:liczba_regulatorow 
        u_final(t) = u_final(t) + u{nr}(t) * w(nr); 
    end
    u_final(t) = u_final(t) / sum(w);
    
    delta=u_final(t)-u_final(t-1);
    if delta>delta_u_max
         delta=delta_u_max; 
    elseif delta<(-delta_u_max)
         delta=-delta_u_max;
    end
    
    u_final(t)=u_final(t-1)+delta;
    %ograniczenie sygna�u steruj�cego
%     if u_final(t)>u_max
%         u_final(t)=u_max;
%     elseif u_final(t)<u_min
%         u_final(t)=u_min;
%     end
     U(t)=u_final(t)+Upp;
     F1ster(t) = U(t);
end
wskaznik_jakosci=sum(e.^2);

% for nr=1:liczba_regulatorow
%     wskaznik_jakosci=sum(e.^2);
%     yzad=yzad(1:timespan-N(nr))+Ypp;
%     Y=Y(1:timespan-N(nr));
%     U=U(1:timespan-N(nr));
% end

    %%%%%%%%prezentacja wynik�w symulacji%%%%%%%%
figure;
stairs(U);hold on; xlim([0, timespan-N(nr)]);
title('Sygna� sterowania DMC'); xlabel('t');ylabel('warto�� sygna�u');

figure;
stairs(Y); hold on; 
stairs(yzad,':'); xlim([0, timespan-N(nr)]); 
title('Wyj�cie regulatora DMC'); xlabel('t'); ylabel('warto�� sygna�u');
legend('wyj�cie y(t)','warto�� zadana');
