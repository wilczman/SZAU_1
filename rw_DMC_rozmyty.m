clear all

%%%%%%%    inicjalizacja    %%%%%%
load('parameters.mat')
liczba_regulatorow = 5; %iloœæ rozmytych
timespan=17000;
F1ster(1:timespan) = 98.5;
Fd(1:timespan) = 14.2;
%[h2_pocz, h2_koniec] - zakres w jakim dziala model rozmyty
h2_pocz=15;
h2_koniec=45; 
centra = linspace(h2_pocz,h2_koniec,liczba_regulatorow);

%warunki pocz¹tkowe
warpoczv1=C1*h1lin^2; 
warpoczv2=C2*h2lin^3;
V1roz(1:timespan)=warpoczv1;
V2roz(1:timespan)=warpoczv2;

%parametry DMC
D=5000; N=200; Nu=50; lambda = 5;
params=[D N Nu];
params_DMC{1} = params;
params_DMC{2} = [ params;  params];
params_DMC{3} = [ params; params; params ];
params_DMC{4} = [ params; params; params; params ] ;
params_DMC{5} = [ params; params ;params;params; params ];
lambda=[lambda,lambda,lambda,lambda,lambda]; %parametr lambda np. 1
 
%%%%%%%%Punkt Pracy%%%%%%%%
Upp=98.5;
Ypp=0%h2lin;

%%%%%%%%%Ograniczenia%%%%%%%
u_min=0-Upp;
u_max=180-Upp;
delta_u_max=2;
delta_u_min=-delta_u_max;


%%%%%%deklaracja wektorów sygna³ów oraz b³êdów%%%%%%
for nr=1:liczba_regulatorow
    
    D(nr)=params_DMC{liczba_regulatorow}(nr,1); 
    N(nr)=params_DMC{liczba_regulatorow}(nr,2);
    Nu(nr)=params_DMC{liczba_regulatorow}(nr,3);  %parametry dobrane funkcj¹ ga

end

Y = zeros(1,timespan) + h2lin;
y = zeros(1,timespan)+ h2lin;

%Y(1:51)=Ypp;
U=zeros(1, timespan);  
U(:)=Upp;
e=zeros(1, timespan);
u_final=zeros(1,timespan);


% yzad=zeros(1, timespan)+28.63;
% yzad(round(1*timespan/6):round(1*timespan/3))=30;
% yzad(round(1*timespan/3):round(2*timespan/3))=26;
% yzad(round(2*timespan/3):round(3*timespan/3))=27;
% yzad=yzad-Ypp;

yzad=zeros(1, timespan)+28.63;
yzad(round(1*timespan/7):round(2*timespan/7))=h2lin - 9;
yzad(round(2*timespan/7):round(3*timespan/7))=h2lin - 4;
yzad(round(3*timespan/7):round(4*timespan/7))=h2lin + 2;
yzad(round(4*timespan/7):round(5*timespan/7))=h2lin + 14;
yzad(round(5*timespan/7):round(6*timespan/7))=h2lin + 1;
yzad(round(6*timespan/7):round(7*timespan/7))=h2lin - 6;
Fd=zeros(1, timespan);
Fd(round(1*timespan/14):round(2*timespan/7))=10;
Fd(round(3*timespan/14):round(3*timespan/7))=15;
Fd(round(5*timespan/14):round(4*timespan/7))=0;
Fd(round(7*timespan/14):round(5*timespan/7))=5;
Fd(round(9*timespan/14):round(6*timespan/7))=25;
Fd(round(12*timespan/14):round(7*timespan/7))=22.5;


%%%%%%%   tycie odpowiedzi skokowe  %%%%%%% 
for nr=1:liczba_regulatorow
    u{nr}=zeros(1,timespan);
    %odp_skok{nr} = gotowa_odp_skokowa;
    odp_skok{nr} = fun_odp_skok(centra(nr),2,'liniowy',max(D)+max(N));
    odp_skok_fd{nr} = fun_odp_skok_fd(centra(nr),2,'liniowy',max(D)+max(N));
end

%%%%%%% Macierze M, K, Mp %%%%%%%
for nr=1:liczba_regulatorow
    Mp{nr}=zeros(N(nr),D(nr)-1);        %macierz ma wymiary Nx(D-1)
    MZp{nr}=zeros(N(nr),D(nr)-1);
    for i=1:D(nr)-1 %wypelnianie macierzy Mp
       Mp{nr}(1:N(nr), i)=odp_skok{nr}(i+1:N(nr)+i)-odp_skok{nr}(i);
       MZp{nr}(1:N(nr), i)=odp_skok_fd{nr}(i+1:N(nr)+i)-odp_skok_fd{nr}(i);
    end
    %macierz wspó³czynników odpowiedzi skokowej wymiary(NxNu)
    M=zeros(N(nr), Nu(nr));  
    i=0;
    for j=1:Nu(nr)  %wypelnianie macierzy trojkatnej dolnej M
       M(j:N(nr),j)=odp_skok{nr}(1:N(nr)-i).';  
       i=i+1;
    end

    I=eye(Nu(nr));  %tworzenie macierzy jednostkowej o wymiarach NuxNu
    K{nr}=inv(M.'*M+lambda(nr)*I)*M.';   %macierz K
    KZ{nr}=K{nr}(1,:)*MZp{nr};

    deltaUP{nr}(1:D(nr)-1,1)=0;
    deltaU{nr}=0;
    deltaZp{nr}(1:D(nr)-1,1)=0;
end

%%%%%%%%% Algorytm DMC %%%%%%%%%
for t=tau+2:timespan-max(N) 
    %symulacja obiektu i regulatora
    [V1roz(t), V2roz(t)] = object(t-1,h,V1roz,V2roz,F1ster(t-1-tau),Fd(t-1-tau),alfa1,alfa2,C1,C2);
    Y(t) = nthroot(V2roz(t)/C2, 3);%(V2roz(t)/C2-h2lin^3)/(3*h2lin^2)+h2lin;
    y(t)=Y(t)-Ypp;
    e(t)=yzad(t)-y(t);
        
    for nr=1:liczba_regulatorow
        deltaUP{nr}(2:D(nr)-1)=deltaUP{nr}(1:D(nr)-2);
        deltaUP{nr}(1) = u_final(t-1)-u_final(t-2);
        
        deltaZp{nr}(2:D(nr)-1)=deltaZp{nr}(1:D(nr)-2);
        deltaZp{nr}(1) = Fd(t-1)-Fd(t-2);
        
        Y0{nr}=Mp{nr}*deltaUP{nr}+y(t);
        Yzad=yzad(t+1:t+N(nr));
        deltaU{nr}=K{nr}*(Yzad-Y0{nr});	
        delta_u{nr}=deltaU{nr}(1);

        u{nr}(t)=u_final(t-1)+delta_u{nr};%+KZ{nr}*deltaZp{nr};
    end
    
    %%%%%%%------ROZMYCIE------%%%%%%%%
    for nr = 1:liczba_regulatorow
        w(nr) = gaussmf(y(t), [gausy(liczba_regulatorow, h2_pocz, h2_koniec) centra(nr)]);
    end
    
    %œrednia wa¿ona z f.przynale¿noœci
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
    
    %ograniczenie sygna³u steruj¹cego
    if u_final(t)>u_max
        u_final(t)=u_max;
    elseif u_final(t)<u_min
        u_final(t)=u_min;
    end
    
    U(t)=u_final(t)+Upp;
    F1ster(t) = U(t);
end
wskaznik_jakosci=sum(e.^2)


%%%%%%%%prezentacja wyników symulacji%%%%%%%%
figure;title('Regulator DMC - sterowanie');hold on
ylabel('F_1_i_n [cm^3/s]');
xlabel('t [s]');
xlim([1 timespan-max(N)])
plot(F1ster(1:timespan-max(N)));plot(Fd,'--');legend('U','Fd')
hold off

str=sprintf('E=%f', wskaznik_jakosci);
figure;title({'Rozmyty regulator DMC - wyjœcie obiektu',str});hold on;
ylabel('h_2 [cm]');
xlabel('t [s]');
xlim([1 timespan-max(N)])
plot(Y(1:timespan-max(N))-Ypp);plot(yzad(1:timespan-max(N)));
legend('Y','yzad');hold off;

save('dmc_rozmyty_y.mat','Y','yzad','F1ster')