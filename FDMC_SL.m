clear all

%%%%%%%    inicjalizacja    %%%%%%
load('parameters.mat')
liczba_regulatorow = 3; %iloœæ rozmytych
timespan=3000;
F1ster(1:timespan) = 98.5;
Fd = 14.2;
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
D=5000; N=200; Nu=50; lambda = 1;

%%%%%%%%Punkt Pracy%%%%%%%%
Upp=98.5;
Ypp=0;%h2lin;

%%%%%%%%%Ograniczenia%%%%%%%
u_min=70;%-1-Upp;

u_max=120;%1-Upp;

delta_u_max=10;
delta_u_min=-delta_u_max;
lb=ones(Nu,1)*delta_u_min;
ub=ones(Nu,1)*delta_u_max;
ymin=15;
ymax=45;


Y = zeros(1,timespan) + h2lin;
y = zeros(1,timespan)+ h2lin;

%Y(1:51)=Ypp;
U=zeros(1, timespan);  
U(:)=Upp;
u=zeros(1, timespan);
e=zeros(1, timespan);


% yzad=zeros(1, timespan)+28.63;
% yzad(round(1*timespan/6):round(2*timespan/6))=26;
% yzad(round(2*timespan/6):round(3*timespan/6))=29;
% yzad(round(3*timespan/6):round(4*timespan/6))=30;
% yzad(round(4*timespan/6):round(5*timespan/6))=27;
% yzad(round(5*timespan/6):round(6*timespan/6))=25;
%do krotszej symulacji
yzad=zeros(1, timespan)+28.63;
yzad(round(1*timespan/6):round(1*timespan/3))=30;
yzad(round(1*timespan/3):round(2*timespan/3))=26;
yzad(round(2*timespan/3):round(3*timespan/3))=27;
yzad=yzad-Ypp;

deltaUP(1:D-1,1)=0;
deltaU=zeros(Nu,1);
w=zeros(1,3);
%%%%%%%%% Algorytm DMC %%%%%%%%%
for t=tau+2:timespan-max(N)
    [V1roz(t), V2roz(t)] = object(t-1,h,V1roz,V2roz,F1ster(t-1-tau),Fd,alfa1,alfa2,C1,C2);
    Y(t) = nthroot(V2roz(t)/C2, 3);%(V2roz(t)/C2-h2lin^3)/(3*h2lin^2)+h2lin;
    y(t)=Y(t)-Ypp;
    e(t)=yzad(t)-y(t);
    
    %%%%%%%   tycie odpowiedzi skokowe  %%%%%%% 
    s_suma=0;
    for nr=1:liczba_regulatorow
        odp_skok{nr} = fun_odp_skok(centra(nr),2,'liniowy',y(t),D+N);
        w(nr) = gaussmf(y(t), [gausy(liczba_regulatorow, h2_pocz, h2_koniec) centra(nr)]);
        s_suma = s_suma + odp_skok{nr}*w(nr);
    end
    s_fala=s_suma/sum(w); %otrzymanie s z falka, czyli liniowa aproksymacja modelu nieliniowego
    
    
    Mp=zeros(N,D-1);        %macierz ma wymiary Nx(D-1)
    for i=1:D-1 %wypelnianie macierzy Mp
       Mp(1:N, i)=s_fala(i+1:N+i)-s_fala(i);
    end
    %macierz wspó³czynników odpowiedzi skokowej wymiary(NxNu)
    M=zeros(N, Nu);  
    i=0;
    for j=1:Nu  %wypelnianie macierzy trojkatnej dolnej M
       M(j:N,j)=s_fala(1:N-i).';  
       i=i+1;
    end
    
    deltaUP(2:D-1)=deltaUP(1:D-2);
    deltaUP(1) = u(t-1)-u(t-2);  
    Y0=Mp*deltaUP+y(t);
    Yzad=yzad(t)*ones(N,1); %yzad(t+1:t+N)';
    fun=@(delta_U) (Yzad-Y0-M*delta_U)'*(Yzad-Y0-M*delta_U)+lambda*delta_U'*delta_U;
    deltaU=fmincon(fun, deltaU,[],[],[],[],lb,ub);
    %I=eye(Nu);
    %deltaU=inv(M'*M+lambda*I)*M'*(Yzad-Y0);	
    delta_u=deltaU(1);
    
    u(t)=u(t-1)+delta_u;
    
    U(t)=u(t)+Upp;
    F1ster(t) = U(t);   
end
wskaznik_jakosci=sum(e.^2);

% for nr=1:liczba_regulatorow
%     wskaznik_jakosci=sum(e.^2);
%     yzad=yzad(1:timespan-N(nr))+Ypp;
%     Y=Y(1:timespan-N(nr));
%     U=U(1:timespan-N(nr));
% end

    %%%%%%%%prezentacja wyników symulacji%%%%%%%%
figure;
stairs(U);hold on; xlim([0, timespan-N]);
title('Sygna³ sterowania DMC'); xlabel('t');ylabel('wartoœæ sygna³u');

figure;
stairs(Y); hold on; 
stairs(yzad,':'); xlim([0, timespan-N]); 
title('Wyjœcie regulatora DMC'); xlabel('t'); ylabel('wartoœæ sygna³u');
legend('wyjœcie y(t)','wartoœæ zadana');

