clear all

%%%%%%%    inicjalizacja    %%%%%%
load('parameters.mat')
load('TS_param.mat')
liczba_regulatorow = 5; %iloœæ rozmytych
timespan=20000;
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

%%%%%%%%Punkt Pracy%%%%%%%%
Upp=98.5;
Ypp=0;%h2lin;

%%%%%%%%%Ograniczenia%%%%%%%
u_min=0-Upp;
u_max=160-Upp;

delta_u_max=2;
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


yzad=zeros(1, timespan)+28.63;
yzad(round(1*timespan/7):round(2*timespan/7))=h2lin - 9;
yzad(round(2*timespan/7):round(3*timespan/7))=h2lin - 4;
yzad(round(3*timespan/7):round(4*timespan/7))=h2lin + 2;
yzad(round(4*timespan/7):round(5*timespan/7))=h2lin + 14;
yzad(round(5*timespan/7):round(6*timespan/7))=h2lin + 1;
yzad(round(6*timespan/7):round(7*timespan/7))=h2lin - 6;
Fd=zeros(1, timespan)+28.63;
Fd(round(1*timespan/14):round(2*timespan/7))= 12;
Fd(round(3*timespan/14):round(3*timespan/7))=15;
Fd(round(5*timespan/14):round(4*timespan/7))=12;
Fd(round(7*timespan/14):round(5*timespan/7))=18;
Fd(round(9*timespan/14):round(6*timespan/7))=14;
Fd(round(11*timespan/14):round(7*timespan/7))=12;

% %do krotszej symulacji
% yzad=zeros(1, timespan)+28.63;
% yzad(round(1*timespan/6):round(3*timespan/6))=28.63;
% yzad(round(3*timespan/6):round(5*timespan/6))=26;
% yzad(round(5*timespan/6):round(timespan)) = 27;
% yzad=yzad-Ypp;
% 
% Fd(round(1*timespan/6):round(4*timespan/6))=0;
% Fd(round(4*timespan/6):round(11*timespan/12))=15;
% Fd(round(11*timespan/12):round(timespan))=10;

deltaUP(1:D-1,1)=0;
deltaU=zeros(Nu,1);
w=zeros(1,3);

for nr=1:liczba_regulatorow
        s{nr} = fun_odp_skok(centra(nr),2,'liniowy', D+N);
end
%%%%%%%%% Algorytm DMC %%%%%%%%%
for t=tau+2:timespan-max(N)
    [V1roz(t), V2roz(t)] = object(t-1,h,V1roz,V2roz,F1ster(t-1-tau),Fd(t-1),alfa1,alfa2,C1,C2);
    Y(t) = nthroot(V2roz(t)/C2, 3);%(V2roz(t)/C2-h2lin^3)/(3*h2lin^2)+h2lin;
    y(t)=Y(t)-Ypp;
    e(t)=yzad(t)-y(t);
    
    %%%%%%%   liczenie wag  %%%%%%% 
    for nr=1:liczba_regulatorow
        w(nr) = gaussmf(y(t), [sigma(nr) centra(nr)]);
    end
    w_fala=w/sum(w);
    
    %otrzymanie s z falka
    s_fala=0;
    for nr=1:liczba_regulatorow
        s_fala = s_fala + s{nr}*w_fala(nr);
    end
    
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
    options=optimoptions('fmincon','Display','off');
    deltaU=fmincon(fun, deltaU,[],[],[],[],lb,ub,[],options);
    %I=eye(Nu);
    %deltaU=inv(M'*M+lambda*I)*M'*(Yzad-Y0);	
    delta_u=deltaU(1);
    
    %ograniczenie ró¿nic sygna³u steruj¹cego 
%     if delta_u>delta_u_max
%          delta_u=delta_u_max; 
%     elseif delta_u<(-delta_u_max)
%          delta_u=-delta_u_max;
%     end

    u(t)=u(t-1)+delta_u;
    %ograniczenie sygna³u steruj¹cego
    if u(t)>u_max
        u(t)=u_max;
    elseif u(t)<u_min
        u(t)=u_min;
    end
    
    U(t)=u(t)+Upp;
    F1ster(t) = U(t);   
end
wskaznik_jakosci=sum(e.^2);



    %%%%%%%%prezentacja wyników symulacji%%%%%%%%
figure;title('Regulator FDMC-SL - sterowanie');hold on
ylabel('F_1_i_n [cm^3/s]');
xlabel('t [s]');
xlim([1 timespan-max(N)])
plot(F1ster(1:timespan-max(N)));plot(Fd,'--');legend('U','Fd');hold off

str=sprintf('E=%f', wskaznik_jakosci);
figure;title({'Regulator FDMC-SL - wyjœcie obiektu',str});hold on;
ylabel('h_2 [cm]');
xlabel('t [s]');
xlim([1 timespan-max(N)])
plot(Y(1:timespan-max(N))-Ypp);plot(yzad(1:timespan-max(N)));legend('Y','yzad');hold off;

save('fdmc_sl_y.mat','Y','yzad','F1ster')