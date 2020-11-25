clear all;
%% Sta³e
C1 = 0.55;
C2 = 0.75;
alfa1 = 21;
alfa2 = 21;
h1zero = 24;
h2zero = 24;
h1lin = 28.63;
h2lin = 28.63;
timespan = 17000;
tau = 50;

%% I
Fd(1:timespan) = 0; %dop³yw zak³ócaj¹cy
F1ster(1:4000) = 98.5; %dop³yw wody do zbiornika
h=1;  %ustawienie kroku
t=tau+1;
% starttime=1; %pocz¹tek przedzia³u%
warpoczv1=C1*h1lin^2; %warunki pocz¹tkowe
warpoczv2=C2*h2lin^3;
V1(1:timespan)=warpoczv1;
V2(1:timespan)=warpoczv2;

%% DMC

load('rw_odp_skok.mat','s')
load('odp_skok_fd.mat')

final_odp_skok=s;
D=5000; N=200; Nu=50;lambda = 5;

Upp=98.5;
Ypp=0;%24;

%Ograniczenia
u_min=0-Upp;
u_max=180-Upp;
delta_u_max=2;
delta_u_min=-delta_u_max;

%%%%%%deklaracja wektorów sygna³ów oraz b³êdów%%%%%%
U=zeros(1, timespan);
u=zeros(timespan);
deltaUP(1:D-1,1)=0;
deltaU=0;
deltaZP(1:D-1,1)=0;

Y = zeros(1,timespan) + h2lin;
y = zeros(1,timespan)+ h2lin;
e = zeros(1, timespan);

yzad=zeros(1, timespan)+28.63;
yzad(round(1*timespan/7):round(2*timespan/7))=h2lin - 9;
yzad(round(2*timespan/7):round(3*timespan/7))=h2lin - 4;
yzad(round(3*timespan/7):round(4*timespan/7))=h2lin + 2;
yzad(round(4*timespan/7):round(5*timespan/7))=h2lin + 14;
yzad(round(5*timespan/7):round(6*timespan/7))=h2lin + 1;
yzad(round(6*timespan/7):round(7*timespan/7))=h2lin - 6;
Fd=zeros(1, timespan)+28.63;
Fd(round(1*timespan/14):round(2*timespan/7))= 10;
Fd(round(3*timespan/14):round(3*timespan/7))=15;
Fd(round(5*timespan/14):round(4*timespan/7))=5;
Fd(round(7*timespan/14):round(5*timespan/7))=25;
Fd(round(9*timespan/14):round(6*timespan/7))=10;
Fd(round(11*timespan/14):round(7*timespan/7))=12;

%%%%%---Tworzenie macierzy M oraz Mp do DMC---%%%%%%%
Mp=zeros(N,D-1);        %macierz ma wymiary Nx(D-1)
MZp=zeros(N,D-1);        %macierz ma wymiary Nx(D-1)
for i=1:D-1 %wypelnianie macierzy Mp
    Mp(1:N, i)=final_odp_skok(i+1:N+i)-final_odp_skok(i);
    MZp(1:N, i)=odp_skok_fd(i+1:N+i)-odp_skok_fd(i);
end
%macierz wspó³czynników odpowiedzi skokowej wymiary(NxNu)
M=zeros(N, Nu);
i=0;
for j=1:Nu  %wypelnianie macierzy trojkatnej dolnej M
    M(j:N,j)=final_odp_skok(1:N-i).';
    i=i+1;
end
I=eye(Nu);              %tworzenie macierzy jednostkowej o wymiarach NuxNu
K=inv(M.'*M+lambda*I)*M.';   %macierz K
KZ=K*MZp;

%%%%%%%%% Algorytm DMC %%%%%%%%%
for t=tau+2:timespan-N %symulacja obiektu i regulatora
    
    
    [V1(t), V2(t)] = object(t-1,h,V1,V2,F1ster(t-1-tau),Fd(t-1),alfa1,alfa2,C1,C2);
    Y(t) = nthroot(V2(t)/C2, 3);
    
    y(t)=Y(t)-Ypp;
    e(t)=yzad(t)-y(t);
    deltaUP(2:D-1)=deltaUP(1:D-2);
    deltaUP(1) = u(t-1)-u(t-2);
    
    deltaZP(2:D-1)=deltaZP(1:D-2);
    deltaZP(1) =Fd(t-1)-Fd(t-2);  
    
    Y0=Mp*deltaUP+y(t);
    Yzad=yzad(t)*ones(N,1);
    deltaU=K*(Yzad-Y0);
    delta_u=deltaU(1) ;%- KZ(1,1:D-1)*deltaZP;
    
    %ograniczenie ró¿nic sygna³u steruj¹cego
    if delta_u>delta_u_max
        delta_u=delta_u_max;
    elseif delta_u<(-delta_u_max)
        delta_u=-delta_u_max;
    end
    u(t)=u(t-1)+delta_u ;
    
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

figure;title('Regulator DMC - sterowanie');hold on
ylabel('F_1_i_n [cm^3/s]');
xlabel('t [s]');
xlim([1 timespan-N])
plot(F1ster(1:timespan-N));hold off

str=sprintf('E=%f', wskaznik_jakosci);
figure;title({'Regulator DMC - wyjœcie obiektu',str});hold on;
ylabel('h_2 [cm]');
xlabel('t [s]');
xlim([1 timespan-N])
plot(Y(1:timespan-N)-Ypp);plot(yzad(1:timespan-N));legend('Y','yzad');hold off;

save('dmc_lin_y.mat','Y','yzad','F1ster')