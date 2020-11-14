clear all

%%%%%%%    inicjalizacja    %%%%%
load('gotowa_odp_skokowa.mat');
transpose(gotowa_odp_skokowa);

%%%%%%%%Punkt Pracy%%%%%%%%
Upp=0;
Ypp=0;

%%%%%%%%%Ograniczenia%%%%%%%
u_min=-1-Upp;
u_max=1-Upp;

%%%%%%%%%-----Parametry DMC------%%%%%%%
lambda=1; %parametr lambda np. 1
D=37; %horyzont dynamiki (D)
N=30;%horyzont predykcji (N)
Nu=1; %horyzont sterowania (Nu)(ilosc przyszlych przyrostow wartosci sterowania)

kk=450;
% load('optymalne_parametry_DMC.mat');
%parametry dobrane przez funkcjê fmincon()
% D=150; N=nastawy_DMC_fmincon(1); Nu=nastawy_DMC_fmincon(2); lambda=nastawy_DMC_fmincon(3);
%D=150; N=35; Nu=5; lambda=1 %NAJLEPSZE EKSPERYMENTALNE PARAMETRY SWIATA
%D=200; N=50; Nu=10; % robocze parametry

%%%%%%deklaracja wektorów sygna³ów oraz b³êdów%%%%%%
U=zeros(1, kk);
u=zeros(1,kk);
U(:)=Upp;
e=zeros(1, kk);
Y=zeros(1, kk);
y=zeros(1, kk);
Y(1:12)=Ypp;
yzad=zeros(1, kk);

yzad(round(1*kk/7):round(2*kk/7))=4;
yzad(round(2*kk/7):round(3*kk/7))=0;
yzad(round(3*kk/7):round(4*kk/7))=0.2;
yzad(round(4*kk/7):round(5*kk/7))=1.1;
yzad(round(5*kk/7):round(7*kk/7))=3;

yzad=yzad-Ypp;

Mp=zeros(N,D-1);        %macierz ma wymiary Nx(D-1)
for i=1:D-1 %wypelnianie macierzy Mp
   Mp(1:N, i)=gotowa_odp_skokowa(i+1:N+i)-gotowa_odp_skokowa(i);
end
%macierz wspó³czynników odpowiedzi skokowej wymiary(NxNu)
M=zeros(N, Nu);  
i=0;
for j=1:Nu  %wypelnianie macierzy trojkatnej dolnej M
   M(j:N,j)=gotowa_odp_skokowa(1:N-i).';  
   i=i+1;
end

I=eye(Nu);              %tworzenie macierzy jednostkowej o wymiarach NuxNu
K=inv(M.'*M+lambda*I)*M.';   %macierz K

deltaUP(1:D-1,1)=0;
deltaU=0;

%%%%%%%%% Algorytm DMC %%%%%%%%%
for k=12:kk-N %symulacja obiektu i regulatora
    Y(k)=symulacja_obiektu5y(U(k-5),U(k-6),Y(k-1),Y(k-2));
    y(k)=Y(k)-Ypp;
    e(k)=yzad(k)-y(k);
    deltaUP(2:D-1)=deltaUP(1:D-2);
    deltaUP(1) = u(k-1)-u(k-2);  
    Y0=Mp*deltaUP+y(k);
    Yzad=yzad(k+1:k+N);
    deltaU=K*(Yzad-Y0);	
    delta_u=deltaU(1);

    u(k)=u(k-1)+delta_u;
    %ograniczenie sygna³u steruj¹cego
    if u(k)>u_max
        u(k)=u_max;
    elseif u(k)<u_min
        u(k)=u_min;
    end
     U(k)=u(k)+Upp;

end
wskaznik_jakosci=sum(e.^2);
yzad=yzad(1:kk-N)+Ypp;
Y=Y(1:kk-N);
U=U(1:kk-N);

%%%%%%%%prezentacja wyników symulacji%%%%%%%%
figure; stairs(U);hold on
title('Sygna³ sterowania DMC'); xlabel('k');
figure; stairs(Y); hold on; 
stairs(yzad,':'); xlim([0, kk-N]);  % ylim([-1, 7]); 
title('Wyjœcie regulatora DMC'); xlabel('k'); ylabel('wartoœæ sygna³u');
legend('wyjœcie y(k)','wartoœæ zadana');


% %%%ZAPIS DO PLIKU%%%%%%%
% % odkomentuj by w³¹czyæ
j=0; k=0;
for k=1:kk       
       j=j+1;
       table_U(1,k)=double(j);
       table_Y(1,k)=double(j);
       table_yzad(1,k)=double(j);
end
%macierze wykorzystane do zapisu
table_U(2,1:420)=U;
table_Y(2,1:420)=Y;
table_yzad(2,1:420)=yzad;

fname_U = sprintf('outU.txt');
fname_Y = sprintf('outY.txt');
fname_yzad = sprintf('outYzad.txt');
mkdir_status_Y=mkdir(sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_DMC\\'));
if mkdir_status_Y
    savdir_Y =  sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_DMC\\');
    
   fileID = fopen([savdir_Y fname_U],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_U);
   fclose(fileID);
   fileID = fopen([savdir_Y fname_Y],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_Y);
   fclose(fileID);
   fileID = fopen([savdir_Y fname_yzad],'w');
   fprintf(fileID,'%6.3f %6.3f\r\n',table_yzad);
   fclose(fileID);
else 
   disp('Nie udalo sie stworzyc folderów')
end
warning('on','all') %wlaczenie warningow