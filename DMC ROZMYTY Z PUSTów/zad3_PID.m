clear all;
%inicjalizacja
%Punkt Pracy
Upp=0;
Ypp=0;

%Ograniczenia
u_min=-1;
u_max=1;
delta_umax=2;
u_max=u_max-Upp;
u_min=u_min-Upp;

<<<<<<< HEAD
kk=1450; %d³ugoœæ symulacji
=======
kk=400; %d³ugoœæ symulacji
>>>>>>> b50d5f5d2b0370571bfdf8d1fbd682ad81a06cbf

%deklaracja wektorów sygna³ów oraz b³êdów
U=zeros(1, kk);
u=zeros(1, kk);
U(:)=Upp;
Y=zeros(1, kk);
y=zeros(1, kk);
Y(1:12)=Ypp;
e=zeros(1, kk);
yzad=zeros(1, kk);

yzad(round(1*kk/7):round(2*kk/7))=4;
yzad(round(2*kk/7):round(3*kk/7))=0;
yzad(round(3*kk/7):round(4*kk/7))=0.2;
yzad(round(4*kk/7):round(5*kk/7))=1.1;
yzad(round(5*kk/7):round(7*kk/7))=3;
<<<<<<< HEAD

=======
>>>>>>> b50d5f5d2b0370571bfdf8d1fbd682ad81a06cbf
yzad=yzad-Ypp;

load('optymalne_parametry_PID_zad3.mat');
%K=K_kryt*0.5; Ti=22; Td=7; %NAJLEPSZE EKSPERYMENTALNE PARAMETRY

K=nastawy_PID_fmincon(1); Ti=nastawy_PID_fmincon(2); Td=nastawy_PID_fmincon(3);  %parametry dobrane funkcj¹ fmincon

T=0.5;

r2=K*Td/T;
r1=K*(T/(2*Ti)-2*Td/T-1);
r0=K*(1+T/(2*Ti)+Td/T);

for k=12:kk %g³ówna pêtla symulacyjna
    %symulacja obiektu    
    Y(k)=symulacja_obiektu5y(U(k-5),U(k-6),Y(k-1),Y(k-2));
    y(k)=Y(k)-Ypp;
    %uchyb regulacji
    e(k)=yzad(k)-y(k);
    %sygna³ steruj¹cy regulatora PID
    %u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    delta_u=r2*e(k-2)+r1*e(k-1)+r0*e(k);
   
     u(k)=delta_u+u(k-1);
    if u(k)>u_max
        u(k)=u_max;
    elseif u(k)<u_min
        u(k)=u_min;
    end
     U(k)=u(k)+Upp;
end
yzad=yzad+Ypp;
wskaznik_jakosci=sum(e.^2);
%wyniki symulacji
figure; stairs(U);
title('Sterowanie regulatora PID'); xlabel('k'); ylabel('wartoœæ sygna³u');
figure; stairs(Y); hold on; stairs(yzad,':'); legend('wyjœcie y(k)','wartoœæ zadana');
title('Wyjœcie regulatora PID'); xlabel('k'); ylabel('wartoœæ sygna³u');

% %%%ZAPIS DO PLIKU%%%%%%%
j=0; k=0;
for k=1:kk       
       j=j+1;
       table_U(1,k)=double(j);
       table_Y(1,k)=double(j);
       table_yzad(1,k)=double(j);
       table_Z(1,k)=double(j);
end
%macierze wykorzystane do zapisu
table_U(2,:)=U;
table_Y(2,:)=Y;
table_yzad(2,:)=yzad;

fname_U = sprintf('outU.txt');
fname_Y = sprintf('outY.txt');
fname_yzad = sprintf('outYzad.txt');
mkdir_status_Y=mkdir('strojenie_PID\fmincon\Y');
if mkdir_status_Y
    savdir_Y =  sprintf('C:\\Users\\Kuba\\Desktop\\GIT\\PUST_3\\strojenie_PID\\');
    
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