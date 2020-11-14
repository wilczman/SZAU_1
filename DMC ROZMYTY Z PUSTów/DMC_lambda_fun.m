function [wskaznik_jakosci]=DMC_lambda_fun(parameters, liczba_regulatorow)
lambda=parameters;
    
centra = linspace(-1,1,liczba_regulatorow);
% centra = [1];
params_DMC{1} = [ 29  10  5];
params_DMC{2} = [ 20 26 5;  85 12 6];
params_DMC{3} = [ 20 26 5; 29  10  5; 85 12 6 ];
params_DMC{4} = [ 20 26 5; 40 22 5 ; 60 8 71; 85 12 6 ] ;
params_DMC{5} = [ 20 26 5; 40 22 5 ; 29  10  5; 60 8 71; 85 12 6 ];

%%%%%%%%Punkt Pracy%%%%%%%%
Upp=0;
Ypp=0;

%%%%%%%%%Ograniczenia%%%%%%%
u_min=-1-Upp;
u_max=1-Upp;


%%%%%%%%%-----Parametry DMC------%%%%%%%
D=[50,50,50,50,50]; %horyzont dynamiki (D)
N=[50,50,50,50,50];%horyzont predykcji (N)
Nu=[30,30,30,30,30]; %horyzont sterowania (Nu)(ilosc przyszlych przyrostow wartosci sterowania)

%load('optymalne_parametry_DMC.mat');
%parametry dobrane przez funkcjê fmincon()
%D=150; N=nastawy_DMC_fmincon(1); Nu=nastawy_DMC_fmincon(2); lambda=nastawy_DMC_fmincon(3);
%D=150; N=35; Nu=5; lambda=1 %NAJLEPSZE EKSPERYMENTALNE PARAMETRY SWIATA
%D=200; N=50; Nu=10; % robocze parametry


%%%%%%deklaracja wektorów sygna³ów oraz b³êdów%%%%%%
for nr=1:liczba_regulatorow
    
    D(nr)=params_DMC{liczba_regulatorow}(nr,1); 
    N(nr)=params_DMC{liczba_regulatorow}(nr,2);
    Nu(nr)=params_DMC{liczba_regulatorow}(nr,3);  %parametry dobrane funkcj¹ ga
    
end
kk=450;


y=zeros(1, kk);
Y=zeros(1, kk); 
Y(1:12)=Ypp;
U=zeros(1, kk);  
U(:)=Upp;
e=zeros(1, kk);
u_final=zeros(1,kk);

yzad=zeros(1, kk);
yzad(round(1*kk/7):round(2*kk/7))=4;
yzad(round(2*kk/7):round(3*kk/7))=0;
yzad(round(3*kk/7):round(4*kk/7))=0.2;
yzad(round(4*kk/7):round(5*kk/7))=1.1;
yzad(round(5*kk/7):round(7*kk/7))=3;
yzad=yzad-Ypp;

%%%%%%%   tycie odpowiedzi skokowe  %%%%%%% 
for nr=1:liczba_regulatorow
    u{nr}=zeros(1,kk);
    %odp_skok{nr} = gotowa_odp_skokowa;
    odp_skok{nr} = fun_odp_skok(centra(nr),0.1);
end

%%%%%%% Macierze M, K, Mp %%%%%%%
for nr=1:liczba_regulatorow
    Mp{nr}=zeros(N(nr),D(nr)-1);        %macierz ma wymiary Nx(D-1)
    for i=1:D(nr)-1 %wypelnianie macierzy Mp
       Mp{nr}(1:N(nr), i)=odp_skok{nr}(i+1:N(nr)+i)-odp_skok{nr}(i);
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

    deltaUP{nr}(1:D(nr)-1,1)=0;
    deltaU{nr}=0;
end

%%%%%%%%% Algorytm DMC %%%%%%%%%
for k=12:kk-max(N) 
    %symulacja obiektu i regulatora
    Y(k)=symulacja_obiektu5y(U(k-5),U(k-6),Y(k-1),Y(k-2));
    y(k)=Y(k)-Ypp;
    e(k)=yzad(k)-y(k);
    for nr=1:liczba_regulatorow
        deltaUP{nr}(2:D(nr)-1)=deltaUP{nr}(1:D(nr)-2);
        deltaUP{nr}(1) = u_final(k-1)-u_final(k-2);
        Y0{nr}=Mp{nr}*deltaUP{nr}+y(k);
        Yzad=yzad(k+1:k+N(nr));
        deltaU{nr}=K{nr}*(Yzad-Y0{nr});	
        delta_u{nr}=deltaU{nr}(1);

        u{nr}(k)=u_final(k-1)+delta_u{nr};
    end
    
    %%%%%%%------ROZMYCIE------%%%%%%%%
    for nr = 1:liczba_regulatorow
        w(nr) = gaussmf(u_final(k-1), [gausy(liczba_regulatorow) centra(nr)]);
    end
    
    %œrednia wa¿ona z f.przynale¿noœci
    for nr=1:liczba_regulatorow 
        u_final(k) = u_final(k) + u{nr}(k) * w(nr); 
    end
    u_final(k) = u_final(k) / sum(w);
    
    %ograniczenie sygna³u steruj¹cego
    if u_final(k)>u_max
        u_final(k)=u_max;
    elseif u_final(k)<u_min
        u_final(k)=u_min;
    end
     U(k)=u_final(k)+Upp;
     
end

wskaznik_jakosci=sum(e.^2);
yzad=yzad(1:kk-N(nr))+Ypp;
Y=Y(1:kk-N(nr));
U=U(1:kk-N(nr));
save('temp_lambda_DMC.mat','U','Y','yzad')
end
    