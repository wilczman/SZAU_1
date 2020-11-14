function [wskaznik_jakosci]=DMC_fun(parameters)
    lambda=1;
    
%     
%     if parameters(1)<parameters(2)
%        parameters(2)=parameters(1);
%     end
%     
%     if parameters(2)<parameters(3)
%        parameters(3)=parameters(2);
%     end
        
    D=parameters(1);
    N=parameters(2);
    Nu=parameters(3);
    load('skoki.mat', 'start_u', 'start_y', 'end_y', 'skok_odp')
    %%%WKLEJONE
    %%%%%%%    inicjalizacja    %%%%%%
    kk=400;
    

    %%%%%%%%Punkt Pracy%%%%%%%%
    Upp=0;
    Ypp=0;

    %%%%%%%%%Ograniczenia%%%%%%%
    u_min=-1-Upp;
    u_max=1-Upp;
    step = 40;
    
    U=zeros(1, kk);
    u=zeros(1, kk);
    u(:)=start_u;
    U(:)=Upp+start_u;
    Y=zeros(1, kk) + start_y;
    y=zeros(1, kk) + start_y;
%     Y(1:7)=y(1:7)+Ypp;

    e=zeros(1, kk);
    yzad=zeros(1, kk);
    yzad(1:kk)=start_y;
    yzad(step:kk)=end_y;
%     yzad(round(2*kk/5):kk)=end_y;
    yzad=yzad-Ypp;
    
    gotowa_odp_skokowa=fun_odp_skok(start_u, skok_odp);
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
%      Yzad=zeros(N,1);
%      Y0=zeros(N,1);
    %%%%%%%%% Algorytm DMC %%%%%%%%%
    for k=12:kk-N %symulacja obiektu i regulatora
        
        Y(k)=symulacja_obiektu5y(U(k-5),U(k-6),Y(k-1),Y(k-2));
        y(k)=Y(k)-Ypp;
        e(k)=yzad(k)-y(k);
        deltaUP(2:D-1)=deltaUP(1:D-2);
        deltaUP(1) = u(k-1)-u(k-2);  
        Y0=Mp*deltaUP+y(k);
        Yzad=yzad(k+1:k+N)';
        deltaU=K*(Yzad-Y0);	
        delta_u=deltaU(1);

        u(k)=u(k-1)+delta_u;
        %ograniczenie sygna³u steruj¹cego
        if u(k)>u_max
            u(k)=u_max;fprintf('%d %d %d \n',k,u(k),u_max)
        elseif u(k)<u_min
            u(k)=u_min;fprintf('%d %d %d \n',k,u(k),u_min)
        end
         U(k)=u(k)+Upp;

    end
    wskaznik_jakosci=sum(e.^2);
    yzad=yzad(1:kk-N)+Ypp;
    Y=Y(1:kk-N);
    U=U(1:kk-N);
    save('temp_DMC.mat','U','Y','yzad')
end
    