function wskaznik_jakosci = PID_fun(parameters,start_y,end_y,start_u)
    K=parameters(1);
    Ti=parameters(2);
    Td=parameters(3);

    %inicjalizacja
    %Punkt Pracy
    Upp=0;
    Ypp=0;

    %Ograniczenia
    u_min=-1;
    u_max=1;
    u_max=u_max-Upp;
    u_min=u_min-Upp;
    delta_u_max=2;
    
    kk=400; %d³ugoœæ symulacji

    %deklaracja wektorów sygna³ów oraz b³êdów
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

    T=0.5;

    r2=K*Td/T;
    r1=K*(T/(2*Ti)-2*Td/T-1);
    r0=K*(1+T/(2*Ti)+Td/T);

    for k=step:kk %g³ówna pêtla symulacyjna
        %symulacja obiektu    
        Y(k)=symulacja_obiektu5y(U(k-5),U(k-6),Y(k-1),Y(k-2));
        y(k)=Y(k)-Ypp;
        %uchyb regulacji
        e(k)=yzad(k)-y(k);
        %sygna³ steruj¹cy regulatora PID
        %u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
        delta_u=r2*e(k-2)+r1*e(k-1)+r0*e(k);
        
        %ograniczenie ró¿nic sygna³u steruj¹cego 
        if delta_u>delta_u_max
             delta_u=delta_u_max; 
        elseif delta_u<(-delta_u_max)
             delta_u=-delta_u_max;
        end
        u(k)=delta_u+u(k-1);
        %ograniczenie wielkosci sygnalu
        if u(k)>u_max
            u(k)=u_max;
        elseif u(k)<u_min
            u(k)=u_min;
        end
         
        
        U(k)=u(k)+Upp;
    end
    wskaznik_jakosci=sum(e.^2);
    save('temp_to_del.mat','U','Y','yzad')
end