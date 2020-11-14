function wskaznik_jakosci=zad3_PID_fun(parameters)
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

    kk=400; %d³ugoœæ symulacji

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
    yzad=yzad-Ypp;


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
    wskaznik_jakosci=sum(e.^2);
end