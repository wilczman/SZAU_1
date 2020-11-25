clear all;

liczba_regulatorow = 2;
%zakres wysokosci
figure;
hold on;
i=1
for liczba_regulatorow=2:5
    h2_pocz=15;
    h2_koniec=45;
    x = (h2_pocz:0.1:h2_koniec)';
    kk=101;
    center=linspace(h2_pocz,h2_koniec,liczba_regulatorow);

    
    for nr=1:liczba_regulatorow
        y{nr}=gaussmf(x, [gausy(liczba_regulatorow, h2_pocz, h2_koniec) center(nr)]);
        
        subplot(4,1,i);
        str=sprintf('Liczba modeli = %d', liczba_regulatorow);
        title(str);
        xlabel('h_2')
        hold on
        plot(x, y{nr});
        
    end
    
    i=i+1;
end