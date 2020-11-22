function sigma = gausy(liczba_regulatorow, pocz, koniec)
if liczba_regulatorow >1 
    x=linspace(pocz,koniec,1000);
    centra = linspace(pocz,koniec,liczba_regulatorow);

    d_skok=abs(pocz-koniec)/(liczba_regulatorow-1);    %odleglosc miedzy pkt przeciêcia
    d0=pocz+d_skok/2;  %pierwszy z punktow przeciecia 

    x_przeciecia=d0;     % punkt przeciecia
    sigma=zeros(liczba_regulatorow,2);
    for i=1:liczba_regulatorow

       sigma(i,:)=roots([2 0  ((x_przeciecia-centra(i)).^2)/(log(0.5))]);
       x_przeciecia=x_przeciecia+d_skok;
       gs{i}=gaussmf(x, [sigma(i,1),centra(i)]);
    end
    sigma=sigma(1,1);
    
    
else
    sigma=1;
end
