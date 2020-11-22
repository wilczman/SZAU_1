clear all;

liczba_regulatorow = 5;
%zakres wysokosci
h2_pocz=15;
h2_koniec=45;
x = (h2_pocz:0.1:h2_koniec)';
kk=101;
center=linspace(h2_pocz,h2_koniec,liczba_regulatorow);

figure;
hold on;
for nr=1:liczba_regulatorow
    y{nr}=gaussmf(x, [gausy(liczba_regulatorow, h2_pocz, h2_koniec) center(nr)]);
    plot(x, y{nr});
end
