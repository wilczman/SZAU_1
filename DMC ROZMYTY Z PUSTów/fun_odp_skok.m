function odp_skok = fun_odp_skok(center,wymuszenie)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Dlugosc symulacji
kk=400;
%Ograniczenia
u_min=-1;
u_max=1;

U=zeros(1, kk);
U(1:99)=U(1:99)+center; %pozycja startowa
Y=zeros(1, kk);
if center+wymuszenie>1
    wymuszenie=-wymuszenie;
end
U(100:kk)=center+wymuszenie; %skok
for k=1:kk
   if k>7
        Y(k)=symulacja_obiektu5y(U(k-5),U(k-6),Y(k-1),Y(k-2));  
   end
end
odp_skok=(Y(101:kk)-Y(100))/(wymuszenie);

% figure
% hold on
% stairs(odp_skok);
% 
% hold off
% title('Odpowiedz skokowa');
% xlabel('k');
% ylabel('wartoœæ sygna³u');
% 
% figure
% hold on
% stairs(Y);
% stairs(U);
% title('Wyjœcie y(k) i u(k)');
% xlabel('k');
% ylabel('wartoœæ sygna³u');
% end

