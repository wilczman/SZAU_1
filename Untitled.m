turn = 0;
final_odp_skok(1:500)=0;
for i = 1:5000
    if mod(i,10) == 1
        turn = turn+1;
        final_odp_skok(turn)=odp_skok_scaled(i);
    end
end
plot(final_odp_skok)