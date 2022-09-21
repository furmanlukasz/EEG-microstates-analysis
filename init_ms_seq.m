function [ms_seq] = init_ms_seq(ms_at_peaks,locs)
% INIT_MS_SEQ estende il vettore ms_at_peaks contenente, per ogni picco
% della GFP, l'indice del microstato attivo in quell'istante, in modo da
% far si che il vettore descriva la sequenza dei microstati per ogni
% istante del segnale eeg originale (e non soltanto per i picchi).
    ms_sequence = zeros(9632,1);
    ms_sequence(locs,1) = ms_at_peaks;
    if ms_sequence(1,1) == 0
       i = 2;
       while ms_sequence(i,1) == 0
           i=i+1;
       end
       for j = 1 : i 
          ms_sequence(j,1) = ms_sequence(i,1); 
       end
    end
    prev = ms_sequence(j,1);
    last = j;
    for i = j : 9632
        if ms_sequence(i,1)~=0 && ms_sequence(i,1)~=prev
            distance = i - last;
            border = int16(distance/2);
            for t = last : i%+distance
                if t<border && t <= 9632
                    ms_sequence(t,1) = prev;
                end
                if t>=border && t <= 9632
                    ms_sequence(t,1) = ms_sequence(i,1);
                end
            end
            last = i;
            prev = ms_sequence(i,1);
        elseif ms_sequence(i,1)==0
            ms_sequence(i,1) = prev;
        end
    end
    ms_seq = ms_sequence;
end

