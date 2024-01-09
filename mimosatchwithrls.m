clc
clear all
close all
%%%%%%%%%%%%in the name of ALLAH%%%%%%%%%

%%%%%%%%%%%%initial values%%%%%%%%%
L = 1200;
ts = 1/10;
fd = 0;
M = 4;
m = log2(M);
SNR_range = 0:1:30;
k = 3;
nRx = 2;
nTx = 2;
L_chan = 1;
p01 = .028;
p10 = .02;
% u = 0.0001;
lplant = 2;
delta = 0.1;
lambda = 1;
numSNR = length(SNR_range);
env_factor = 20;
%%%%%%%%%

data = randi([0,M-1] ,1 ,L);
x = reshape(data ,[nTx ,L/nTx]);
xmod = pskmod(x ,M);

desired = xmod;

berVec = zeros(1 ,numSNR);

for i = 1:nRx
    for j = 1:nTx
        hchan1 (i, j) = rayleighchan(ts ,fd);
        hchan1 (i, j).storePathGains = 1;
        hchan1 (i, j).ResetBeforeFiltering = 0;
        
        hchan2 (i, j) = ricianchan(ts ,fd ,k);
        hchan2 (i, j).storePathGains = 1;
        hchan2 (i, j).ResetBeforeFiltering = 0;
    end;
end;



for SNR_i = 1 :numSNR

    W0 = randn (nRx ,nTx);
    W1 = randn (nRx ,nTx);
    yW = zeros(nRx ,L/nTx);
%%%%%%%%%%%%%
a = 0;
 
for l = 1:(L/nTx)
    true_state_sequence (l) = a;
    if (a == 0)
        for i = 1:nRx
            y (i, l) = 0;
            for j = 1:nTx
                y (i, l) = y (i, l) + filter (hchan1 (i, j), xmod (j,l));
            end;
        end;
        r = rand;
        if (r < p01)
            a = 1;
        end;
    else
        for i = 1:nRx
            y (i, l) = 0;
            for j = 1:nTx
                y (i, l) = y (i, l) + filter (hchan2 (i, j), xmod(j,l));
            end;
        end;
        r = rand;
        if (r < p10)
            a = 0;
        end;
    end;
    
end;

ynoisy = awgn (y ,SNR_range(SNR_i));

   ch_state_sequence = chanel_state_sequence (ynoisy, env_factor);
   
%         figure;
%         hold on;
%         plot (true_state_sequence, 'b');
%         plot (ch_state_sequence, 'r');
    
%%%equalizer%%%

psy_inv0 = diag((ones(nRx ,1)/delta));
psy_inv1 = diag((ones(nRx ,1)/delta));

for kk = 1:(L/nRx)
    if (1-true_state_sequence (kk))
        
        u0 = psy_inv0 * ynoisy(: ,kk);
        k0 = u0/(lambda + ynoisy(: ,kk)'*u0);
        yW(1:nRx, kk)= W0' * ynoisy(:, kk);
        e(:, kk) = desired(1:nRx, kk) - yW(1:nRx, kk);
        W0 = W0 + k0 * e(: ,kk)';
        psy_inv0 = (psy_inv0 - (k0 *ynoisy(: ,kk)' *psy_inv0));
    else
        u1 = psy_inv1 * ynoisy(: ,kk);
        k1 = u1/(lambda + ynoisy(: ,kk)'*u1);
        yW(1:nRx, kk)= W1' * ynoisy(:, kk);
        e(:, kk) = desired(1:nRx, kk) - yW(1:nRx, kk);
        W1 = W1 + k1 * e(: ,kk)';
        psy_inv1 = (psy_inv1 - (k1 *ynoisy(: ,kk)' *psy_inv1));
    end
end
%     figure
%     plot (sum (abs (e)));
    
    ydemod = pskdemod(yW ,M);
    berVec(SNR_i) = sum (sum (xor(x ,ydemod )));

end

BER_sim = berVec / L;


% Compute theoretical performance results, for comparison
BERtheory = berfading(SNR_range,'psk',M,1);

% Plot BER results
semilogy(SNR_range,BERtheory,'b-',SNR_range,BER_sim,'r-');
grid on
legend('Theoretical BER','Empirical BER');
xlabel('SNR (dB)'); ylabel('BER');
title(['MIMO Sattelite Channel (fd = ' num2str(fd) ') ']);
