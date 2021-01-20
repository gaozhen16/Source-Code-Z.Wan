clear all;
clc;
SNR_dB = 10;  SNR_linear=10.^(SNR_dB/10.);
N_iter = 1000;
N1 = 16;
N2 = 8;
maxNum = 16;32;
Q = 64;
L = 3;
IndexN1 = 0:1:N1-1;
IndexN2 = 0:1:N2-1;
NMSE1 = zeros(N_iter,length(SNR_dB));
NMSE2 = zeros(N_iter,length(SNR_dB));
NMSE3 = zeros(N_iter,length(SNR_dB));
cout2 = zeros(N_iter,length(SNR_dB));
cout3 = zeros(N_iter,length(SNR_dB));

U1 = zeros(N1,N1);
for i = 1:N1
    for j = 1:N1
        U1(i,j) = sqrt(1/N1)*exp(1i*2*pi/N1*(i-1-(N1-1)/2)*(j-1-(N1-1)/2));
    end
end

U2 = zeros(N2,N2);
for i = 1:N2
    for j = 1:N2
        U2(i,j) = sqrt(1/N2)*exp(1i*2*pi/N2*(i-1-(N2-1)/2)*(j-1-(N2-1)/2));
    end
end
t1 = clock;
for i_snr=1:length(SNR_dB) 
    i_snr
    sigma2=1/SNR_linear(i_snr);
    for iter = 1:N_iter
        tic
        AoD1 = -pi/2 + rand * pi;(2*rand())/N1*pi;
        AoD2 = -pi/2 + rand * pi;(2*rand())/N2*pi;
        Abh1 = sqrt(1/N1)*exp(1j*IndexN1*AoD1)';
        Abh2 = sqrt(1/N2)*exp(1j*IndexN2*AoD2)';
        Abhh = kron(U1 * sqrt(N1)*Abh1,U2 * sqrt(N2)*Abh2);
        for i = 2:L
            AoD1 = (2*rand())/N1*pi;
            AoD2 = (2*rand())/N2*pi;
            Abh1 = sqrt(1/N1)*exp(1j*IndexN1*AoD1)';
            Abh2 = sqrt(1/N2)*exp(1j*IndexN2*AoD2)';
            Abh = sqrt(N1*N2)*kron(Abh1,Abh2);
            Abhh = Abhh + 0.01*kron(U1 * sqrt(N1)*Abh1,U2 * sqrt(N2)*Abh2);
        end
        
        W = (1/sqrt(Q))*(2*(rand(Q,N1*N2)<0.5) - 1);
%         W = OVSF(N1*N2);
%         W = W(:,1:Q);
%         W = 2 * W';
        
        noise = sqrt(sigma2)*(randn(Q,1)+1i*randn(Q,1))/sqrt(2);
        yy = W*Abhh + noise;
 
%         h_hat1 = OMP_new(yy, W, maxNum*L, maxNum*L);
%         NMSE1(iter,i_snr) = norm(Abhh - h_hat1,2)^2/norm(Abhh,2)^2;
%         [h_hat2 ,cout2(iter,i_snr)] = SSD(yy,W,L,N1,N2);
%         NMSE2(iter,i_snr) = norm(Abhh - h_hat2,2)^2/norm(Abhh,2)^2;
        [h_hat3 ,cout3(iter,i_snr)] = mySSD2_2(yy,W,N1,N2,L,maxNum);
        NMSE3(iter,i_snr) = norm(Abhh - h_hat3,2)^2/norm(Abhh,2)^2;
        
        toc
        disp(['  DC, overhead = ' num2str(Q) ', SNR = ' num2str(SNR_dB(i_snr)) ', iter_max = ' num2str(N_iter) ', iter_now = ' num2str(iter) ...
            ', NMSE_now = ' num2str(10*log10(NMSE3(iter,i_snr))) 'dB, NMSE_avg = ' num2str(10*log10(mean(NMSE3(1:iter,i_snr),1))) 'dB.']);
    end
end
t2 = clock;
Time = etime(t2,t1)/N_iter;
NMSE1avg = mean(NMSE1,1);
NMSE2avg = mean(NMSE2,1);
NMSE3avg = mean(NMSE3,1);
cout2avg = mean(cout2,1);
cout3avg = mean(cout3,1);

figure();
semilogy(SNR_dB,NMSE1avg,'-ob','Linewidth',1.5);
hold on;
semilogy(SNR_dB,NMSE2avg,'-*r','Linewidth',1.5);
semilogy(SNR_dB,NMSE3avg,'-+k','Linewidth',1.5);
grid on;
legend('OMP-based channel estimation','ASD-based channnel estimation','DC-based channel estimation');
xlabel('SNR (dB)');
ylabel('NMSE (dB)');
title(['N1 \times N2 = ', num2str(N1), '\times', num2str(N2), ', L = ', num2str(L), ', CR = ', num2str(Q/N1/N2)])