%% 
clear all;
clc;
% close all;

%% LED Bandwidth - 20-30 MHz

%Global variables
N = 64;                      % Number of OFDM subcarriers 
Ns = 64;                      % Number of data symbols
roll_off = 1;              % Window roll off
M = 16;                       % Constellation order - M-QAM
m = log2(M);                 % Number of bits per symbol
Ncp = 0;                     % Cyclic prefix
span = 64;

Bsignal = 1e6;               % Signal bandwidth (in Hz)
fc = 2e6;
deltaf=Bsignal/N;
fs = 16e6;
NumSamplesPerSymbol = ceil(fs/Bsignal);



% Bit Rate
%Modulação
hMod = comm.RectangularQAMModulator(M,'BitInput',true);
hDem = comm.RectangularQAMDemodulator(M,'BitOutput',true); 

MaxBitsTransmitted = 10^5;

%% Transmissor

% Transmitted bits
X_bit_OFDM = randi([0 1],N,(Ns*m));

% Modulation
X_data_OFDM = reshape(hMod(X_bit_OFDM(:)),N,Ns);

% Normalization Type
% % Normalization by the longest symbol
E = sqrt(2/3*(M-1));  
%Normalization by the Eavg
% E=sqrt(2/3*(M-1));
% 
X_OFDM = X_data_OFDM;
Group = 8;
% X_OFDM(1:8:end,:) = zeros(8,Ns);
% Hermitian Symmetry Criteria
% X_OFDM_hmc = X_OFDM;
% X_OFDM_hmc(1,:) = 0;
% X_OFDM_hmc(N/2+1,:)= 0;

% for i=1:N/2
%     X_OFDM_hmc(N-i+1,:) = conj(X_OFDM_hmc(i+1,:));
% end
X_OFDM_hmc = [zeros(1,Ns) ;X_OFDM ; zeros(1,Ns); conj(flipud(X_OFDM))];


% IFFT
x_OFDM_hmc = ifft(X_OFDM_hmc);

x_CE = x_OFDM_hmc(:);

rrcFilter = rcosdesign(roll_off,span,NumSamplesPerSymbol,'sqrt');
% 
x_up = upsample(x_CE,NumSamplesPerSymbol);

x_OFDM_filter = conv(x_up,rrcFilter);

x_CE_filter = x_OFDM_filter(span/2*NumSamplesPerSymbol+1:end-span/2*NumSamplesPerSymbol);

CN = sqrt(2/((N*2+2)*NumSamplesPerSymbol*rms(x_CE_filter)^2));
% IH = [1 1 1 1 1 1 1.3 1.3 1.3 1.3];
% VPW = [0.5 0.6 0.7 0.8 0.9 1 0.7 0.83 0.95 1];
% IH = [ 1 1 1 1 1.3 1.3 1.3];
% VPW = [0.7 0.8 0.9 1 0.7 0.83 0.95];
IH = [0.5 0.7 0.8 0.9 1 1.1 1.2 1.3];
VPW = [1 1 1 1 1 1 1 1];

% IH = [1 1 1 1];
% VPW = [0.96 0.97 0.98 1];

% IH = 1.3;
% VPW = 0.95;
IL = 0.1;
n_mat = [];
for ii=1 : length(IH) 
 
    for z=1:2
        if(z==1)
            K=1;  
        else
            K=2; % only change this K
        end
    
        TimeWindow=(2*N+2)*(Ns)/Bsignal;
        dt = 1/fs;                   % seconds per sample
        t= (0:dt:TimeWindow-dt).';     % seconds
        
        % 
%         x_test = cos(2*pi*fc*t);
%         Ps_cos = var(x_test);
    
    
        
        x = cos(2*pi*fc*t + x_CE_filter*K*CN);

        x_test2 = cos(2*pi*fc*t).*cos(x_CE_filter*K*CN);
        x_test3 = sin(2*pi*fc*t).*sin(x_CE_filter*K*CN);
        x_test = x_test2 - x_test3;

%         x_test2 = cos(2*pi*fc*t);
%         x = exp(1j*x_CE_filter*K*CN);
        
        
        % Add CP
        Ncp_index = length(x)*Ncp;
        x_cp = [x(end+1-Ncp_index:end); x];
    
    %%  DIMMING CONTROL    
    %INTENSITY-DOMAIN
        s_dimming = reshape(x_cp,[],Ns);
        [l,c] = size(s_dimming);
        index_dimming = VPW(ii)*l;
        
        s_dimming(index_dimming+1:end,:) = s_dimming(index_dimming+1:end,:)*IL;
        s_dimming(1:index_dimming,:) = s_dimming(1:index_dimming,:)*IH(ii);
        
        s_dimming = s_dimming(:);
        n_new = rms(s_dimming).^2;

    %% LED 
    %100 mA
    % p = [-0.0154 0.2957 0.7057 -0.1459];
    
    %200 mA
%     p= [-0.0617 0.2067 1.06 -0.1045];
    
    %250 mA
    % p = [-0.0437 0.0218 0.7955 -0.0238];
    
%     current = 0.2 + s_dimming;
%     sigmax = rms(s_dimming+0.2);
%     s2 = polyval(p,current);

    %%
        %     EbN0_db=0:20;
            Ber_vector = [];
            %channel
            if(z==1)
                EbN0_db = 20;
            else
                n_mat = [n_mat n_new];
                EbN0_db=6:24;
%                 EbN0_db = 30;
                SNR=EbN0_db+10*log10(log2(M)/(1+roll_off))-10*log10(NumSamplesPerSymbol);
            end
            for i=1:length(EbN0_db)
                TotalError = 0;
                TotalBits = 0;
                while (TotalBits<MaxBitsTransmitted)
            %          y = awgn(x_cp,SNR(i),'measured');
                    if(z~=1)
                         y = ElectricChannel_CE_OFDM(s_dimming,EbN0_db(i),M,NumSamplesPerSymbol,Ps_CE);
%                          TotalBitRs = 10e6;
                    else
                         y = s_dimming;
                         TotalBits = 10e6;
                    end
%                      y = s_dimming;

                     %Receiver
                     % COMPENSAR O CLIPPING FEITO, um x o outro tem de dar
                    y_rec = reshape(y,[],Ns);

  

                    y_rec(1:index_dimming,:) = y_rec(1:index_dimming,:)/IH(ii);
                    y_rec(index_dimming+1:end,:) = y_rec(index_dimming+1:end,:)/IL;    


                     %remove CP
                     r = y_rec(Ncp_index+1:end);
            
%                      r2 = exp(-i* r);
                      r2 = cos(-2*pi*fc*t + r.');
                     
                     y_OFDM_filter = conv(r2,rrcFilter);
                     SignalOut = y_OFDM_filter(span*NumSamplesPerSymbol/2+1:end-span*NumSamplesPerSymbol/2);
                     Signal_down = downsample(SignalOut,NumSamplesPerSymbol); 
                     
        
                     if(z==1)
                         Ps_CE = var(SignalOut);
                     end
        
                     test = (acos(Signal_down));
%                      real_part = real(Signal_down);
%                      imag_part = imag(Signal_down);
        
%                      test = atan(imag_part./real_part);
            
                     r_par = reshape(test,[],Ns);
            
                     Test = fft(r_par);
                    Test2 = Test(2:N+1,2:end-1);
                    Eavi=sqrt(sum((abs(Test2(:)).^2))/(N*(Ns-2))); 
                    I_Eq=Test2./Eavi;
            
                    Kdesn=sqrt(2/3*(M-1)); 
                    SymRec=Kdesn*I_Eq;
%                     if(z==2)
%                         scatterplot(SymRec(:))
%                     end
            %%
                    BitsRec = hDem(SymRec(:)); 
                        
                    SymTransmitted = X_OFDM(:,2:end-1);
            
                    test =[SymRec(:) SymTransmitted(:)];
                    BitsTransmitted = hDem(SymTransmitted(:));
                
                    test = xor(BitsTransmitted,BitsRec);
                    Error = length(find(test == 1));
            
                    TotalError =TotalError+Error;
                    TotalBits = TotalBits + length(BitsRec);
                end
                
                Ber = TotalError/TotalBits;
                Ber_vector = [Ber_vector; Ber];
            end
            
    end
    semilogy(EbN0_db,Ber_vector);
    hold on
    ylim([10^-4 1])
    xlabel('EbN0')
    ylabel('Ber')
    title('CE-OFDM')
end
legendStrings = "n_{normalized} = " + string(round(n_mat,2));
legend(legendStrings)
grid