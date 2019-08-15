%Thesis work: Estimating HRV features with missing data (for normal
%patients)
% Author: Md. Surat-E-Mostafa

%% importes data information

% normal patients data, 360 sampling rate.
% 100.m data has 2 lead , first lead don need to do filter , peak 1100 but
% second lead value need to do filter, peak value 55.

rng(86757658) % make random non negative number

%% MIT-BIH arrythmia database(mitdb), frquency 360, 2 leads data.
load('100m.mat')

Signal_lead1 = val(1,:);

%%if the dat need to do smooth
% % Ecg1_noise=Signal_lead1(1:1000*60*8);
% Ecg1_noise=Signal_lead1(1:360*60*10);
% Ecg1_sm= smooth(Ecg1_noise, 20*4, 'sgolay',3);
% Ecg1= Ecg1_noise'-Ecg1_sm;


% Ecg1=Signal_lead1(1:360*60*10); % 10 mins data
Ecg1=Signal_lead1(360*60*1:360*60*6); % 10 mins data
%   Ecg1=Signal_lead1(1:360*60*1);    % 1 mins data
t1= 0:1/360 :(length(Ecg1)-1)/360 ;
 
figure(501);
subplot(2,1,1);
plot(t1,Ecg1);
title('10 Mins MIT-BIH arrythmia data)');
xlabel('Time in seconds');
ylabel('Amplitude');
hold on
 %% Base line wandering figure part
 
Ecg1_noise=Signal_lead1(360*60*15:360*60*15.3);
t1= 0:1/360 :(length(Ecg1_noise)-1)/360 ;
% subplot(211)
% plot (t1,Ecg1_noise)
% legend('Raw ECG');

Ecg1_sm= smooth(Ecg1_noise, 20*4, 'sgolay',3);
Ecg1= Ecg1_noise'-Ecg1_sm;
xlabel('Time(s)')
ylabel( 'Amplitude')

figure (1)
subplot(211)
plot(t1,Ecg1);

xlabel('Time(s)')
ylabel( 'Amplitude')
title('ECG');
hold on
[qrspeaks,locs] = findpeaks(Ecg1,t1,'MinPeakHeight',125 ,...
    'MinPeakDistance',0.05);
 plot(locs,qrspeaks,'o');
 legend('Baseline wander reduced ECG', 'R-peaks annotation');
 
 

R_interval = diff(locs);% for normal command
t_RR= locs(2:end);

subplot(212)
plot(t_RR , R_interval*1000, '*');
title('RR interval for complete data ');
hold on
t_int= t_RR(1): .12 : t_RR(end);
R_interval_intp = interp1(t_RR,R_interval*1000,t_int,'spline');
plot(t_int,R_interval_intp, 'o');
title('R-R Interval ')
xlabel('Time(s)'); ylabel('Time(ms)');
legend('Tachogram', 'Interpolated RRI signal at 4Hz');



%% R peak calculation
%% RR interval 1100 for 10 mins data,1700 af

[qrspeaks,locs] = findpeaks(Ecg1,t1,'MinPeakHeight',1100 ,...
    'MinPeakDistance',0.05);

 plot(locs,qrspeaks,'o');

R_interval = diff(locs);% for normal command
t_RR= locs(2:end);

subplot(212)
plot(t_RR, R_interval);
title('RR interval for complete data ');





%% ectoppic beat detection and replace it by NaN

dRRI=diff(locs,2);  % nth differnce, 2 difference like [1 2 3], so 2-1= , 3-2=
ectopic_inds = dRRI > mad(dRRI,1)*6;   % y is median(abs(X-median(X)))
ectopic_inds( find(ectopic_inds)+1 ) = 1;
R_interval( ectopic_inds ) = nan;

figure(502);
plot(dRRI);
hold on;
plot( find(ectopic_inds), dRRI(ectopic_inds), 'rd')

figure(503)
subplot(311);
plot(R_interval*1000);

R_interval_full = R_interval;




%%  Interpolation 
%sfg :  generalized moving average with filter coefficients determined by an unweighted linear 
%least-squares regression and a polynomial model of specified degree (default is 2). 
%The method can accept nonuniform predictor data.

t_int= t_RR(1): .25 : t_RR(end);
R_interval_intp = interp1(t_RR( ~ectopic_inds ),R_interval( ~ectopic_inds ),t_int,'spline');

subplot(3,1,2);
plot(t_int,R_interval_intp);
title('Tachogram(Normal Patients) ')
xlabel('Time(sec'); ylabel('RR(s)');

sgf_signal= smooth(R_interval_intp, 20*4, 'sgolay',3);
detrended_RR_interv = R_interval_intp'-sgf_signal


subplot(313);
plot(t_int,detrended_RR_interv,'.-');
title('Smoothing figure');

%LS calculation
[pxx_lsorg F_lsorg]= plomb( R_interval( ~ectopic_inds ), t_RR( ~ectopic_inds ) );
 

%fft calculation
L = 360;  % 360 for 10min daata,1000
noverlap = 180;  %180 for 10min data,500
[signal_org,f_norm_fftorg] = pwelch(detrended_RR_interv,hann(L),noverlap,F_lsorg,4,...
    'onesided','ConfidenceLevel',0.95);
% f_norm_fftorg = ((f_org-min(f_org)) ./ (max(f_org)-min(f_org))) ;
ylim([0 0.1]);


 
figure(504)
subplot(311)
plot(f_norm_fftorg(f_norm_fftorg>0.05),signal_org(f_norm_fftorg>0.05));
plot(f_norm(f_norm>0.05),signal(f_norm>0.05), 'DisplayName', ['Missing data percentage=' num2str(p*100)]);
title('Power spectrum estimation via FFT base methodfor normal patient');
xlabel('frequency(Hz)');
ylabel('PSD-FFT');
legend('Original Data Spectrum');
ax = gca;
ax.YAxis.Exponent = 0;
% ylim([0 0.1]);
% set(gca, 'YScale', 'log')



% ar calculation
[Pxx_ar_org,F_ar_norm_org] = pburg(detrended_RR_interv,11,F_lsorg,4);
% [Pxx_ar_org,F_ar_norm_org] = pburg(detrended_RR_interv,10,length(detrended_RR_interv),360);

subplot(3,1,2);
% F_ar_norm_org = ((F_ar_org-min(F_ar_org)) ./ (max(F_ar_org)-min(F_ar_org))) ;
plot(F_ar_norm_org(F_ar_norm_org>0.05),Pxx_ar_org(F_ar_norm_org>0.05));
title('Power spectrum estimation via AR based method for normal patient');
% legend show
xlabel('frequency(Hz)');
ylabel('PSD-AR');
legend('Original Data Spectrum');
ax = gca;
ax.YAxis.Exponent = 0;
% ylim([0 0.1]);
% set(gca, 'YScale', 'log')


%LS  figure 
subplot(313);
plot(F_lsorg(F_lsorg>0.05), pxx_lsorg(F_lsorg>0.05));
% F_norm = ((F_lsorg-min(F_lsorg)) ./ (max(F_lsorg)-min(F_lsorg))) ;
%  plot(F_norm(F_norm>0.05),pxx_lsorg(F_norm>0.05));


title('Power spectrum estimation vis LS(Lomb-Scargle) based method for normal patient');
% legend show
xlabel('frequency(Hz)');
ylabel('PSD-LS');
legend('Original Data Spectrum')
ax = gca;
ax.YAxis.Exponent = 0;
% ylim([0 0.1]);
% set(gca, 'YScale', 'log')






%% entropy calculation for original data mean reference data

%entropy for normal rythmfft
 P = signal_org./ sum( signal_org );
 P_nan = nan( size(P) );
 P_nan( P>0 ) = P( P>0 );
%  P_nan = P;
%  P_nan( P<0.00001 ) = NaN;
 Entropy_fft_originalN = - nansum( P.*log2(P) );  %nansum(X) is the sum of X, computed after removing NaN values.
 
 %entropy for normal rythm ar
 P = Pxx_ar_org./ sum( Pxx_ar_org );
 P_nan = nan( size(P) );
 P_nan( P>0 ) = P( P>0 );
%  P_nan = P;
%  P_nan( P<0.00001 ) = NaN;
 Entropy_ar_originalN = - nansum( P.*log2(P) );
 
 %entropy for normal rythm ls
 P = pxx_lsorg./ sum( pxx_lsorg );
 P_nan = nan( size(P) );
 P_nan( P>0 ) = P( P>0 );
%  P_nan = P;
%  P_nan( P<0.00001 ) = NaN;
 Entropy_ls_originalN = - nansum( P.*log2(P) );
 
 [Entropy_normaldata_fft_ar_ls]= [Entropy_fft_originalN Entropy_ar_originalN  Entropy_ls_originalN];

dlmwrite('Entropy_normaldata_fft_ar_ls.txt',Entropy_normaldata_fft_ar_ls, 'delimiter','\t','newline','pc')  %Line terminator, specified as the comma-separated pair consisting of 'newline' and either 'pc' to use a carriage return/line feed (CR/LF), or 'unix' to use a line feed (LF).





%% gap generation  based on the poisson distribution

w = Ecg1;
% amount of samples
N = length(w);
v = w;
L = length(v);
% estimate distibution parameters
lambda_gap = 70*rand(1000,1);
lambda_ok  = 700*rand(1000,1)+1100/45*lambda_gap ; %mean(lens_ok);1100/45 70 700
lambdagapandok=[lambda_gap lambda_ok];


figure(505)
scatter(lambda_gap,lambda_ok);
xlabel('lambda gap parameters'); ylabel('lambda ok parameters');
title('Gap percentage distribution parameters based on the VitalSens data ');


for s= 1: 1: length(lambda_gap)
    p= lambda_ok(s);
    q= lambda_gap(s);

% probability of a gap
p_gap = sum(v) / L;

% generate simulated gaps using exponential distribution for inter event
% intervals as the number of events is Poisson distributed

% result vector, true = gap, false = ok
u = zeros(size(v),'logical');

% starting state for the random generation
i = 1;
state = rand() < p_gap; % do we start with ok or gap?

% generate either ok or gap until enough samples as we wanted
while i < L
    if state == 0 % let's generate a ok run
        len = ceil( exprnd(p) );  % ceil(X) rounds each element of X and exprnd , generate exponential distrvution of random.
        j = min([i+len-1 L]);
        u(i:j) = false;
    else % let's generate a gap run
        len = ceil( exprnd(q) );
        j = min([i+len-1 L]);
        u(i:j) = true;
%         us(s,:)= u(i:j);
        
        
    end    
    
    i = j + 1;
    state = ~state;
end

% s(i,:)=double(u)
% one ble packet contains 4 samples so we need to upsample 

uu(s,:) = interp1(double(u) , 1:N, 'previous');


% calculate original signal gap lenghts % w diasilo kai
CW_w = bwconncomp( v, 4);
lens_w = cellfun( @length, CW_w.PixelIdxList );

% calculate simulated signal gap lenghts
CW_uu = bwconncomp( uu, 4);  %bwconncomp(BW) returns the connected components CC found in the binary image BW. 
lens_uu = cellfun( @length, CW_uu.PixelIdxList );



% show histograms

figure;
histogram(lens_w, 'facealpha',.5,'edgecolor','none')
hold on
histogram(lens_uu, 'facealpha',.5,'edgecolor','none')
set(gca,'YScale','log','XScale','log')
legend('original','generated');
title('Gap length distribution');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gap lenth and implement of PSD estimation method
size_uu= size(uu);
size_uu=size_uu(1,1);
ecg= Ecg1;
 for l=1 :1:size_uu
    uunan=find(uu(l,:)==1);
    
    %gap calculation
    length_uunan= length(uunan);
    length_uu= length(uu);
    gap_percentage(l,:)= length_uunan/ length(uu) *100 ;
    ecg(uunan)= nan;
    gapSignal(l,:) = ecg ;
    
    
%     figure
%     plot(gapSignal(l,:));
%     ecg=Ecg1;

    
 %%RR interval
 [qrspeaks,locs] = findpeaks(gapSignal(l,:),t1,'MinPeakHeight', 1100,...   % 1100 original 10 min
    'MinPeakDistance',0.05);

% figure
% subplot(211)
% plot(t1,gapSignal);
% hold on
% plot(locs,qrspeaks,'o');

R_interval = diff(locs);% for normal command
t_RR= locs(2:end);


% subplot(212)
% plot(t_RR, R_interval);
% title('RR interval for complete data ');


% ectoppic beat detection and replace it by NaN
dRRI=diff(locs,2);
ectopic_inds = dRRI > mad(dRRI,1)*6;
ectopic_inds( find(ectopic_inds)+1 ) = 1;
R_interval( ectopic_inds ) = nan;

% figure;
% plot(dRRI);
% hold on;
% plot( find(ectopic_inds), dRRI(ectopic_inds), 'rd')  % index wrong some 

% 
% figure
% subplot(311);
% plot(R_interval);


% Interpolation
t_int= t_RR(1): .25 : t_RR(end);
R_interval_intp = interp1(t_RR( ~ectopic_inds ),R_interval( ~ectopic_inds ),t_int,'spline');


% subplot(3,1,2);
% plot(t_int,R_interval_intp);
% title('Figure after interpolation to 4 Hz ')

sgf_signal= smooth(R_interval_intp, 20*4, 'sgolay',3);
detrended_RR_interv = R_interval_intp'-sgf_signal


% subplot(313);
% plot(t_int,detrended_RR_interv,'.-');
% title('Smoothing figure');
 
 
%LS calculation 
[pxx_lm F_lm]= plomb( R_interval( ~ectopic_inds ), t_RR( ~ectopic_inds ),F_lsorg );

%fft calculation
L = 360;  % 360 for 10min daata(234)
noverlap = 180;  %180 for 10min data
[signal,f_norm,pxxc] = pwelch(detrended_RR_interv,hann(L),noverlap,F_lsorg,4,...
    'onesided','ConfidenceLevel',0.95);
% f_norm = ((f_fft-min(f_fft)) ./ (max(f_fft)-min(f_fft))) ;


figure
subplot(3,1,1)
%plot(f_norm_fftorg(f_norm_fftorg>0.05),signal_org(f_norm_fftorg>0.05),'DisplayName','Original data spectrum');
plot(f_norm_fftorg(f_norm_fftorg>0.05 ),signal_org(f_norm_fftorg>0.05),'DisplayName','Original data spectrum');

hold on
%plot(f_norm(f_norm>0.05),signal(f_norm>0.05), 'DisplayName', ['Missing data percentage=' num2str(p*100)]);
plot(f_norm(f_norm>0.05 ),signal(f_norm>0.05 ),'DisplayName', ['Missing data percentage=' num2str( gap_percentage(l,:))]);

title('Power spectrum estimation via FFT base method For Normal Patient');
xlabel('frequency(Hz)');
ylabel('PSD-FFT');
% legend('Missing data spectrum');
legend show
ax = gca;
ax.YAxis.Exponent = 0;

f_fft_norm_perc(l,:) = (f_norm(f_norm>.05));
diff_perc_fft_signal(l,:) = (signal(f_norm>.05))

% f_fft_norm_perc(l,:)= f_norm((f_norm>.05)& (f_norm<.4))
% diff_perc_fft_signal(l,:) = signal((f_norm>.05) & (f_norm<.4))




%ar calculation
[Pxx_ar,F_ar_norm] = pburg(detrended_RR_interv,11,F_lsorg,4);
% F_ar_norm = ((F_ar-min(F_ar)) ./ (max(F_ar)-min(F_ar))) ;
 

subplot(3,1,2);
plot(F_ar_norm_org(F_ar_norm_org>0.05 ),Pxx_ar_org(F_ar_norm_org>0.05 ),'DisplayName','Original data spectrum');
hold on
plot(F_ar_norm(F_ar_norm>0.05 ),Pxx_ar(F_ar_norm>0.05 ),'DisplayName', ['Missing data percentage=' num2str( gap_percentage(l,:))]);
title('Power spectrum estimation via AR based method For Normal Patient');
xlabel('frequency(Hz)');
ylabel('PSD-AR');
% legend('Missing data spectrum');
legend show
ax = gca;
ax.YAxis.Exponent = 0;



for pp=1
F_ar_norm_perc{pp,l} = (F_ar_norm(F_ar_norm>.05));
diff_perc_ar_signal{pp,l} = (Pxx_ar(F_ar_norm>.05));
end

% for pp=1
% F_ar_norm_perc{pp,l} = F_ar_norm((F_ar_norm>.05) & (F_ar_norm<.4));
% diff_perc_ar_signal{pp,l} = Pxx_ar((F_ar_norm>.05) & (F_ar_norm<.4));
% end


%LS calculation
[pxx_lm F_lm]= plomb( R_interval( ~ectopic_inds ), t_RR( ~ectopic_inds ),F_lsorg );
 
 
subplot(313);
plot(F_lsorg(F_lsorg>0.05 ), pxx_lsorg(F_lsorg>0.05 ),'DisplayName','Original data spectrum');
hold on
plot(F_lm(F_lm>0.05 ), pxx_lm(F_lm>0.05 ) ,'DisplayName', ['Missing data percentage=' num2str( gap_percentage(l,:))]);
title('Power spectrum estimation vis LS(Lomb-Scargle) based method For Normal Patient');
xlabel('frequency(Hz)');
ylabel('PSD-LS');
% legend('Missing data spectrum');
legend show
ax = gca;
ax.YAxis.Exponent = 0;


F_lm_perc(l,:) = (F_lm(F_lm>.05));
diff_perc_lm_signal(l,:) = (pxx_lm(F_lm>.05));

% F_lm_perc(l,:) = F_lm((F_lm>.05) & (F_lm<.4));
% diff_perc_lm_signal(l,:) = pxx_lm((F_lm>.05) & (F_lm<.4));

ecg=Ecg1;

 end
 
 

 
 
 %% Differenc percentage calculation
 % fft method:
 
 for R=1:size_uu
    ref_signal= signal_org';
    size_didd_fft = size(diff_perc_fft_signal);
        if(length(ref_signal)>   size_didd_fft(1,2));
        p_diff_fft(R,:)= diff_perc_fft_signal(R,:)- ref_signal(1:length(diff_perc_fft_signal(R,:))); %;af er jonno boro soto look
     else 
         p_diff_fft(R,:)= diff_perc_fft_signal(R,(1:length(ref_signal)))- ref_signal ;%(1:length(diff_perc_fft_signal(R,:))) ;af er jonno boro soto look
     end
        p_diff_fft_sum(R,:)= sum(abs(p_diff_fft(R,:)));

    
     P_diff_fft_perc(R,:) = p_diff_fft_sum(R,:)/sum(ref_signal)*100; %perc store culumn wise for every rpeation for fixed percentage editdata
 end
 
 
 %ar method
   for R=1:size_uu
    ref_signal_ar= Pxx_ar_org';
    p_diff_ar{1,R}= diff_perc_ar_signal{1,R}- ref_signal_ar(1:length(diff_perc_ar_signal{1,R}))' ;
    
    p_diff_ar_sum{1,R}= sum(abs( p_diff_ar{1,R}));

    
     P_diff_ar_perc_cell{1,R} = (p_diff_ar_sum{1,R}/sum(ref_signal_ar)*100);
       
   end
 P_diff_ar_perc_cell1 =cell2mat(P_diff_ar_perc_cell) ;
 P_diff_ar_perc= P_diff_ar_perc_cell1';
 
  
 
 %ls method
  for R=1:size_uu
    ref_signal_ls= pxx_lsorg';
    p_diff_ls(R,:)= diff_perc_lm_signal(R,:)- ref_signal_ls(1:length(diff_perc_lm_signal(R,:))) ;
    
    p_diff_ls_sum(R,:)= sum(abs(p_diff_ls(R,:)));

    
     P_diff_ls_perc(R,:) = p_diff_ls_sum(R,:)/sum(ref_signal_ls)*100; %perc store culumn wise for every rpeation for fixed percentage editdata
 end
 
 
 
 
 %% entropy calculation
 
%fft entropy calculation

 P = diff_perc_fft_signal ./ sum( diff_perc_fft_signal, 2 );
 P_nan = nan( size(P) );
 P_nan( P>0 ) = P( P>0 );
 %  P_nan = P;
 %  P_nan( P<0.00001 ) = NaN;
 pse_rr_fft = - nansum( P.*log2(P), 2 );
 AV = sum( P.*f_fft_norm_perc, 2);
 

 
 
%ar entropy calculation
 for ka= 1:size_uu
 P = diff_perc_ar_signal {1,ka} ./ sum( diff_perc_ar_signal {1,ka} );
 
 P_nan = nan( size(P) );
 P_nan( P>0 ) = P( P>0 );


%  P_nan = P;
%  P_nan( P<0.00001 ) = NaN;
   pse_rr_ar(ka,:) = - nansum( P.*log2(P));
   
 end
 
 %ls entropy calculation
 P = diff_perc_lm_signal ./ sum( diff_perc_lm_signal, 2 );
 P_nan = nan( size(P) );
 P_nan( P>0 ) = P( P>0 );
%  P_nan = P;
%  P_nan( P<0.00001 ) = NaN;
   pse_rr_ls = - nansum( P.*log2(P), 2 );
 AV = sum( P.*f_fft_norm_perc, 2);
 
 
%  for kk= 1:size_uu
%  fresult1 =diff_perc_fft_signal(kk,:);
%     sum_fresult1=0.0;
%     for i=1:1:length(fresult1)-1
%      sum_fresult1= sum_fresult1 + (abs(fresult1(i)));
%     end
%     fresult1=fresult1/sum_fresult1;
%     entropy1=0.0;
%     for i=1:1:length(fresult1)-1
%         if  abs(fresult1(i)) > 0
%             entropy1_fft = entropy1 - abs(fresult1(i))*log2(abs(fresult1(i)));
%         end
%      pse_rr_fft(kk,:)=entropy1_fft;
%     end
%  end


%  for ka= 1:size_uu
%  fresult1=diff_perc_ar_signal{1,ka};
%     sum_fresult1=0.0;
%     for i=1:1:length(fresult1)-1
%      sum_fresult1= sum_fresult1 + (abs(fresult1(i)));
%     end
%     fresult1=fresult1/sum_fresult1;
%     entropy1=0.0;
%     for i=1:1:length(fresult1)-1
%      if  abs(fresult1(i)) > 0
%             entropy1_ar = entropy1 - abs(fresult1(i))*log2(abs(fresult1(i)));
%       end
%      pse_rr_ar(ka,:)= entropy1_ar;
%     end
%  end
 
 %ls entropy calculation
%   for kl= 1:size_uu
%  fresult1 =diff_perc_lm_signal(kl,:);
%     sum_fresult1=0.0;
%     for i=1:1:length(fresult1)-1
%      sum_fresult1= sum_fresult1 + (abs(fresult1(i)));
%     end
%     fresult1=fresult1/sum_fresult1;
%     entropy1=0.0;
%     for i=1:1:length(fresult1)-1
%       if  abs(fresult1(i)) > 0
%             entropy1_ls = entropy1 - abs(fresult1(i))*log2(abs(fresult1(i)));
%       end
%      pse_rr_ls(kl,:)=entropy1_ls;
%     end
%  end
 
  summary =[lambda_gap lambda_ok gap_percentage P_diff_fft_perc  P_diff_ar_perc P_diff_ls_perc pse_rr_fft pse_rr_ar pse_rr_ls ];


%% Figure betwwen gap percentage and signal difference

[percentage_fft  ii_fft]= sort(summary(:,3));
singal_difference_1=summary(:,4);
singal_difference_fft= singal_difference_1(ii_fft);

figure(506)
subplot(311)
plot(percentage_fft,singal_difference_fft,'*');
hold on
p_fft = polyfit(percentage_fft,singal_difference_fft,1);
f_fft = polyval(p_fft,percentage_fft);
plot(percentage_fft,f_fft,'r')
ylim([0 4000])
xlabel(' Missing data (%)'); ylabel('Data distribution');
title('Data distribution for normal patient based on the FFT method (1000 times gap generation)');
legend('FFT Method','Fitted line');

%FFT box plot part

i=0;
jj=1;
for ii= 200:200: length(singal_difference_fft)
%     for jj= 1:4 
    X_fft(:,jj)= singal_difference_fft(i+1:ii);
    X_fft_median_normal(:,jj) = median(X_fft(:,jj));
    
%        end
     i=ii;
     jj=jj+1;
    
end


dlmwrite(' X_fft_median_normal.txt', X_fft_median_normal, 'delimiter','\t','newline','pc')  %Line terminator, specified as the comma-separated pair consisting of 'newline' and either 'pc' to use a carriage return/line feed (CR/LF), or 'unix' to use a line feed (LF).



[percentage_ar  ii_ar]= sort(summary(:,3));
singal_difference_2=summary(:,5);
singal_difference_ar= singal_difference_2(ii_ar);
subplot(312)
plot(percentage_ar,singal_difference_ar,'*');
hold on
p_ar = polyfit(percentage_ar,singal_difference_ar,1);
f_ar = polyval(p_ar,percentage_ar);
plot(percentage_ar,f_ar,'r')
ylim([0 4000])
xlabel(' Missing data (%)'); ylabel('Data distribution');
title('Data distribution for normal patient based on the AR method (1000 times gap generation)');
legend('AR (11th order) Method','Fitted line');

%AR box plot part
i=0;
jj=1;
for ii= 200:200: length(singal_difference_ar)
%     for jj= 1:4 
    X_ar(:,jj)= singal_difference_ar(i+1:ii);
    X_ar_median_normal(:,jj) = median(X_ar(:,jj));
    
%        end
     i=ii;
     jj=jj+1;
    
end


dlmwrite('   X_ar_median_normal.txt',  X_ar_median_normal, 'delimiter','\t','newline','pc')  %Line terminator, specified as the comma-separated pair consisting of 'newline' and either 'pc' to use a carriage return/line feed (CR/LF), or 'unix' to use a line feed (LF).



subplot(313)
[percentage_ls  ii_ls]= sort(summary(:,3));
singal_difference_3=summary(:,6);
singal_difference_ls= singal_difference_3(ii_ls);
plot(percentage_ls,singal_difference_ls,'*');
hold on
p_ls = polyfit(percentage_ls,singal_difference_ls,1);
f_ls = polyval(p_ls,percentage_ls);
plot(percentage_ls,f_ls,'r')
ylim([0 4000])
xlabel(' Missing data (%)'); ylabel('Data distribution');
title('Data distribution for normal patient based on the LS method (1000 times gap generation)');
legend('LS Method','Fitted line');

%ls box plot part
i=0;
jj=1;
for ii= 200:200: length(singal_difference_ls)
%     for jj= 1:4 
    X_ls(:,jj)= singal_difference_ls(i+1:ii);
    
     X_ls_median_normal(:,jj) = median(X_ls(:,jj));
    
%        end
     i=ii;
     jj=jj+1;
    
end


 dlmwrite(' X_ls_median_normal.txt',  X_ls_median_normal, 'delimiter','\t','newline','pc')  %Line terminator, specified as the comma-separated pair consisting of 'newline' and either 'pc' to use a carriage return/line feed (CR/LF), or 'unix' to use a line feed (LF).


%% Box plot for mismatch calculation


figure(507)
subplot(311)
%boxplot([X_fft(:,1),X_fft(:,2),X_fft(:,3),X_fft(:,4),X_fft(:,5)],{'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'});
boxplot([X_fft(:,1),X_fft(:,2),X_fft(:,3),X_fft(:,4),X_fft(:,5)],{'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'});

ylim([0 3000])
% ylabel('Mismatch between specctrum'); 
ylabel('Error power (%)')
xlabel(' Missing data (%)')
title('Data distribution for normal patient based on the FFT method')
hLegend = legend(findall(gca,'Tag','Box'), {'FFT Method'});
% title('Data Distribution Comparison for Normal patients based on the FFT method (1000 times gap generation)')


subplot(312)
% boxplot([X_ar(:,1),X_ar(:,2),X_ar(:,3),X_ar(:,4),X_ar(:,5)],{'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'});
boxplot([X_ar(:,1),X_ar(:,2),X_ar(:,3),X_ar(:,4),X_ar(:,5)],{'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'});

ylim([0 3000])
% ylabel('Mismatch between specctrum'); 
ylabel('Error power (%)')
xlabel(' Missing data (%)')
% title('Data Distribution Comparison for Normal patients based on the AR method(1000 times gap generation)')
title('Data distribution for normal patient based on the AR method')
hLegend = legend(findall(gca,'Tag','Box'), {'AR(11th order) Method'});

subplot(313)
%boxplot([X_ls(:,1),X_ls(:,2),X_ls(:,3),X_ls(:,4),X_ls(:,5)],{'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'});
boxplot([X_ls(:,1),X_ls(:,2),X_ls(:,3),X_ls(:,4),X_ls(:,5)],{'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'});

ylim([0 3000])
% ylabel('Mismatch between specctrum'); 
ylabel('Error power (%)')
xlabel(' Missing data (%)')
% title('Data Distribution Comparison for Normal patients based on the LS method (1000 times gap generation)')
title('Data distribution for normal patient based on the LS method')

hLegend = legend(findall(gca,'Tag','Box'), {'LS Method'});


%% median calculation and represent it using  figure

figure(508)
y_values = [1 2 3 4 5]
y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};
plot(y_values,X_fft_median, 'Display', 'FFT method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};
plot(y_values,X_ar_median, 'Display', 'AR method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};
plot(y_values,X_ls_median, 'Display', 'LS method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);
legend('show');
xlabel('Missing data (%)');
ylabel('Medians');
title('Comparison between spectral based methods based on the median' );
ylim([80 180]);

%% statistical significance analysis using Ktest
% ktest % ranksumtest ranksum(x,y) also returns a logical value indicating the test decision. The result h = 1 indicates a rejection of the null hypothesis, 
%and h = 0 indicates a failure to reject the null hypothesis at the 5% significance level.

[h1,p1] = ranksum(X_fft(:,1),X_ar(:,1),'Alpha',0.05);
[h2,p2] = ranksum(X_fft(:,1),X_ls(:,1),'Alpha',0.05);
[h3,p3] = ranksum(X_ar(:,1),X_ls(:,1),'Alpha',0.05);

% [ktest_normalPatientsBin1] = [p1 p2; NaN p3];
[rank_normalPatientsBin1] = [h1 h2; NaN h3];


[h4,p4] =ranksum(X_fft(:,2),X_ar(:,2),'Alpha',0.05);
[h5,p5] = ranksum(X_fft(:,2),X_ls(:,2),'Alpha',0.05);
[h6,p6] = ranksum(X_ar(:,2),X_ls(:,2),'Alpha',0.05);

% [ktest_normalPatientsBin2] = [p4 p5; NaN p6];
[rank_normalPatientsBin2] = [h4 h5; NaN h6];

[h7,p7] = ranksum(X_fft(:,3),X_ar(:,3),'Alpha',0.05);
[h8,p8] = ranksum(X_fft(:,3),X_ls(:,3),'Alpha',0.05);
[h9,p9] = ranksum(X_ar(:,3),X_ls(:,3),'Alpha',0.05);

[rank_normalPatientsBin3] = [h7 h8; NaN h9];

[h10,p10] = ranksum(X_fft(:,4),X_ar(:,4),'Alpha',0.05);
[h11,p11] = ranksum(X_fft(:,4),X_ls(:,4),'Alpha',0.05);
[h12,p12] = ranksum(X_ar(:,4),X_ls(:,4),'Alpha',0.05);

[rank_normalPatientsBin4] = [h10 h11; NaN h12];


[h13,p13] = ranksum(X_fft(:,5),X_ar(:,5),'Alpha',0.05);
[h14,p14] = ranksum(X_fft(:,5),X_ls(:,5),'Alpha',0.05);
[h15,p15] = ranksum(X_ar(:,5),X_ls(:,5),'Alpha',0.05);

[rank_normalPatientsBin5] = [h13 h14; NaN h15];
 
 

%% % figure for entropy vs percentage ( ascending order)

lambda_gap_s=summary(:,1);
lambda_gap_sort= lambda_gap_s(ii_fft);

lambda_ok_s=summary(:,2);
lambda_ok_sort= lambda_ok_s(ii_fft);

pse_rr_fft_s=summary(:,7);
pse_rr_fft_sort= pse_rr_fft_s(ii_fft);

pse_rr_ar_s=summary(:,8);
pse_rr_ar_sort= pse_rr_ar_s(ii_fft);

pse_rr_ls_s=summary(:,9);
pse_rr_ls_sort= pse_rr_ls_s(ii_fft);

summary_sort =[lambda_gap_sort lambda_ok_sort percentage_fft singal_difference_fft  singal_difference_ar singal_difference_ls pse_rr_fft_sort pse_rr_ar_sort pse_rr_ls_sort ];
% % fid=fopen('abc.xls','wt');
% % fprintf(fid, '%d\n',  summary_sort);
% % fclose(fid);
dlmwrite('myFile_summary sort for boxentropy_normal.txt',summary_sort, 'delimiter','\t','newline','pc')  %Line terminator, specified as the comma-separated pair consisting of 'newline' and either 'pc' to use a carriage return/line feed (CR/LF), or 'unix' to use a line feed (LF).
% % dlmwrite('myFile.xls', summary_sort,'newline' ,'pc')


% figure for entropy vs percentage
figure(509)
subplot(311)
plot(summary_sort(:,3),summary_sort(:,7),'*');
hold on
p_fft_entr = polyfit(summary_sort(:,3),summary_sort(:,7),1);
f_fft_entr = polyval(p_fft_entr,summary_sort(:,3));
plot(summary_sort(:,3),f_fft_entr,'r')
xlabel('Percentage of Missing data'); ylabel('Entropy ');
title('Gap Percentage Vs Entropy for Normal patient');
legend('FFT Method','Fitted Line');

ax = gca;
ax.YAxis.Exponent = 0;
% ylim([0 0.04])



% %box plot
% boxplot([summary_sort(:,7)],{'Bin1( 0.1806% to 4.9639%)'});
% xlabel('Percentage of Missing data'); ylabel('Entropy ');
% title('Gap Percentage Vs Entropy for AF patient');
% legend('FFT Method');
% ax = gca;
% ax.YAxis.Exponent = 0;



subplot(312)
plot(summary_sort(:,3),summary_sort(:,8),'*');
hold on
p_ar_entr = polyfit(summary_sort(:,3),summary_sort(:,8),1);
f_ar_entr = polyval(p_ar_entr,summary_sort(:,3));
plot(summary_sort(:,3),f_ar_entr,'r')
xlabel('Percentage of Missing data'); ylabel('Entropy ');
title('Gap Percentage Vs Entropy for Normal patient');
legend('AR Method','Fitted Line');

ax = gca;
ax.YAxis.Exponent = 0;
% ylim([0 0.04])


subplot(313)
plot(summary_sort(:,3),summary_sort(:,9),'*');
hold on
p_ls_entr = polyfit(summary_sort(:,3),summary_sort(:,9),1);
f_ls_entr = polyval(p_ls_entr,summary_sort(:,3));
plot(summary_sort(:,3),f_ls_entr,'r')
xlabel('Percentage of Missing data'); ylabel('Entropy ');
title('Gap Percentage Vs Entropy for Normal patient');
legend('LS Method','Fitted Line');

ax = gca;
ax.YAxis.Exponent = 0;
% ylim([0 0.04])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Md. Surat-E- Mostafa %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




  


