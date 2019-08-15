X_fft_median_normal= importdata('X_fft_median_normal.txt');
X_ar_median_normal= importdata('X_ar_median_normal.txt');
X_ls_median_normal= importdata('X_ls_median_normal.txt');

X_fft_median_af= importdata('X_fft_median_af.txt');
X_ar_median_af= importdata('X_ar_median_af.txt');
X_ls_median_af= importdata('X_ls_median_af.txt');



X_fft_median_tachy= importdata('X_fft_median_tachy.txt');
X_ar_median_tachy= importdata('X_ar_median_tachy.txt');
X_ls_median_tachy= importdata('X_ls_median_tachy.txt');

X_fft_median_brady= importdata('X_fft_median_brady.txt');
X_ar_median_brady= importdata('X_ar_median_brady.txt');
X_ls_median_brady= importdata('X_ls_median_brady.txt');



figure(508)
subplot(411)
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};
% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};


plot(y_values,X_fft_median_normal,'-v', 'Display', 'FFT method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};

plot(y_values,X_ar_median_normal,'-v', 'Display', 'AR method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_ls_median_normal,'-v', 'Display', 'LS method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);
legend('show');
xlabel('Missing data (%)');
ylabel('Medians');
% title('Comparison between spectral based methods based on the median value (normal patient)' );
title(' Normal patient' );

ylim([20 200]);


%
subplot(412)
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_fft_median_af,'-v', 'Display', 'FFT method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_ar_median_af,'-v', 'Display', 'AR method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_ls_median_af,'-v', 'Display', 'LS method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);
% legend('show');
xlabel('Missing data (%)');
ylabel('Medians');
% title('Comparison between spectral based methods based on the median value (AF patient)' );

title(' Atrial fibrillation patient' );


ylim([20 200]);

%
subplot(413)
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_fft_median_tachy,'-v', 'Display', 'FFT method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_ar_median_tachy,'-v', 'Display', 'AR method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};
y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_ls_median_tachy,'-v', 'Display', 'LS method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);
% legend('show');
xlabel('Missing data (%)');
ylabel('Medians');
% title('Comparison between spectral based methods based on the median value (tachycardia patient)' );
title(' Tachycardia patient' );

ylim([20 200]);
%
subplot(414)
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};
y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_fft_median_brady,'-v', 'Display', 'FFT method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};
y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_ar_median_brady,'-v', 'Display', 'AR method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);

hold on
y_values = [1 2 3 4 5]
% y_labels ={'Bin1(0.1954% to 2.1046% )','Bin2(2.1102% to 2.5972%)','Bin3(2.6018% to 3.0620%)','Bin4(3.0639% to 3.5926%)','Bin5(3.5990% to 5.5999%)'};

% y_labels ={'Bin1(0.20% to 2.10% )','Bin2(2.11% to 2.59%)','Bin3(2.60% to 3.06%)','Bin4(3.07% to 3.59%)','Bin5(3.60% to 5.60%)'};

y_labels ={'Bin1','Bin2','Bin3','Bin4','Bin5'};
plot(y_values,X_ls_median_brady,'-v', 'Display', 'LS method' );
set(gca, 'Xtick',y_values,'XTickLabel',y_labels);
% legend('show');
xlabel('Missing data (%)');
ylabel('Medians');
% title('Comparison between spectral based methods based on the median value (bradycardia patient)' );
title('Bradycardia patient' );

ylim([20 200]);


