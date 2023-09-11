%%%%%PROBLEM_1 
a=1+mod(282,3);
t = 0 : (2/120) : 2 - (2/120);
x = sin(2*pi*15*a*t);
figure
y1 = fft(x,120);
m1= abs(y1);
f1=  (0:length(y1)-1);
stem(f1,m1,'b')
hold on 
y2 = fft(x,130);
m2= abs(y2);
f2 =  (0:length(y2)-1)*120/length(y2);
stem(f2,m2,'r')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('FFT of Sine Wave')
legend('120 Samples', '130 Samples')
hold off
hold on
y3 = fft(x,60);
m3= abs(y2);
f3 =  (0:length(y2)-1)*120/length(y2);
stem(f3,m3,'r')




% %%%%% PROBLEM_2
A_values = 160;
B_values = 166;
fs = 200; % Sampling rate (samples per second)
duration = 10; % Duration of the signal (seconds)

% Time vector
t = 0:1/fs:duration;

% Generate the signal xa(t)
xa = 0.1 * sin(A_values * pi * t) + cos(B_values * pi * t);

% Calculate the DFT for different sample sizes
sample_sizes = [215, 415, 1115, 1515, 1915];

for j = 1:length(sample_sizes)
   N = sample_sizes(j)
   X = fft(xa(1:N), N); % Compute the DFT

   % Plot the DFT
   figure;
   stem(abs(X));
   title(['DFT of xa(t) for N = ' num2str(N)]);
   xlabel('Frequency Bin');
   ylabel('Magnitude');
end



%%%%PROBLEM_3


duration = 10;             % Duration in seconds
sampling_rate = 200;       % Sampling rate in samples/second
A = 160; % A values
B = 166; % B values
% Time vector
t = 0:1/sampling_rate:duration;
% Generate signal xa(t)
xa = 0.1 * sin(A * pi * t) + cos(B * pi * t);
%Sample the signal
sample_indices = [215, 415, 1115, 1515, 1915];
for sample_idx = 1:length(sample_indices)
    num_samples = sample_indices(sample_idx);
    sampled_xa = xa(1:num_samples);

    % window = hann(num_samples);
    window = blackman(num_samples);
    windowed_xa = sampled_xa .* window';

    % Compute the DFT
    X = fft(windowed_xa);

    % Plot the DFT
    figure;
    stem(abs(X));
    title(['DFT of windowed xa(t) for \alpha = ', num2str(3), ', A = ', num2str(A), ', B = ', num2str(B), ', ', num2str(num_samples), ' samples']);
    xlabel('Frequency Bin');
    ylabel('Magnitude');
    xlim([0, num_samples]);
end

%%%% % PROBLEM 4
load("Exp4Data1.txt")
x = Exp4Data1 ;
m = abs(fft(x,10000));
f = (0:length(m)-1)*1/length(m);
figure(1)
stem(f,m)
title("With Rectangular Window")
xlabel("Frequnecy");
ylabel("Amplitude");

w = hamming(500);
x1 = Exp4Data1.* w' ; 
m1 = abs(fft(x1,10000));
f1 = (0:length(m1)-1)*1/length(m1);
figure(2)
stem(f1,m1)
title("With hamming")
xlabel("Frequency");
ylabel("Amplitude");