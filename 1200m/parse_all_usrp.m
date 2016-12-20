clear all;
close all
clc

% fileID = fopen('testio.bin','rb');
% A = fread(fileID,'float');
quadrant = 0;
count = 0;
average = [];
output = [];
all_outputs = [];

output_1200m_68mv_250k  = read_complex_binary('250k_68mv.bin',12288);
output_1200m_68mv_1M    = read_complex_binary('1M_68mv.bin',12288);
output_1200m_68mv_1770k = read_complex_binary('1770k_68mv.bin',12288);

output_1200m_60mv_2300k = read_complex_binary('2300k_60mv.bin',12288);
output_1200m_62mv_2300k = read_complex_binary('2300k_62mv.bin',12288);
output_1200m_64mv_2300k = read_complex_binary('2300k_64mv.bin',12288);
output_1200m_66mv_2300k = read_complex_binary('2300k_66mv.bin',12288);
output_1200m_68mv_2300k = read_complex_binary('2300k_68mv.bin',12288);
output_1200m_70mv_2300k = read_complex_binary('2300k_70mv.bin',12288);

noise_250k  = read_complex_binary('noise.bin',12288);
noise_1M    = read_complex_binary('1M_noise.bin',12288);
noise_1770k = read_complex_binary('1770k_noise.bin',12288);
noise_2300k = read_complex_binary('2300k_noise.bin',12288);


[max_value_output_1200m_68mv_250k  ,max_index_output_1200m_68mv_250k]  = max(abs(output_1200m_68mv_250k));
[max_value_output_1200m_68mv_1M    ,max_index_output_1200m_68mv_1M]    = max(abs(output_1200m_68mv_1M));
[max_value_output_1200m_68mv_1770k ,max_index_output_1200m_68mv_1770k] = max(abs(output_1200m_68mv_1770k));

[max_value_output_1200m_60mv_2300k ,max_index_output_1200m_60mv_2300k] = max(abs(output_1200m_60mv_2300k));
[max_value_output_1200m_62mv_2300k ,max_index_output_1200m_62mv_2300k] = max(abs(output_1200m_62mv_2300k));
[max_value_output_1200m_64mv_2300k ,max_index_output_1200m_64mv_2300k] = max(abs(output_1200m_64mv_2300k));
[max_value_output_1200m_66mv_2300k ,max_index_output_1200m_66mv_2300k] = max(abs(output_1200m_66mv_2300k));
[max_value_output_1200m_68mv_2300k ,max_index_output_1200m_68mv_2300k] = max(abs(output_1200m_68mv_2300k));
[max_value_output_1200m_70mv_2300k ,max_index_output_1200m_70mv_2300k] = max(abs(output_1200m_70mv_2300k));

[max_value_noise_250k              ,max_index_noise_250k]   = max(abs(noise_250k));
[max_value_noise_1M                ,max_index_noise_1M]     = max(abs(noise_1M));
[max_value_noise_1770k             ,max_index_noise_1770k]  = max(abs(noise_1770k));
[max_value_noise_2300k             ,max_index_noise_2300k]  = max(abs(noise_2300k));


max_complex_values = zeros(1,13);
max_complex_values(1) = output_1200m_68mv_250k(max_index_output_1200m_68mv_250k);
max_complex_values(2) = output_1200m_68mv_1M(max_index_output_1200m_68mv_1M);
max_complex_values(3) = output_1200m_68mv_1770k(max_index_output_1200m_68mv_1770k);

max_complex_values(4) = output_1200m_60mv_2300k(max_index_output_1200m_60mv_2300k);
max_complex_values(5) = output_1200m_62mv_2300k(max_index_output_1200m_62mv_2300k);
max_complex_values(6) = output_1200m_64mv_2300k(max_index_output_1200m_64mv_2300k);
max_complex_values(7) = output_1200m_66mv_2300k(max_index_output_1200m_66mv_2300k);
max_complex_values(8) = output_1200m_68mv_2300k(max_index_output_1200m_68mv_2300k);
max_complex_values(9) = output_1200m_70mv_2300k(max_index_output_1200m_70mv_2300k);

max_complex_values(10) = noise_250k  (max_index_noise_250k);
max_complex_values(11) = noise_1M    (max_index_noise_1M);
max_complex_values(12) = noise_1770k (max_index_noise_1770k);
max_complex_values(13) = noise_2300k (max_index_noise_2300k);


all_outputs(1,:) = output_1200m_68mv_250k;
all_outputs(2,:) = output_1200m_68mv_1M;
all_outputs(3,:) = output_1200m_68mv_1770k;

all_outputs(4,:) = output_1200m_60mv_2300k;
all_outputs(5,:) = output_1200m_62mv_2300k;
all_outputs(6,:) = output_1200m_64mv_2300k;
all_outputs(7,:) = output_1200m_66mv_2300k;
all_outputs(8,:) = output_1200m_68mv_2300k;
all_outputs(9,:) = output_1200m_70mv_2300k;

all_outputs(10,:) = noise_250k;
all_outputs(11,:) = noise_1M;
all_outputs(12,:) = noise_1770k;
all_outputs(13,:) = noise_2300k;

% scatter(real(noise_1770k),imag(noise_1770k))

%%
Vpk = [];
Pdbm = [];
SNR = [];

for p = 1:13
        average = [];
        count = 0;
        output = all_outputs(p,:);
        max_complex = max_complex_values(p);
        phase_of_max = angle(max_complex);
        if (phase_of_max > 0 && phase_of_max <= pi/2)
            quadrant = 1;
            for b = 1:length(output)
                if (real(output(b)) > 0 && imag(output(b)) > 0)
                    count = count+1;
                    average(b) = output(b);
                end
            end
            
            
        elseif (phase_of_max > pi/2 && phase_of_max <= pi)
            quadrant = 2;
            for b = 1:length(output)
                if (real(output(b)) < 0 && imag(output(b)) > 0)
                    count = count+1;
                    average(b) = output(b);
                 end
            end

        elseif (phase_of_max > -pi/2 && phase_of_max <= 0)
            quadrant = 4;
            for b = 1:length(output)
                if (real(output(b)) > 0 && imag(output(b)) < 0)
                    count = count+1;
                    average(b) = output(b);
                 end
            end
        else
            quadrant = 3;
                for b = 1:length(output)
                if (real(output(b)) < 0 && imag(output(b)) < 0)
                    count = count+1;
                    average(b) = output(b);
                 end
            end
        end
        average = average(abs(average)~=0);        
        average_single(p) = sum(average)/count;
        noise_magnitude(p) = sum(abs(average-average_single(p)))/count;
        amplitude(p) = abs(average_single(p));
        Vpk(p) = abs(average_single(p));
        Pdbm(p) = 20*log10(Vpk(p))-20;
        SNR(p) = 10*log10(Vpk(p).^2/noise_magnitude(p).^2);
end
Pdbm = transpose(Pdbm);
Vpk = transpose(Vpk);



%plot IQ plots for different bias levels
for a = 1:9
figure(a)
scatter(1e6*real(all_outputs(a,:)),imag(1e6*all_outputs(a,:)))
hold on
scatter(1e6*real(average_single(a)),imag(1e6*average_single(a)))
end

%plot noise
figure(11)
scatter(1e6*real(all_outputs(13,:)),imag(1e6*all_outputs(13,:)))
hold on
scatter(1e6*real(average_single(13)),imag(1e6*average_single(13)))


