%%
clear all
close all
clc

 
% fileID = fopen('testio.bin','rb');
% A = fread(fileID,'float');
quadrant = 0;
count = 0;
average = [];
output = [];
all_outputs = [];

path = 'C:\data\outdoor_meas_oct2016\Data_processing_for_Figures\650m';

output_650m_60mV_250khz = read_complex_binary('650m_250khz_60mV.bin',12182);
output_650m_60mV_250khz(abs(output_650m_60mV_250khz) > 12e-6) = [];
output_650m_60mV_250khz = output_650m_60mV_250khz(1:10000);

 
output_650m_62mV_250khz = read_complex_binary('650m_250khz_62mV.bin',12182);
output_650m_62mV_250khz(abs(output_650m_62mV_250khz) > 12e-6) = [];
output_650m_62mV_250khz = output_650m_62mV_250khz(1:10000);

 

 
output_650m_64mV_250khz =   read_complex_binary('650m_250khz_64mV.bin',12182);
output_650m_64mV_250khz(abs(output_650m_64mV_250khz) > 12e-6) = [];
output_650m_64mV_250khz = output_650m_64mV_250khz(1:10000);

 
output_650m_66mV_250khz = read_complex_binary('650m_250khz_66mV.bin',12182);
output_650m_66mV_250khz(abs(output_650m_66mV_250khz) > 12e-6) = [];
output_650m_66mV_250khz = output_650m_66mV_250khz(1:10000);

 

 
output_650m_68mV_250khz = read_complex_binary('650m_250khz_68mV.bin',12182);
output_650m_68mV_250khz(abs(output_650m_68mV_250khz) > 12e-6) = [];
output_650m_68mV_250khz = output_650m_68mV_250khz(1:10000);

 

 
output_650m_70mV_250khz = read_complex_binary('650m_250khz_70mV.bin',12182);
output_650m_70mV_250khz(abs(output_650m_70mV_250khz) > 12e-6) = [];
output_650m_70mV_250khz = output_650m_70mV_250khz(1:10000);

 

 
output_650m_72mV_250khz = read_complex_binary('650m_250khz_72mV.bin',12182);
output_650m_72mV_250khz(abs(output_650m_72mV_250khz) > 12e-6) = [];
output_650m_72mV_250khz = output_650m_72mV_250khz(1:10000);

 

 
output_650m_70mV_1Mhz = read_complex_binary('650m_1Mhz_70mV.bin',12182);
output_650m_70mV_1Mhz(abs(output_650m_70mV_1Mhz) > 12e-6) = [];
output_650m_70mV_1Mhz = output_650m_70mV_1Mhz(1:10000);

 

 
output_650m_noise_250khz = read_complex_binary('250k_650m_noise.bin',12182);
output_650m_noise_250khz(abs(output_650m_noise_250khz) > 12e-6) = [];
output_650m_noise_250khz = output_650m_noise_250khz(1:10000);

 

 
output_650m_noise_1Mhz =   read_complex_binary('650m_noise.bin',12182);
output_650m_noise_1Mhz(abs(output_650m_noise_1Mhz) > 12e-6) = [];
output_650m_noise_1Mhz = output_650m_noise_1Mhz(1:10000);

 

 
output_2_650m_70mV_250khz = read_complex_binary('2_650m_250khz_70mV.bin',12182);
output_2_650m_70mV_250khz(abs(output_2_650m_70mV_250khz) > 12e-6) = [];
output_2_650m_70mV_250khz = output_2_650m_70mV_250khz(1:10000);

 

 

 

 
[max_value_output_650m_60mV_250khz,max_index_output_650m_60mV_250khz] = max(abs(output_650m_60mV_250khz));
[max_value_output_650m_62mV_250khz,max_index_output_650m_62mV_250khz] = max(abs(output_650m_62mV_250khz));
[max_value_output_650m_64mV_250khz,  max_index_output_650m_64mV_250khz] = max(abs(output_650m_64mV_250khz));
[max_value_output_650m_66mV_250khz,max_index_output_650m_66mV_250khz] = max(abs(output_650m_66mV_250khz));
[max_value_output_650m_68mV_250khz,  max_index_output_650m_68mV_250khz] = max(abs(output_650m_68mV_250khz));
[max_value_output_650m_70mV_250khz,max_index_output_650m_70mV_250khz] = max(abs(output_650m_70mV_250khz));
[max_value_output_650m_72mV_250khz,max_index_output_650m_72mV_250khz] = max(abs(output_650m_72mV_250khz));
[max_value_output_650m_70mV_1Mhz,max_index_output_650m_70mV_1Mhz] = max(abs(output_650m_70mV_1Mhz));
[max_value_output_650m_noise_250khz,max_index_output_650m_noise_250khz] = max(abs(output_650m_noise_250khz));
[max_value_output_650m_noise_1Mhz,  max_index_output_650m_noise_1Mhz] = max(abs(output_650m_noise_1Mhz));
[max_value_output_2_650m_70mV_250khz,max_index_output_2_650m_70mV_250khz] = max(abs(output_2_650m_70mV_250khz));

 
max_complex_values = zeros(1,11);
max_complex_values(1) = output_650m_60mV_250khz(max_index_output_650m_60mV_250khz);
max_complex_values(2) = output_650m_62mV_250khz(max_index_output_650m_62mV_250khz);
max_complex_values(3) = output_650m_64mV_250khz(max_index_output_650m_64mV_250khz);
max_complex_values(4) = output_650m_66mV_250khz(max_index_output_650m_66mV_250khz);
max_complex_values(5) = output_650m_68mV_250khz(max_index_output_650m_68mV_250khz);
max_complex_values(6) = output_650m_70mV_250khz(max_index_output_650m_70mV_250khz);

 
max_complex_values(7) = output_650m_72mV_250khz(max_index_output_650m_72mV_250khz);
max_complex_values(8) = output_650m_70mV_1Mhz(max_index_output_650m_70mV_1Mhz);
max_complex_values(9) = output_650m_noise_250khz(max_index_output_650m_noise_250khz);
max_complex_values(10) = output_650m_66mV_250khz(max_index_output_650m_noise_1Mhz);
max_complex_values(11) = output_2_650m_70mV_250khz(max_index_output_2_650m_70mV_250khz);

 
all_outputs(1,:) = output_650m_60mV_250khz;
all_outputs(2,:) = output_650m_62mV_250khz;
all_outputs(3,:) = output_650m_64mV_250khz;
all_outputs(4,:) = output_650m_66mV_250khz;
all_outputs(5,:) = output_650m_68mV_250khz;
all_outputs(6,:) = output_650m_70mV_250khz;

 
all_outputs(7,:) = output_650m_72mV_250khz;
all_outputs(8,:) = output_650m_70mV_1Mhz;
all_outputs(9,:) = output_650m_noise_250khz;
all_outputs(10,:) = output_650m_noise_1Mhz;
all_outputs(11,:) = output_2_650m_70mV_250khz;
%%
outputs_processed=[];
for p = 1:11
   output= all_outputs(p,:); 
   max_complex = max_complex_values(p);
   phase_of_max = angle(max_complex);   
        if (phase_of_max > 0 && phase_of_max <= pi/2)
            quadrant = 1;
            for b = 1:length(output)
                if ((real(output(b)) > 0 && imag(output(b)) > 0) || (real(output(b)) < 0 && imag(output(b)) < 0) )
                    outputs_processed(p,b) = output(b);
                end
            end

            

            
        elseif (phase_of_max > pi/2 && phase_of_max <= pi)
            quadrant = 2;
            for b = 1:length(output)
                if ((real(output(b)) < 0 && imag(output(b)) > 0) || (real(output(b)) > 0 && imag(output(b)) < 0) )
                    outputs_processed(p,b) = output(b);
                end
            end

 
        elseif (phase_of_max > -pi/2 && phase_of_max <= 0)
            quadrant = 4;
            for b = 1:length(output)
                if ((real(output(b)) < 0 && imag(output(b)) < 0) || (real(output(b)) > 0 && imag(output(b)) > 0) )
                    outputs_processed(p,b) = output(b);
                end
            end
        else
            quadrant = 3;
            for b = 1:length(output)
                if ((real(output(b)) > 0 && imag(output(b)) < 0) || (real(output(b)) < 0 && imag(output(b)) > 0) )
                    outputs_processed(p,b) = output(b);
                end
            end
        end

           
end

 

 
%%
temp = 0;
Vpk = [];
Pdbm = [];
SNR = [];

 
for p = 1:11
        average = [];
        count = 0;
        output = outputs_processed(p,:);
        max_complex = max_complex_values(p);
        phase_of_max = angle(max_complex);

 
        if (phase_of_max > 0 && phase_of_max <= pi/2)
            quadrant = 1;
            for b = 1:length(output)
                if (real(output(b)) > 0 && imag(output(b)) > 0)
                    temp = abs(output(b));
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

 
%%
for a = 1:11
figure(a)
scatter(1e6*real(outputs_processed(a,:)),1e6*imag(outputs_processed(a,:)))
hold on
scatter(1e6*real(average_single(a)),1e6*imag(average_single(a)))
axis([-20 20 -20 20])
end
%%
figure(9)
scatter(1e6*real(outputs_processed(9,:)),1e6*imag(outputs_processed(9,:)))
hold on
scatter(1e6*real(average_single(9)),1e6*imag(average_single(9)),'g')
scatter(1e6*real(outputs_processed(5,:)),1e6*imag(outputs_processed(5,:)))
scatter(1e6*real(average_single(5)),1e6*imag(average_single(5)))
axis([-20 20 -20 20])



