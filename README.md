# 1kmRFID
Long range Quantum RFID with tunnel diodes


Long Range measurements data analysis with RFID quantum tunneling tag (QTT). Data acquisition and processing. 
Authors: Francesco Amato (f.amato@gatech.edu) and Hakki Mert Torun (htorun3@gatech.edu) Georgia Institute of Technology
December 2016

Readme file structure
1. *The Experimental Setups*
2. *Folders and Files descriptions*
3. *Figures*



1. *The Experimental Setups*

The *.bin data have been acquired during the measured campaign held at the Georgia Institute of Technology during the Fall2016 semester.

Several experimental setups and settings have been used.

Exp. Setup 1:
Reader transmitting 0 dBm EIRP power
Transmitting antenna gain: 6 dBi
Receiving antenna gain: 24 dBi
QTT antenna gain: 6 dBi


2. *Folders and Files descriptions*

- 25TO50m: 
this folder contains all the *.bin files directly acquired from the reader (Exp. Setup 1) during the field measurement on the Van Leer rooftop at the Georgia Insitute of Technology the 19th of October 2016 and the *.m files to process and plot the data.
The file names structure is the following: distance_bias_modulation.bin.
For example, the wireless measurement taken at 45 meters away from the reader with the QTT being biased at 62 mV and modulating at 250 kHz is saved as 45m_62mv_250khz.bin

- 70TO160m: 
this folder contains all the *.bin files directly acquired from the reader (Exp. Setup 1) during the field measurement on the Tech Green Walk at the Georgia Insitute of Technology the 20th of October 2016.
The file names structure is the following: distance_bias_modulation.bin.
For example, the wireless measurement taken at 70 meters away from the reader with the QTT being biased at 60 mV and modulating at 1 MHz is saved as 70m_60mv_1Mhz.bin

- parse_all_usrp.m:
this m.file processes the data contained in the folders listed above and generates four plots: 
Plot 1 shows the I and Q data received on the reader from the QTT at 25 m and 45 m, (+ avarages) at 250 kHz 
Plot 2 shows the Power received on the reader vs distances
Plot 3 shows I and Q from 100 and 160 m, (+ avarages) at 250 kHz 
Plot 4 shows the I and Q from 25 and 45 m and 160 m (+ avarages) at 1 MHz 

- Figures:
this folder contains the .eps files generate by the m.file above.

3. *Figures*

- Fig_IQ25-45-250khz.eps: IQ diagrams for the received backscattered signal at 250 kHz modulation. QTT at 25 and 45 meters away from the reader
- Fig_IQ100-160-250khz.eps: IQ diagrams for the received backscattered signal at 250 kHz modulation. QTT at 100 and 160 meters away from the reader
- Fig_IQ25-45-160-1mhz.eps: IQ diagrams for the received backscattered signal at 1 MHz modulation. QTT at 25, 45 and 160 meters away from the reader
- Fig_Wireless_ranges.eps: Received signal strengths as function of distances, biasing voltages and modulation speeds.








