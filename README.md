# vertical_greenhouse

execute main_xx.m files using matlab. Those files;
1-Perform light-energy calculations
2-create idf (EnergyPlus) files by using template idf files
3-run idf files using EnergyPlus which calculate thermal load of the greenhouse
4-read output htm files created by EnergyPlus and finally create output Tables for latitudes xx (22 41 and 60).

Solar data is 112 Mb so it is not possible to upload to GitHub. 
Please download them from:
https://www.refetyalcin.com/files/solar_data_22.mat
https://www.refetyalcin.com/files/solar_data_41.mat
https://www.refetyalcin.com/files/solar_data_60.mat
and put into the same folder with the files

prerequisites: -energyplus v24.1 -python v3.8 -install eppy package to python
if you do not want to use python wrapper, you can manually run generated idf files using EnergyPlus and read the EnergyPlus output data from .htm files 
do not hesitate to ask if you encounter a problem: refetali@gmail.com
