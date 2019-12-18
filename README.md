# Attention and optical imaging

By Ru-Yuan Zhang

Last updated: 12/



## Experiment and data details

* Two monkeys, Blofeld (B) and Kanaga (K).
* The experiment calls "Deviant"
* All data locate in 
  * /home/range1-raid2/sgwarren-data/Deviant/Blofeld/Data/Deviant/B
  * /home/range1-raid2/sgwarren-data/Deviant/Kananga/Data/Deviant/K

* *PositionTuning_AlignTo_stimulusOnTime.mat* is the data for **localizer** experiment
* *DeviantOff_AlignTo_stimulusOnTime.mat*，*DeviantOn_AlignTo_stimulusOnTime.mat* are the data for the **main** experiment.
* In all datasets, the images in different trials have already been preprocessed and aligned to stimulus onset.
* To access the data, you should at least be able to access to CMRR server 



## Other notes

* To run the code, you need to add utility functions by Ruyuan Zhang, https://github.com/ruyuanzhang/RZutil 
* To run the code, you also need the helper function by Scott, *ImageHelper.m*
* *Localizer_PCA.m* has been carefully documented. This is the program that performs PCA on localizer data.
* *mainExp_PCA\*\*.m* are programs that perform PCA on main experiment dataset
* ON and OFF indicate the main experiment and retrain datasets
* The *localizer_PCA.m* has been carefully documented. The mainExp* programs are almost the same.



