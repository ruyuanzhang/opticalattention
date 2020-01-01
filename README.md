# Attention and optical imaging

By Ru-Yuan Zhang

Last updated: 01/01/2020



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
* Relevant codes and figures locate at /home/stone/ruyuan/dropbox/stonesync/19attentionprobV1optimaging/opticalattention



## Figure

* \*localizer* refer to the localizer experiment. \*main\*  refer to the main experiment.

* \*map\* refer the data images. \*PC\* refer to the plot of two principle components.

* *Figure2.pdf* aims to replace the original figures in the manuscript. It contains the localizer maps for monkey B, and four PC plots for the two monkeys and the localizer and mainON experiments.

  

## Materials

* Two pdf files are printed emails from Scott. Scott explains the data details in these two emails.

## Other notes

* To run the code, you need to add utility functions by Ruyuan Zhang, https://github.com/ruyuanzhang/RZutil 
* To run the code, you also need the helper function by Scott, *ImageHelper.m*
* *Localizer_PCA.m* has been carefully documented. This is the program that performs PCA on localizer data.
* *mainExp_PCA\*\*.m* are programs that perform PCA on main experiment dataset
* ON and OFF indicate the main experiment and retrain datasets
* The *localizer_PCA.m* has been carefully documented. The mainExp* programs are almost the same.
* **DeviantON** is the dataset where the probability gradient was ON and the two positions 0 and 8 (1 or 9 in MATLAB) had a higher probability of being shown. **DeviantOFF** in Kananga was a dataset where the probability was even for all (e.g. the probability gradient was OFF). In Blofeld, DeviantOFF actually is a different control set where the probability gradient was FLIPPED to try to see larger training effects.



