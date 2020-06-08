#This program is written in R, which is for prediction of abiotic stress-responsive proteins
# associated with eight abiotic stresses i.e., cold, drought, heat, light, salt, oxidative, phyto and wound.

Following are the steps for running the code.

1. This program is dependent upon window operating sytem (Windows 7 or above)

2. create a folder in your system.

3. Download all the files and put them in the created folder.

4. Download compatibele version of ModPred.exe from http://www.modpred.org/ and place in that folder (kindly put only.exe file)

5. Install suitable version of MCRInstaller available at http://www.modpred.org/.

6. Download R-3.3.0 window version in your system.

7. Install the required packages i.e., Biostrings, R2HTML, RSNNS, e1071

8. Set the working directory at three different places as marked in asrpro.R file.

9. Create the test dataset and name it as test.fasta, whre the sequences would be in fasta format in that file.

10. Kindly Unzip pscdnc.zip file and put them in the same folder name.

10. Provide sequences only with standard amino acid residues.

11. Open R-GUI and run the r-code given in asrpro.R

12. A result file asrpro.txt will be craeted in the folder, which provides the prediction probability of 8 different abiotic stresses categories.

13. If you are facing any problem, donot hesiste to write at meherprabin@yahoo.com

14. all the datasets used in thsi study are provided in training_dataset.zip and test_dataset.zip.