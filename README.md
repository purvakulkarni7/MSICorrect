# MSICorrect

MSICorrect is a Java-based commandline tool that performs recalibration of mass spectrometry imaging datasets without the need of specifying a lock-mass value. The computational approach implemented in the tool generates a reference spectrum, termed here as a _consensus spectrum_, using all the spectra present in the dataset. This _consensus spectrum_ is used as a reference spectrum to recalibrate the input mass spectra.

## Literature ##

[1] Purva Kulkarni, Filip Kaftan, Philipp Kynast, Aleš Svatoš and Sebastian Böcker  
**Correcting mass shifts: A lock-mass-free recalibration procedure for mass spectrometry imaging data.**  
_[Anal Bioanal Chem, 405(25):7603-7613, 2015.](http://link.springer.com/article/10.1007/s00216-015-8935-4)_  

[2] Sebastian Böcker and Veli Mäkinen  
**Combinatorial Approaches for Mass Spectra Recalibration.**  
_[IEEE/ACM Trans Comput Biology Bioinform, 5(1):91-100, 2008.](https://dx.doi.org/10.1109/tcbb.2007.1077)_  

## Download link ##

MSICorrect v1.0.0 can be downloaded as a jar file here

The complete source code can be found at [GitHub](https://github.com/purvakulkarni7/MSICorrect).

## Running MSICorrect

You can run the jar file using the commandline on any operating system with java installed. To run the tool, use the below command:

    java -jar /path/to/MSICorrect/MSICorrect.jar
    
MSICorrect does not require any input arguments while running the jar file. All the mandatory user input is taken one by one at the command prompt. You can always use the `--help` option to get a documentation about the necessary user input. To access the help, simply use the below command:

    java -jar /path/to/MSICorrect/MSICorrect.jar --help
    
Below is the information on mandatory user input:

1. **File path for input mass spectra**  
Diretory path that contains a list of peaklists as .txt files that have to be recalibrated. The peaklists should have coordinates values x and y separated by tab in the first line and then the mass and intensity value pairs separated by a tab from the next line and onwards.

2. **File path to store recalibrated spectra**  
Diretory path that  will store the recalibrated peaklists. Make sure that you have permissions to write the files
in the specified directory.

3. **Mass threshold for distance calculation between two peaklists**    
Suitable mass window (in Da) depending on the dataset, to calculate pairwise similarity between the two peaklists
(example: 0.5).

4. **Peaklist ordering approach (use TG or CG)**    
Before performing recalibration, it isimportant that all the peaklists are ordered in aspecific manner.
This especially is helpful to minimize the error within the growing consensus spectrum. MSICorrect implements 
two ordering approaches: type TG for topological greedy approach and CG for Crytal growth ordering approach.

5. **Line distance to perform line-pair stabbing recalibration**    
The Maximum Line-Pair Stabbing (MLS) algorithm performs outlier detection and recalibration of a pair of peaklists.
This parameter is a distance measure that separates the computational geometry interpretation of two parallel lines
in geometrical space (example: 0.5).

6. **Mass threshold for for line-pair stabbing recalibration**   
Suitable mass window (in Da) depending on the dataset, to perform recalibration of two peaklists (example: 0.8).

# Input file type
As an pinput, MSICorrect requires a directory path containing mass spectra. The spectra should be be preprocessed and peak-picked before being submitted for recalibration. Preprocessing and peak-detection of imaging MS data be performed using several packages available in Matlab or R. The peaklists generated should be exported as individual tab-separated text files with coordinate values for each spectra comprising the first line of the first. A short example of such a file is provided below:

    34	18
    101.7920	112.253425
    103.6095	1278.583291
    104.7244	145.525605
    105.6506	263.697254
    106.7781	129.368985
    107.7910	105.195891  
    108.7951	127.642998
    109.9640	157.623051
    110.7048	166.225015
    111.8753	118.473873
    112.8081	115.421191

# Sample data

Sample data to run MSICorrect can be found [here]().
    





    
    





