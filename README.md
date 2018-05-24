# MSICorrect

MSICorrect is a Java-based commandline tool that performs recalibration of mass spectrometry imaging datasets without the need of specifying a lock-mass value. The computational approach implemented in the tool generates a reference spectrum, termed here as a _consensus spectrum_, using all the spectra present in the dataset. This _consensus spectrum_ is used as a reference spectrum to recalibrate the input mass spectra.

---

## Literature ##

[1] Purva Kulkarni, Filip Kaftan, Philipp Kynast, Aleš Svatoš and Sebastian Böcker  
**Correcting mass shifts: A lock-mass-free recalibration procedure for mass spectrometry imaging data.**  
_[Anal Bioanal Chem, 405(25):7603-7613, 2015.](http://link.springer.com/article/10.1007/s00216-015-8935-4)_  

[2] Sebastian Böcker and Veli Mäkinen  
**Combinatorial Approaches for Mass Spectra Recalibration.**  
_[IEEE/ACM Trans Comput Biology Bioinform, 5(1):91-100, 2008.](https://dx.doi.org/10.1109/tcbb.2007.1077)_  

---

## Download link ##

MSICorrect v1.0.0 can be downloaded as a jar file here

The complete source code can be found at [GitHub](https://github.com/purvakulkarni7/MSICorrect).

---

## Running MSICorrect

You can run the jar file using the commandline on any operating system with java installed. To run the tool, use the below command:

    java -jar /path/to/MSICorrect/MSICorrect.jar
    
    
    
    





