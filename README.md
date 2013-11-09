# Maxillaria

A light-weight nested sampler. The code is based on the example given in
Sivia and Skilling 2006, Data Analysis, A Bayesian Tutorial, Second Edition.
I have made the interface similar to the widely used nested-sampler MultiNest (Feroz et al 2008).

This is a work in progress... The code will have many bugs. If you find one please
report using issue tracking.

## Author

1. Sreekumar Thaithara Balan (tbs1980@gmail.com)

## Compilation
Unpack the Maxillaria files

	$ tar xvzf Maxillaria.tar.gz
	
Make a direcotry to build the files. Let's call it `build`. Then change directory to `build`

	$ cd Maxillaria
	$ mkdir build
	$ cd build
	
Now we are ready to build, type

	$ cmake ../
	$ make
	
This will create a static library called `libmaxillaria.a` and an example excutable `unit_gauss.exe`.
	
## Examples
A unit-gaussian posterior distribution case can be found in examples direcotry. 
First make a direcotry called chains and the run the excutable.
	$ mkdir chains
	$ ./examples/unit_gauss.exe
	
The output will look like

	-------------------------------------------------
		     Sree's Nested Sampler
		        Version 1.0
	   Sreekumar Thaithara Balan (tbs1980@gmail.com)
	-------------------------------------------------
	Starting sampling from scratch...

	Number of iterates    = 100
	Evidence,log(Z)        = -144.73 +- 0.0688597
	Dlog(Z)                = 0.25315

	Number of iterates    = 200
	Evidence,log(Z)        = -121.72 +- 0.0659767
	Dlog(Z)                = 0.14761

	Number of iterates    = 300
	Evidence,log(Z)        = -110.338 +- 0.0640854
	Dlog(Z)                = 0.216004

	Number of iterates    = 400
	Evidence,log(Z)        = -98.297 +- 0.06388
	Dlog(Z)                = 0.106163
	
	................................................
	................................................
	................................................

	Number of iterates    = 7300
	Evidence,log(Z)        = -5.15544 +- 0.0634252
	Dlog(Z)                = 0.000105041

	Number of iterates    = 7400
	Evidence,log(Z)        = -5.14546 +- 0.0634269
	Dlog(Z)                = 9.48654e-05

	-------------------------------------------------
	Sampling finished!
	Iteraions          = 7400
	Evidence:log(Z)    = -5.14546 +- 0.0634269
	Information: H     = 4.022972 nats(= 5.803922 bits)
	-------------------------------------------------
	
## Output files
The output files can be found in chains direcotry. 

	$ ls -l chains/
	-rw-r--r-- 1 sree sree 592079 Nov  9 10:03 test.extract.txt
	-rw-r--r-- 1 sree sree  79002 Nov  9 10:03 test.phys_live_points.txt
	-rw-r--r-- 1 sree sree 289270 Nov  9 10:03 test.post_equal_weights.txt
	-rw-r--r-- 1 sree sree    676 Nov  9 10:03 test.stats.txt

The *.stats.txt contains the posterior estimates.
	





