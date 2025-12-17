# Homogeneity_Test_Combined_Data
The source code is contained in folder <homogeneity_test> written in Fortran, assuming default compiler <gfortran>

Assume the current workdirectroy is: path-to/homogeneity_test/

In order to make the code fully compilable, the external LAPACK package needs be compiled first. 
Release the zipped file <lapack-3.12.1.tar.zip> and suppose the released folder name <1apack-3.12.1> and is under path-to/homogeneity_test/.
Follow the instructions inside (README.md) to compile the code. 

## To reproduce the simulation study presented in the manuscript, do the following: 
> cd run
> gfortran -o gen_tasks.f -o gen_tasks.exe
> ./gen_tasks.exe

Then the following prompts will appear: 
> generate configurations? -- enter "Y" or "y" for generating all sample size and parameter configurations; "N" or "n" for not doing so
If yes, then
> What is the model?
> type 'Don' for Donner's model; 'Ros' for Rosner's model
> empirical type I error rates or powers?
> type 'H0' for TIE rates; 'H1' for powers -- enter "H0" for TIE; "H1" for power
> Compile excutables w/ corresponding # of groups??? -- enter "Y" or "y" for compiling source code with g=2,...,8, and create symbolic links for each excutable with proper # of groups in folder <run_xxx_model_Hx>

Execute ./gen_tasks.exe again to run the simulation for MLE-based tests meanwhile creating SAS files ready to be run in SAS
> ./gen_tasks.exe
> generate configurations? -- enter "N" or "n"
> ...no generating configurations...
> run excutives? -- enter "Y" or "y" to run simulation and SAS code generation

After running excutables, the simulation results are inside folder <fort> and the SAS codes are inside folder <sas>

## To reproduce the real analysis
go to folder <src> and run
> ./run_real_data_ana.sh

Then run excutable
> real_data_ana.exe
