# SystemsIdentification-2017.2

Platform requirements: These scripts were tested on the operating systems described below.
- Ubuntu 14.04 or 16.04 (Recommended). 
- Ubuntu 18.04 (Eventually should require some additional dependencies).

Run the file "install_dependencies.sh" on the shell to fetch the required packages for the experiments:

		sh install_dependencies.sh

On directory At1, the description of files are showed below:
- The file "evaluation.R" runs the whole experiment. 
- The file "methods.R" implements the identification procedures.
- The file "filters.R" implements the filtering procedures.
- The file "systemblocks.R" implements generic blocks for transfer functions.
- Do the evaluation by typing on the shell the command below:

		Rscript evaluation.R

On directory At2, the description of files are showed below:
- The files "Q1", "Q2", "Q3", "Q4" implements the transfer functions from exercises. Run the evaluations by typing:

		Rscript Q1.R
		Rscript Q2.R
		Rscript Q3.R
		Rscript Q4.R
