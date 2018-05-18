# Run this file by calling on shell: Rscript dependencies_install.R

usePackage <- function(p)  {
	if (!is.element(p, installed.packages()[,1]))
		install.packages(p, dep = TRUE, repos='http://cran.us.r-project.org')
	
	suppressMessages(require(p, character.only = TRUE))
}

p <- c("devtools","Metrics","control")

for(i in p)
	usePackage(i)
