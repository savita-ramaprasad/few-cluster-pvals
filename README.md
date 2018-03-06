#### Calculate Cluster-Robust p-Values and Confidence Intervals Based on Webb (2013)

##### Description
 Calculate cluster robust p-values and confidence intervals using wild cluster bootstrapped t-statistics based on Webb (2013) which is the prefered method to use when the number of clusters are < 15. Webb, M. D. (2013). <Reworking wild bootstrap based inference for clustered errors (No. 1315). Queen's Economics Department Working Paper>
 
##### Usage
cluster.webb.glm(mod, dat, cluster, vars.boot = NULL, 
ci.level = 0.95, impose.null = TRUE, boot.reps = 1000, 
report = TRUE, prog.bar = TRUE, output.replicates = FALSE) 



##### Arguments

  **mod** A linear model estimated using glm.
      
  **dat** the data set used to estimate mod.   
   
  **cluster** The clustering variable.    
  
  **vars.boot** The variables to bootstrap over when null is imposed. Default is p-values for all variables. Displays p-values for all variables with null not imposed.  
  
  **ci.level** The confidence level of the confidence interval. Reported when impose.null = FALSE
    
  **impose.null** Should null be imposed?  
  
  **boot.reps** The number of bootstrap repititions.  
  
  **report** Report the result to the console?  
  
  **prog.bar** Show a progress bar of the bootstrap?  
  
  **output.replicates** Should the cluster bootstrap replicates be outputted as well?  



##### Value

**p.values** A vector of estimated p-values.
	
**ci** A matrix of confidence intervals, reported when null is not imposed.


##### References
	
Esarey, Justin, and Andrew Menger. 2017. "Practical and Effective Approaches to Dealing with Clustered Data." Political Science Research and Methods forthcoming: 1-35. <URL:http://jee3.web.rice.edu/cluster-paper.pdf>.  
	
Reworking wild bootstrap based inference for clustered errors (No. 1315). Queen's Economics Department Working Paper
	
