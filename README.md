* The MCB is a r package which depends on `parallel`,`methods`,`leaps`,`lars`,`MASS`,`glmnet`,`ncvreg`,`parcor`,`flare`,`smoothmest`,`ggplot2` and `reshape2`.
* It includes function `mcb` and `mcb.compare`.

---

### Model Confidence Bound
*Description*  
>When choosing proper variable selection methods, it is important to consider the uncertainty of a certain method. The MCB for variable selection identifies two nested models (upper and lower confidence bound models) containing the true model at a given confidence level. A good variable selection method is the one of which the MCB under a certain confidence level has the shortest width. When visualizing the variability of model selection and comparing different model selection procedures, Model uncertainty curve is a good graphical tool. A good variable selection method is the one of whose MUC will tend to arch towards the upper left corner. This function aims to obtain the MCB and draw the MUC of certain single model selection method under a coverage rate equal or little higher than user-given confidential level.
  
*Usage*
>```mcb(x, y, B=200, lambda=NULL, method='Lasso', level=0.95, seed=122)```
  
*Arguments*
>`x` input matrix; each column is an observation vector of certain independent variable, and will be given a name automatically in the order of x1, x2, x3…  
>`y` y is a matrix of one column which presents the response vector B	number of bootstrap replicates to perform, default value is 200.  
>`lambda` A user supplied lambda value. It is the penalty tuning parameter for the variable selection method tested. The default value is the optimization outcome automatically computed in consideration of the specific case.  
>`method` Default value is ‘Lasso; user can choose from 'aLasso', 'Lasso', 'SCAD', 'MCP', 'stepwise', 'LAD', 'SQRT'.  
>`level` a positive value between 0 and 1, like the concept of confidence level for point estimation; Default value is 0.95.  
>`seed` seed for bootstrap procedures; Default value is 122.  
  
*Values*  
>`mcb` a list containing the bootstrap coverage rate (which is the closest to the user-given confidence level) and the corresponding model confidence bound of the user-chosen variable selection method in the form of lower confidence bound and upper confidence bound.  
>`mucplot` plot of the model uncertainty curve for this specific user-chosen variable selectionmethod.  
>`mcbframe` a dataframe containing all the information about MCBs for the specific variable selectionmethod under all bootstrap coverage rates including width(w), lower confidence bound(lcb) and upper confidence bound(ucb) for each bootstrap coverage rate(bcr)  

*Examples*
>```data(Diabetes) # load data```  
>```x <- Diabetes[,c('S1', 'S2', 'S3', 'S4', 'S5')]```  
>```y <- Diabetes[,c('Y')]```  
>```x <- data.matrix(x)```  
>```y <- data.matrix(y)```  
>```result <- mcb(x=x, y=y)```  
>```result$mucplot # plot of the model uncertainty curve```  
>```result$mcb # a list containing the bootstrap coverage rate and mcb```  
>```result$mcbframe # a dataframe containing all the information about MCBs```  

---

### Comparisons of Model Confidence Bounds for Different Variable selection Methods
*Description*  
>This function is a supplement of the function mcb. It is used to compare different variable selection methods and would return all the MUCs on same canvas. A good variable selection method’s MUC will tend to arch towards the upper left corner.
  
*Usage*
>```mcb.compare(x, y, B=200, lambdas=NULL, methods=NULL, level=0.95, seed=122)```
  
*Examples*
>```data(Diabetes) # load data```  
>```x <- Diabetes[,c('S1', 'S2', 'S3', 'S4', 'S5')]```  
>```y <- Diabetes[,c('Y')]```  
>```x <- data.matrix(x)```  
>```y <- data.matrix(y)```  
>```result <- mcb.compare(x=x, y=y)```  
>```result$mucplot # plot of the model uncertainty curves for all variable selection methods```  
>```result$mcb$Lasso # a list containing the bootstrap coverage rate and mcb which based on Lasso```  
>```result$mcbframe$Lasso # a dataframe containing all the information about MCBs which based on Lasso```  
