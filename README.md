* The MCB is a r package which depends on `parallel`,`methods`,`leaps`,`lars`,`MASS`,`glmnet`,`ncvreg`,`parcor`,`flare`,`smoothmest`,`ggplot2` and `reshape2`.
* It includes function `mcb` and `mcb.compare`.

---

### Model Confidence Bound
*Description*  
>When choosing proper variable selection methods, it is important to consider the uncertainty of a certain method. The MCB for variable selection identifies two nested models (upper and lower confidence bound models) containing the true model at a given confidence level. A good variable selection method is the one of which the MCB under a certain confidence level has the shortest width. When visualizing the variability of model selection and comparing different model selection procedures, Model uncertainty curve is a good graphical tool. A good variable selection method is the one of whose MUC will tend to arch towards the upper left corner. This function aims to obtain the MCB and draw the MUC of certain single model selection method under a coverage rate equal or little higher than user-given confidential level.
  
*Usage*
>```mcb(x, y, B=200, lambda=NULL, method='Lasso', level=0.95, seed=122)```
  
*Examples*
>```data(Diabetes) # load data```  
>```x <- Diabetes[,c('S1', 'S2', 'S3', 'S4', 'S5')```  
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
>```x <- Diabetes[,c('S1', 'S2', 'S3', 'S4', 'S5')```  
>```y <- Diabetes[,c('Y')]```  
>```x <- data.matrix(x)```  
>```y <- data.matrix(y)```  
>```result <- mcb.compare(x=x, y=y)```  
>```result$mucplot # plot of the model uncertainty curves for all variable selection methods```  
>```result$mcb$Lasso # a list containing the bootstrap coverage rate and mcb which based on Lasso```  
>```result$mcbframe$Lasso # a dataframe containing all the information about MCBs which based on Lasso```  
