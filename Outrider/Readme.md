# run.OUTRIDER timing

| code | user | system | elapsed |
| ---- | ---- | ------ | ------- |
|```ctsTable <- tibble::rownames_to_column(ctsTable, "geneID")```         |0.004| 0.000| 0.005|
|```ctsTable <- mutate(ctsTable, geneID = sub("\\..*$", "", geneID))```   |0.067| 0.000| 0.067|
|```ensembl <- useEnsembl()```|  1.86 |  0.026 | 48.025 |
|```queryResult <- getBM()``` | 1.962  | 0.086 | 58.715 |
|```ctsMatrix <- data.matrix(ctsTable)[, -1]```| 0.389 | 0.016 | 0.406 |
|```rownames(ctsMatrix) <- ctsTable[,1]```| 0.001 | 0.000 | 0.001 |
|```ods <- OutriderDataSet(countData=ctsMatrix)```| 0.685 | 0.040 | 0.812 |
|```ods <- filterExpression()```| 0.177 | 0.000 | 0.177 |
|```ods <- OUTRIDER(ods)```| 56058.445 |  819.045 | 4541.116 |

## OUTRIDER() log
Thu Mar 11 23:21:26 2021: SizeFactor estimation ...  
Thu Mar 11 23:21:27 2021: Controlling for confounders ...  
Using estimated q with: 12. 
Thu Mar 11 23:21:27 2021: Using the autoencoder implementation for controlling.  
[1] "Thu Mar 11 23:21:42 2021: Initial PCA loss: 3.2146859818928"  
[1] "Thu Mar 11 23:31:10 2021: Iteration: 1 loss: 2.67005548013455"  
[1] "Thu Mar 11 23:35:49 2021: Iteration: 2 loss: 2.66702081848771"  
[1] "Thu Mar 11 23:40:23 2021: Iteration: 3 loss: 2.66552423038741"  
[1] "Thu Mar 11 23:44:56 2021: Iteration: 4 loss: 2.66451398359503"  
[1] "Thu Mar 11 23:49:29 2021: Iteration: 5 loss: 2.6637467192021"  
[1] "Thu Mar 11 23:54:02 2021: Iteration: 6 loss: 2.66311766287492"  
[1] "Thu Mar 11 23:58:33 2021: Iteration: 7 loss: 2.66258021266031"  
[1] "Fri Mar 12 00:03:09 2021: Iteration: 8 loss: 2.66209903784366"  
[1] "Fri Mar 12 00:07:40 2021: Iteration: 9 loss: 2.66166819227166"  
[1] "Fri Mar 12 00:12:11 2021: Iteration: 10 loss: 2.66127176607511"  
[1] "Fri Mar 12 00:16:45 2021: Iteration: 11 loss: 2.66091296479433"  
[1] "Fri Mar 12 00:21:11 2021: Iteration: 12 loss: 2.66057895106051"  
[1] "Fri Mar 12 00:25:43 2021: Iteration: 13 loss: 2.66026332901442"  
[1] "Fri Mar 12 00:30:16 2021: Iteration: 14 loss: 2.65996577383506"  
[1] "Fri Mar 12 00:34:45 2021: Iteration: 15 loss: 2.65968076026424"  
Time difference of 1.141891 hours  
[1] "Fri Mar 12 00:34:45 2021: 15 Final nb-AE loss: 2.65968076026424"  
Fri Mar 12 00:34:46 2021: Used the autoencoder implementation for controlling.  
Fri Mar 12 00:34:46 2021: P-value calculation ...  
Fri Mar 12 00:37:07 2021: Zscore calculation ...  
