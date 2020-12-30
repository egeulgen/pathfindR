#This implementation is based on genalg R package, https://github.com/egonw/genalg which is Copyright (C) 2005 Egon Willighagen License: GPL 2.0

#TODO
#Parent selection probabilities

asGA<-function(chromosomes=10,
               popSize=200,
               iters=100,
               mutationChance=NA,
               elitism=NA,
               crossover='uniform', #implemented different from genalg
               suggestion=NULL, #implemented different from genalg
               suggestionRatio=1, #implemented different from genalg
               zeroToOneRatio=c(1,1), #implemented different from genalg
               evalFunc=NULL,
               parentProb= dnorm(1:popSize, mean=0, sd=(popSize/3))
){
  
  #Check input parameters
  if (is.null(evalFunc)) { stop("An evaluation function must be provided"); }
  if (chromosomes<1) { stop("The chromosome number must be at least 1."); }
  if (popSize < 5) { stop("The population size must be at least 5."); }
  if (iters < 1) { stop("The number of iterations must be at least 1."); }
  if (!is.na(elitism) & !(elitism < popSize)) { stop("The population size must be greater than the elitism."); }
  if (!crossover %in% c('uniform', 'singlepoint')) { stop("Crossover method parameter must be either \'uniform\' or \'singlepoint\'"); }
  if (length(zeroToOneRatio)!=2 | zeroToOneRatio[1]<=0 | zeroToOneRatio[2]<=0) { stop("zeroToOneRatio parameters must be a vector of two positive numbers, e.g. c(3,5) which means ~ 3 zeros for every 5 ones"); }
  if ((suggestionRatio < 0) | (suggestionRatio>1)) { stop("Suggestion ratio must be in [0,1] range."); }
  
  
  #Initialize parameters if not given
  if (is.na(mutationChance)) {
    mutationChance = 1/(chromosomes+1);
  }
  if (is.na(elitism)) {
    elitism = floor(popSize/5)
  }
  
  #Initialize population
  population = matrix(nrow=popSize, ncol=chromosomes)
  
  if (is.null(suggestion)){
    for (child in 1:popSize) {
      population[child,] = sample(c(rep(0,zeroToOneRatio[1]), rep(1, zeroToOneRatio[2])), chromosomes, replace=TRUE)
      while (sum(population[child,]) == 0) {
        population[child,] = sample(c(rep(0,zeroToOneRatio[1]), rep(1, zeroToOneRatio[2])), chromosomes, replace=TRUE)
      }
    }
  }else{
    initiallyOneNumber=round(length(suggestion)*suggestionRatio)
    for (child in 1:popSize) {
      population[child,] = rep(0, chromosomes)
      onesBasedOnSuggestion=sample(suggestion,initiallyOneNumber,replace=FALSE)
      population[child,onesBasedOnSuggestion]=1
      onesNumber=sum(population[child,])
      onesRequired=round((chromosomes/(zeroToOneRatio[1]+zeroToOneRatio[2]))*zeroToOneRatio[2])-onesNumber
      newOneIndices=sample(which(population[child,]==0), onesRequired, replace=FALSE)
      population[child,newOneIndices]=1
    }
  }

  #print(population)
  
  
  bestEvals = rep(NA, iters+1);
  meanEvals = rep(NA, iters+1);
  evalVals = rep(NA, popSize);
  
  
  for (iter in 1:iters) {
    
    # progress message
    if(iter %% 2 == 0) {
      message(round(iter / iters, 2) * 100, '% of total iterations in genetic algorithm finished')
    }
    
    #Evaluation of individuals
    for (object in 1:popSize) {
      if(is.na(evalVals[object])){
        evalVals[object] = evalFunc(population[object,]);  #sum(population[object,]);
      }
    }
    bestEvals[iter] = max(evalVals);
    meanEvals[iter] = mean(evalVals);
    evalVals = sort(evalVals, index=TRUE, decreasing = TRUE);
    population  = matrix(population[evalVals$ix,], ncol=chromosomes);
    
    
    #print(population)
    
    
    #Creating new population
    newPopulation = matrix(nrow=popSize, ncol=chromosomes);
    newEvalVals = rep(NA, popSize);
    
    #Adding elites to the new population
    if (elitism > 0) {
      newPopulation[1:elitism,] = population[1:elitism,];
      newEvalVals[1:elitism] = evalVals$x[1:elitism]
    }
    
    for (child in (elitism+1):popSize) {
      
      #Crossover
      
      #Pick two random parents
      parentIDs = sample(1:popSize, 2, prob= parentProb)
      parents = population[parentIDs,];
      
      if(crossover=='uniform'){#Uniform crossover
        newPopulation[child, ] = parents[1,]
        for (chr in 1:chromosomes){
          if(runif(1) < 0.5){
            newPopulation[child,chr]=parents[2,chr]
          }
        }
        while (sum(newPopulation[child,]) == 0) {
          newPopulation[child,] = sample(c(rep(0,zeroToOneRatio[1]), rep(1, zeroToOneRatio[2])), chromosomes, replace=TRUE)
        }
      }else{#Single point crossover
        crossOverPoint = sample(0:chromosomes,1);
        if (crossOverPoint == 0) {
          newPopulation[child, ] = parents[2,]
          newEvalVals[child] = evalVals$x[parentIDs[2]]
        } else if (crossOverPoint == chromosomes) {
          newPopulation[child, ] = parents[1,]
          newEvalVals[child] = evalVals$x[parentIDs[1]]
        } else {
          newPopulation[child, ] = c(parents[1,][1:crossOverPoint], parents[2,][(crossOverPoint+1):chromosomes])
          while (sum(newPopulation[child,]) == 0) {
            newPopulation[child,] = sample(c(rep(0,zeroToOneRatio[1]), rep(1, zeroToOneRatio[2])), chromosomes, replace=TRUE)
          }
        } 
      }
      
      #Mutation
      if (mutationChance > 0) {
        for (child in (elitism+1):popSize) {
          
          #For speed this version is chosen
          toMutate=sample(c(1:chromosomes),chromosomes*mutationChance, replace = FALSE)
          for(chr in toMutate){
            if (population[child,chr] == 0) {
              population[child,chr] = 1;
            } else {
              population[child,chr] = 0;
            }
          }
          
          # #This is slower          
          # for (chr in 1:chromosomes) {
          #   if (runif(1) < mutationChance) {
          #     ## OPTION 1: switch bit
          #     if (population[child,chr] == 0) {
          #       population[child,chr] = 1;
          #     } else {
          #       population[child,chr] = 0;
          #     }
          #     # OPTION 2: sample new bit with zeroToOneRatio change
          #     #population[child,chr] = sample(c(rep(0,zeroToOneRatio),1), 1);
          #   }
          # }
          
        }
      }
      

      
    }
    
    population = newPopulation;
    evalVals   = newEvalVals;
    
    
  }
  
  #Evaluation of individuals
  for (object in 1:popSize) {
    if(is.na(evalVals[object])){
      evalVals[object] = evalFunc(population[object,]);  #sum(population[object,]);
    }
  }
  bestEvals[iter+1] = max(evalVals);
  meanEvals[iter+1] = mean(evalVals);
  evalVals = sort(evalVals, index=TRUE, decreasing = TRUE);
  population  = matrix(population[evalVals$ix,], ncol=chromosomes);
  
  
  
  
  #print(population)
  print(bestEvals)
  print(meanEvals)
  
  return(population)
  
}