package ActiveSubnetworkSearchAlgorithms;

import ActiveSubnetworkSearchMisc.ScoreCalculations;
import ActiveSubnetworkSearchMisc.Subnetwork;
import Application.Parameters;
import Network.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author Ozan Ozisik
 */

enum SelectionType {RANKSELECTION, ROULETTEWHEEL}; 
enum CrossoverType {SINGLEPOINT, MULTIPOINT, UNIFORM};

public class GeneticAlgorithm {
    
    ScoreCalculations scoreCalculations;
    ArrayList<Node> networkNodeList;
    Random random;
    int populationSize;
    
     
    
    public ArrayList<Subnetwork> geneticAlgorithm(){
        scoreCalculations=ActiveSubnetworkSearch.scoreCalculations;
        networkNodeList = ActiveSubnetworkSearch.networkNodeList;
        populationSize=Parameters.ga_populationSize;
        ArrayList<GAIndividual> population=new ArrayList<>();
        random=new Random();
        
        initializePopulation(population, populationSize);
        
        printSituation(population);
        
        int addRandomIndividualCount=0;

        boolean running=true;
        int iter=0;
        GAIndividual lastBestIndividual=population.get(0);
        int lastBestRepeatNumber=0;
        while(running){
            
            /**
             * New population created
             */
            ArrayList<GAIndividual> newPopulation=createNewPopulation(population
                    , SelectionType.RANKSELECTION, CrossoverType.UNIFORM);
            
            /**
             * After each 10 steps 10% of the population (worst scoring part)
             * is replaced with random individuals
             */
            if(addRandomIndividualCount==10){
                for(int i=1;i<=(int)Parameters.ga_populationSize*0.1;i++){
                    newPopulation.set(newPopulation.size()-i, createRandomGAIndividual());
                }
                addRandomIndividualCount=0;
            }
            
            
            /**
             * Best scoring individuals are checked to prevent score decrement
             * in the next population. 
             * There is a possibility that this may overwrite one of
             * the randomly added individuals above. It is not a big deal.
            */
            if(Parameters.ga_Elitism){
                if(newPopulation.get(0).compareTo(population.get(0))<0){
                    newPopulation.set(newPopulation.size()-1, population.get(0));
                }
            }
            
            Collections.sort(newPopulation,Collections.reverseOrder());
            
            population=newPopulation;
            
            System.out.println("New Population, iter="+iter);
            
            printSituation(population);
            
            //TODO: Be careful about GAIndividuals with 0 subnetworks, you may 
            //consider adding empty subnetwork with score 0 in subnetwork finder
            //when no subnetwork found
            
            addRandomIndividualCount++;
            
            iter++;
            
            if(population.get(0).compareTo(lastBestIndividual)==0){
                lastBestRepeatNumber++;
            }else{
                lastBestIndividual=population.get(0);
                lastBestRepeatNumber=0;
            }
            
            if(lastBestRepeatNumber>=50){
                running=false;
                System.out.println("The score did not improve in 50 steps");
            }
            if(iter>=Parameters.ga_totalIterations){
                running=false;
            }
            
        }
        
        return population.get(0).getSubnetworkList();

    }
    
    private void printSituation(ArrayList<GAIndividual> population){
        for (int i = 0; i < 10; i++) {
            System.out.println(population.get(i).toString());
        }
    }
    
    
    private void initializePopulation(ArrayList<GAIndividual> population, int populationSize){
        for (int i = 0; i < populationSize; i++) {
            population.add(createRandomGAIndividual());
        }

        /**
         * Creates a GAIndividual that contains all the genes that have 
         * positive scores
         */ 
        if(Parameters.startWithAllPositiveZScoreNodes){
            ArrayList<Boolean> individualPositiveZ=new ArrayList<>();
            for(int i=0;i<networkNodeList.size();i++){
                Node node=networkNodeList.get(i);
                individualPositiveZ.add(scoreCalculations.getZScore(node) > 0);
            }
            population.set(populationSize-1, new GAIndividual(individualPositiveZ));
        }
        
        Collections.sort(population,Collections.reverseOrder());
    }
    
    private GAIndividual createRandomGAIndividual(){
        ArrayList<Boolean> individual = new ArrayList<>();
        for (int i = 0; i < networkNodeList.size(); i++) {
            individual.add(random.nextDouble()<Parameters.geneInitialAdditionProbability);
        }
        return new GAIndividual(individual);
    }
    
    private ArrayList<GAIndividual> createNewPopulation(ArrayList<GAIndividual> population, SelectionType selectionType, CrossoverType crossoverType){
        ArrayList<GAIndividual> newPopulation=new ArrayList<>();
        
        long start=System.nanoTime();
        ArrayList<Thread> threads=new ArrayList();
        for(int i=0;i<Parameters.ga_threadNumber;i++){
            Thread thread=new Thread(new NewPopulationFactory(population, 
                    newPopulation, selectionType, crossoverType));
            threads.add(thread);
            thread.start();
        }
        
        for(Thread thread:threads){
            try {
                thread.join();
            } catch (InterruptedException ex) {
                Logger.getLogger(GeneticAlgorithm.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        long stop=System.nanoTime();
        System.out.println("Time:"+(stop-start)/1000000);

        
        
        
        Collections.sort(newPopulation,Collections.reverseOrder());
        return newPopulation;
    }
    
}

class NewPopulationFactory implements Runnable{

    ArrayList<GAIndividual> newPopulation;
    ArrayList<GAIndividual> population;
    SelectionType selectionType;
    CrossoverType crossoverType;
    
    public NewPopulationFactory(ArrayList<GAIndividual> population, ArrayList<GAIndividual> newPopulation, SelectionType selectionType, CrossoverType crossoverType){
        this.population=population;
        this.newPopulation=newPopulation;
        this.selectionType=selectionType;
        this.crossoverType=crossoverType;
    }
    
    @Override
    public void run() {
        while(newPopulation.size()<population.size()){
            GAIndividual[] parents = selection(population);
            GAIndividual[] offsprings = crossoverAndMutation(parents[0], parents[1]);

            synchronized(newPopulation){
                if(newPopulation.size()<population.size()){
                    for(int i=0;i<offsprings.length;i++){
                        newPopulation.add(offsprings[i]);
                    }
                }
            }
        }
    }
    
    private GAIndividual[] selection(ArrayList<GAIndividual> population) {
        double totalWeight = 0;
        double weights[] = new double[population.size()];

        if (selectionType == SelectionType.RANKSELECTION) {
            for (int i = 0; i < population.size(); i++) {
                weights[i] = population.size()-i;
                totalWeight = totalWeight + weights[i];
            }
        } else {//SelectionType.ROULETTEWHEEL, individuals who have the same 
            //highest scoring subnetwork will have the same weight
            for (int i = 0; i < population.size(); i++) {
                weights[i] = population.get(i).getScore();
                totalWeight = totalWeight + weights[i];
            }
        }

        GAIndividual[] parents=new GAIndividual[2];
        for (int i = 0; i < 2; i++) {
            int randomIndex = -1;
            double rand = ThreadLocalRandom.current().nextDouble() * totalWeight;
            int rr = 0;
            while ((rr < population.size()) && (randomIndex == -1)) {
                rand = rand - weights[rr];
                if (rand <= 0.0) {
                    randomIndex = rr;
                }
                rr++;
            }
            parents[i]=population.get(randomIndex);
        }
        
        return parents;
    }
    
    private GAIndividual[] crossoverAndMutation(GAIndividual parent1, GAIndividual parent2){
        ArrayList<Boolean> parent1Boolean=parent1.getRepresentationBoolean();
        ArrayList<Boolean> parent2Boolean=parent2.getRepresentationBoolean();
        ArrayList<Boolean> child1Boolean=new ArrayList<>();
        ArrayList<Boolean> child2Boolean=new ArrayList<>();
        
        /**
         * Crossover
         */
        if(crossoverType==CrossoverType.SINGLEPOINT){
            int crossoverPoint=ThreadLocalRandom.current().nextInt(parent1Boolean.size());
            for(int i=0;i<crossoverPoint;i++){
                child1Boolean.add(parent1Boolean.get(i));
                child2Boolean.add(parent2Boolean.get(i));
            }
            for(int i=crossoverPoint;i<parent1Boolean.size();i++){
                child1Boolean.add(parent2Boolean.get(i));
                child2Boolean.add(parent1Boolean.get(i));
            }
        }else if(crossoverType==CrossoverType.MULTIPOINT){
            int flag=0, count=0;
            for(int i=0;i<parent1Boolean.size();i++){
                if(flag==0){
                    child1Boolean.add(parent1Boolean.get(i));
                    child2Boolean.add(parent2Boolean.get(i));
                }else{
                    child1Boolean.add(parent2Boolean.get(i));
                    child2Boolean.add(parent1Boolean.get(i));
                }
                count++;
                if(count==10){
                    count=0;
                    flag=(flag+1)%2;
                }             
            }
        }else{//UNIFORM
            for(int i=0;i<parent1Boolean.size();i++){
                if(ThreadLocalRandom.current().nextBoolean()){
                    child1Boolean.add(parent1Boolean.get(i));
                    child2Boolean.add(parent2Boolean.get(i));
                }else{
                    child1Boolean.add(parent2Boolean.get(i));
                    child2Boolean.add(parent1Boolean.get(i));
                }
            }
        }
        
        /**
         * Mutation
         */ 
        if(Parameters.ga_mutation){
            for(int i=0;i<child1Boolean.size();i++){
                if(ThreadLocalRandom.current().nextDouble()<Parameters.ga_mutationRate){
                    child1Boolean.set(i, !child1Boolean.get(i));
                }
                if(ThreadLocalRandom.current().nextDouble()<Parameters.ga_mutationRate){
                    child2Boolean.set(i, !child2Boolean.get(i));
                }
            }
        }
        
        GAIndividual[] offsprings = new GAIndividual[2];
        offsprings[0] = new GAIndividual(child1Boolean);
        offsprings[1] = new GAIndividual(child2Boolean);
        
        return offsprings;
    }
    
}