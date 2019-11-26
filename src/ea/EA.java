package ea;

/***
 * This is an example of an EA used to solve the problem
 *  A chromosome consists of two arrays - the pacing strategy and the transition strategy
 * This algorithm is only provided as an example of how to use the code and is very simple - it ONLY evolves the transition strategy and simply sticks with the default
 * pacing strategy
 * The default settings in the parameters file make the EA work like a hillclimber:
 * 	the population size is set to 1, and there is no crossover, just mutation
 * The pacing strategy array is never altered in this version- mutation and crossover are only
 * applied to the transition strategy array
 * It uses a simple (and not very helpful) fitness function - if a strategy results in an
 * incomplete race, the fitness is set to 1000, regardless of how much of the race is completed
 * If the race is completed, the fitness is equal to the time taken
 * The idea is to minimise the fitness value
 */

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

import teamPursuit.SimulationResult;
import teamPursuit.TeamPursuit;
import teamPursuit.WomensTeamPursuit;

public class EA implements Runnable{
	
	// create a new team with the default settings
	public static TeamPursuit teamPursuit = new WomensTeamPursuit(); 
	
	private ArrayList<Individual> population = new ArrayList<Individual>();
	private int iteration = 0;
	
	public EA() {
		
	}

	
	public static void main(String[] args) {
		EA ea = new EA();
		//ea.testSpecificChromosome();
		ea.run();
	}

	public void testSpecificChromosome() {
		boolean[] transitionStrategy = {true,false,true,false,true,true,true,false,false,false,true,true,true,true,false,true,false,false,false,false,false,false};
		int[] pacingStrategy = {424,372,257,470,251,475,416,480,421,434,414,410,474,483,457,441,316,475,411,455,488,289,388};
		
		try {
			SimulationResult result = teamPursuit.simulate(transitionStrategy, pacingStrategy);
			System.out.println("evaluated: " + result.getFinishTime());
			System.out.println("proportion: " + result.getProportionCompleted());
			double[] velocities = result.getVelocityProfile();
			System.out.print("velocities: ");
			for(double velocity : velocities) {
				System.out.print(velocity + "  ");
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	public void run() {
		//int successrounds = 0;
		Parameters.maxIterations = 500;
		//for(int outerIteration = 0; outerIteration < 25; outerIteration++) {
			initialisePopulation();	
			iteration = 0;
			while(iteration < Parameters.maxIterations){
				
				iteration++;
				//Individual parent1 = tournamentSelection();
				//Individual parent2 = tournamentSelection();
				Individual parent1 = selectionRanking();
				Individual parent2 = selectionRanking();
				Individual parent3 = selectionRanking();
				ArrayList<Individual> parents = new ArrayList<Individual>();
				parents.add(parent1);parents.add(parent2);parents.add(parent3);
				Individual child = crossoverArithmetic(parents);
				//Individual child = crossoverUniform(parent1, parent2);
				child = mutate(child);
				child.evaluate(teamPursuit);
				
				replace(child);
				//printNumNoFinished();
				Individual best = getBest(population);
				best.print();
				printStatsPopulation();
				//printStdDevPop();
//				
			//}
			
			
			//iteration = 0;
//			
		}			
			Individual best = getBest(population);
			best.print();
			printBestStatsEnergy(best);
//		
//		System.out.println("Success rounds: " + successrounds);
	}

	//Debug functions
	
	private void printStdDevPop() {
		double stddev = getStdDevPop();
		System.out.println("stddev pop pacing: " + stddev);
	}
	
	private double getStdDevPop() {
		double stddev = 0.0;
		double mean = 0.0;
		for(Individual i : population) {
			for(int index = 0; index < i.pacingStrategy.length; index++) {
				mean += i.pacingStrategy[index];
			}
		}
		mean /= (population.get(0).pacingStrategy.length * population.size());
		
		for(Individual i : population) {
			double mean_i = 0.0; //mean for individuals
			double stddev_i = 0.0;
			for(int index = 0; index < i.pacingStrategy.length; index++) {
				mean_i += i.pacingStrategy[index];
				
			}
			mean_i /= i.pacingStrategy.length;
			stddev += Math.pow((mean_i - mean), 2);
			
		}
		stddev /= population.size();
		stddev = Math.sqrt(stddev);
		
		return stddev;
	}
	
	private void printBestStats() {
		Individual best = getBest(population);
		System.out.println("best stats: " + "time: " + best.result.getFinishTime() + "\t completed: " + best.result.getProportionCompleted() + "\t");
		printBestStatsEnergy(best);
	}
	
	/*
	 * Called from printBestStats. Have to pass best individual since is already calculated in the calling function
	 */
	private void printBestStatsEnergy(Individual best) {
		for(double energy : best.result.getEnergyRemaining())
			System.out.print("remaining energy: " + energy + " ");
	}
	
	private int getLowFitnessNumPopulation() {
		int lowFitnessCounter = 0;
		for(Individual i : population) {
			if(i.getFitness() > 999)
				lowFitnessCounter++;
		}
		return lowFitnessCounter;
		
	}
	private void printStats() {		
		System.out.println("" + iteration + "\t" + getBest(population) + "\t" + getWorst(population));	
		printBestStats();
		printNumNoFinished();
	}

	/**
	 * print for the whole population different stats (proportion completed and energy remaining)
	 */
	private void printStatsPopulation() {
//		for(Individual i : population) {
//			System.out.println("completed%: " + i.result.getProportionCompleted() + "\t energy: " + i.result.getEnergyRemaining());
//		}
		System.out.println("lowfitnessNum: " + getLowFitnessNumPopulation() + "\t meanFit: " + getMeanFitness() + "\t stddev: " + getStdDevFitness());
	}
	/**
	 * Print the proportion completed for the whole population
	 */
	private void printFinishIndividuals() {
		for(Individual i : population) {
			System.out.println(i.result.getProportionCompleted() + " ");
		}
	}
	private void printNumNoFinished() {
		int numberDidntFinish = 0;
		for(Individual i : population) {
			if(i.result.getProportionCompleted() < 0.99)
				numberDidntFinish++;
		}
		System.out.println("didnt finish: " + numberDidntFinish);
		
		
	}
	
	
	
	private double getMeanFitness() {
		double mean = 0.0;
		for(Individual i : population) {
			mean += i.getFitness();
		}
		mean /= population.size();
		return mean;
	}
	
	private double getStdDevFitness() {
		double stddev = 0.0;
		double mean = getMeanFitness();
		for(Individual i : population) {
			stddev += Math.pow((i.getFitness() - mean), 2);
		}
		stddev /= population.size();
		stddev = Math.sqrt(stddev);
		return stddev;
	}
	
	//End Debug functions
	
	//get absolute fitness of the population
	private double getAbsoluteFitness() {
		double totalFitness = 0.0;
		for(Individual i : population) {
			totalFitness += i.getFitness();
		}
		return totalFitness;		
	}
	
	private void replace(Individual child) {
		Individual worst = getWorst(population);
		if(child.getFitness() < worst.getFitness()){
			int idx = population.indexOf(worst);
			//System.out.println("replaced worst: " + idx); 
			population.set(idx, child);
		}
	}


	private Individual mutate(Individual child) {
		if(Parameters.rnd.nextDouble() > Parameters.mutationProbability){
			return child;
		}
		
		// choose how many elements to alter
		int mutationRate = Parameters.rnd.nextInt(Parameters.mutationRateMax) + 1;
		
		// mutate the transition strategy

		//mutate the transition strategy by flipping boolean value
		for(int i = 0; i < mutationRate; i++){
			int index = Parameters.rnd.nextInt(child.transitionStrategy.length);
			child.transitionStrategy[index] = !child.transitionStrategy[index];
		}
			
		//mutate the pacing strategy
		for(int i = 0; i < mutationRate; i++) {
			int index = Parameters.rnd.nextInt(child.pacingStrategy.length);
			child.pacingStrategy[index] = Parameters.rnd.nextInt(Parameters.WOMENS_PACING_STRATEGY_RANGE_MUTATION[1] - Parameters.WOMENS_PACING_STRATEGY_RANGE_MUTATION[0] + 1)  + Parameters.WOMENS_PACING_STRATEGY_RANGE_MUTATION[0];
		}
		
		return child;
	}

	private Individual crossoverUniform(Individual parent1, Individual parent2) { //have to use array of parents in case multiple parents are used
		Individual child = new Individual();
				
		//pacing
		for(int i = 0; i < parent1.pacingStrategy.length; i++) {
			int rnd = Parameters.rnd.nextInt(Parameters.parentsN);
			if(rnd == 0) {
				child.pacingStrategy[i] = parent1.pacingStrategy[i];
			}else if(rnd == 1) {
				child.pacingStrategy[i] = parent2.pacingStrategy[i];
			}
		}
		
		//transition strategy
		for(int i = 0; i < parent1.transitionStrategy.length; i++) {
			int rnd = Parameters.rnd.nextInt(Parameters.parentsN);
			if(rnd == 0) {
				child.transitionStrategy[i] = parent1.transitionStrategy[i];
			}else if(rnd == 1) {
				child.transitionStrategy[i] = parent2.transitionStrategy[i];
			}
		}
		return child;
	}
	
	private Individual crossover(Individual parent1, Individual parent2) {
		if(Parameters.rnd.nextDouble() > Parameters.crossoverProbability){
			return parent1;
		}
		Individual child = new Individual();
				
		//pacing strategy
		int crossoverPointPacing = Parameters.rnd.nextInt(parent1.pacingStrategy.length);
		for(int i = 0; i < crossoverPointPacing; i++){
			child.pacingStrategy[i] = parent1.pacingStrategy[i];
		}
		for(int i = crossoverPointPacing; i < parent2.pacingStrategy.length; i++){
			child.pacingStrategy[i] = parent2.pacingStrategy[i];
		}
		
		//transition strategy
		int crossoverPointTransition = Parameters.rnd.nextInt(parent1.transitionStrategy.length);
		for(int i = 0; i < crossoverPointTransition; i++){
			child.transitionStrategy[i] = parent1.transitionStrategy[i];
		}
		for(int i = crossoverPointTransition; i < parent2.transitionStrategy.length; i++){
			child.transitionStrategy[i] = parent2.transitionStrategy[i];
		}
		return child;
	}


	private Individual crossoverNPoints(ArrayList<Individual> parents, int nPoints) {
		if(Parameters.rnd.nextDouble() > Parameters.crossoverProbability){
			return parents.get(0);
		}
		Individual child = new Individual();
		//pacing
		int pacingLength = parents.get(0).pacingStrategy.length;
		ArrayList<Integer> points = new ArrayList<Integer>(nPoints);
		for(int i = 0; i < nPoints; i++)
			points.add(Parameters.rnd.nextInt(pacingLength));
		Collections.sort(points);
		//first cut and last cut done outside the main loop
		int parent_id = 0;
		for(int i = 0; i < points.get(0); i++)
			child.pacingStrategy[i] = parents.get(parent_id).pacingStrategy[i];
		parent_id++;
		
		for(int index = 0; index < points.size() - 1; index++) {
			for(int i = points.get(index); i < points.get(index + 1); i++) {
				child.pacingStrategy[i] = parents.get(parent_id).pacingStrategy[i];
			}
			parent_id++;
			if(parent_id % parents.size() == 0)
				parent_id = 0;
		}
		//last index
		for(int i = points.get(points.size() - 1); i < pacingLength; i++){
			child.pacingStrategy[i] = parents.get(parent_id).pacingStrategy[i];
		}
		
		//transition
		int transitionLength = parents.get(0).transitionStrategy.length;
		points.clear();//might be better using same points because transition and pacing might be optimized one with the other
		for(int i = 0; i < nPoints; i++)
			points.add(Parameters.rnd.nextInt(transitionLength));
		Collections.sort(points);
		parent_id = 0;
		for(int i = 0; i < points.get(0); i++)
			child.transitionStrategy[i] = parents.get(parent_id).transitionStrategy[i];
		parent_id++;
		
		for(int index = 0; index < points.size() - 1; index++) {
			for(int i = points.get(index); i < points.get(index + 1); i++) {
				child.transitionStrategy[i] = parents.get(parent_id).transitionStrategy[i];
			}
			parent_id++;
			if(parent_id % parents.size() == 0)
				parent_id = 0;
		}
		
		//last index
		for(int i = points.get(points.size() - 1); i < transitionLength; i++){
			child.transitionStrategy[i] = parents.get(parent_id).transitionStrategy[i];
		}
		
		return child;
	}
	
	//outputting a child having its pacing strategy as a mean of its parents. transition strategy is uniform crossover
	private Individual crossoverArithmetic(ArrayList<Individual> parents) {
		Individual child = new Individual();
		//transition strategy
		int count_true = 0;
		for(int i = 0; i < parents.get(0).transitionStrategy.length; i++) {
			for(int index = 0; index < parents.size(); index++) {
				if(parents.get(index).transitionStrategy[i] == true)
					count_true++; 
			}
			if(count_true == parents.size() - count_true)
				child.transitionStrategy[i] = Parameters.rnd.nextBoolean();
			else
				child.transitionStrategy[i] = (count_true >  (parents.size() - count_true)) ? true : false;
			count_true = 0;
		}
		
		//pacing strategy
		double average = 0.0;
		for(int i = 0; i < parents.get(0).pacingStrategy.length; i++) {
			for(int index = 0; index < parents.size(); index++)
				average += parents.get(index).pacingStrategy[i];
			average /= parents.size();
			child.pacingStrategy[i] = (int) average;
			average = 0.0;
		}
		return child;
	}
	/**
	 * Returns a COPY of the individual selected using tournament selection
	 * @return
	 */
	private Individual tournamentSelection() {
		ArrayList<Individual> candidates = new ArrayList<Individual>();
		for(int i = 0; i < Parameters.tournamentSize; i++){
			candidates.add(population.get(Parameters.rnd.nextInt(population.size())));
		}
		return getBest(candidates).copy();
	}
	
	private int[] getRankingByFitness(ArrayList<Individual> aPopulation) { //very slow
		//rank fitnesses
		int[] ranks = new int[aPopulation.size()];
		ArrayList<Individual> pop = new ArrayList<Individual>(aPopulation.size());
		//create copy of list that will be modified
		for(int index = 0; index < aPopulation.size(); index++) {
			pop.add(aPopulation.get(index).copy());
		}
				
		for(int index = aPopulation.size() - 1; index > 0; index--) {
			int bestIndex = getBestIndex(pop);
			ranks[bestIndex] = index;
			pop.set(bestIndex, null);
		}
		
//		for(int i = 0; i < aPopulation.size(); i++) { //in case two 
//			int count = 0;
//			for(int j = 0; j < aPopulation.size(); j++) {
//				if(aPopulation.get(i).getFitness() > aPopulation.get(i).getFitness())
//					count++;
//			}
//			ranks[i] = count;
//		}
		return ranks;
	}
	
	private Individual selectionRanking() {
		int[] ranks = getRankingByFitness(population);
		double[] probabilities = new double[population.size()];
		for(int i = 0; i < population.size(); i ++)		
			probabilities[i] = (2 - Parameters.scalingFactorRankingSelection) / population.size() + 
			(2 * ranks[i] * (Parameters.scalingFactorRankingSelection - 1) / (population.size() * (population.size() - 1)));
		
		double rnd = Parameters.rnd.nextDouble();
		
		double partial_sum = 0;
    	int parent_id = 0;
 
    	while(rnd > partial_sum) {
    		partial_sum += probabilities[parent_id];
    		parent_id += 1;
    	}
    	parent_id -= 1; //have to get lower parent
    	Individual selected = population.get(parent_id).copy();
    	
    	return selected;
	}
	
	private Individual selectionFitnessProportionate() {
		double absoluteFitness = getAbsoluteFitness();
		
		double partial_sum = 0;
    	int parent_id = 0;
    	double random_fitness = Parameters.rnd.nextDouble() * absoluteFitness;
    	    	
    	while(random_fitness > partial_sum) {
    		partial_sum += population.get(parent_id).getFitness();
    		parent_id += 1;
    	}
    	parent_id -= 1; //have to get lower parent
    	Individual selected = population.get(parent_id).copy();
    	return selected;
	}

	private int getBestIndex(ArrayList<Individual> aPopulation) {
		double bestFitness = Double.MAX_VALUE;
		int idx = 0;
		for(int index = 0; index < aPopulation.size(); index++){
			if(aPopulation.get(index) != null && aPopulation.get(index).getFitness() < bestFitness) {
				bestFitness = aPopulation.get(index).getFitness();
				idx = index;
			}
		}
		//System.out.println("best idx: " + idx);
		return idx;
	}
	
	private Individual getBest(ArrayList<Individual> aPopulation) {
		double bestFitness = Double.MAX_VALUE;
		Individual best = null;
		int idx = 0;
		for(Individual individual : aPopulation){
			if(individual.getFitness() < bestFitness || best == null){
				best = individual;
				bestFitness = best.getFitness();
				idx++;
			}
		}
		System.out.println("best idx: " + idx);
		return best;
	}

	private Individual getWorst(ArrayList<Individual> aPopulation) {
		double worstFitness = 0;
		Individual worst = null;
		int idx = 0;
		for(Individual individual : population){
			if(individual.getFitness() > worstFitness || worst == null){
				worst = individual;
				worstFitness = worst.getFitness();
				idx++;
			}
		}
		//System.out.println("worst: " + idx);
		return worst;
	}
	
	private void printPopulation() {
		for(Individual individual : population){
			System.out.println(individual);
		}
	}

	private void initialisePopulation() {
		System.out.println("cleared population");
		population.clear();
		while(population.size() < Parameters.popSize){
			Individual individual = new Individual();
			individual.initialise();			
			individual.evaluate(teamPursuit);
			population.add(individual);
							
		}		
	}	
}
