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
		boolean[] transitionStrategy = {true,true,true,true,true,true,true,true,false,false,true,true,false,false,true,false,false,false,false,false,false,true};
		int[] pacingStrategy = {382,615,368,661,433,428,322,719,744,246,289,522,274,217,800,462,200,385,386,346,296,542,666};
		
		try {
			SimulationResult result = teamPursuit.simulate(transitionStrategy, pacingStrategy);
			System.out.println("evaluated: " + result.getFinishTime());
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
			int counter = 0;
			int convergeStdDev = 0;
			boolean notfound = true;
			while(iteration < Parameters.maxIterations){
				
				iteration++;
				Individual parent1 = tournamentSelection();
				Individual parent2 = tournamentSelection();
				Individual child = crossover(parent1, parent2);
				child = mutate(child);
				child.evaluate(teamPursuit);
				
				replace(child);
				printNumNoFinished();
//				Individual best = getBest(population);
//				best.print();
				printStatsPopulation();
//				if(getStdDevFitness() < 5.0 && notfound) {
//					notfound = false;
//					convergeStdDev = counter;
//				}
//				counter++;
			//}
			System.out.println("converge stddev: " + convergeStdDev);
			Individual best = getBest(population);
			best.print();
			iteration = 0;
//			if(getMeanFitness() != 1000.0)
//				successrounds++;
		}						
//		
//		System.out.println("Success rounds: " + successrounds);
	}

	//Debug functions
	
	
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
	
	private void replace(Individual child) {
		Individual worst = getWorst(population);
		if(child.getFitness() < worst.getFitness()){
			int idx = population.indexOf(worst);
			System.out.println("replaced worst: " + idx); 
			population.set(idx, child);
		}
	}


	private Individual mutate(Individual child) {
		if(Parameters.rnd.nextDouble() > Parameters.mutationProbability){
			return child;
		}
		
		// choose how many elements to alter
		int mutationRate = 1 + Parameters.rnd.nextInt(Parameters.mutationRateMax);
		
		// mutate the transition strategy

		//mutate the transition strategy by flipping boolean value
		for(int i = 0; i < mutationRate; i++){
			int index = Parameters.rnd.nextInt(child.transitionStrategy.length);
			child.transitionStrategy[index] = !child.transitionStrategy[index];
		}
			
		//mutate the pacing strategy
		for(int i = 0; i < mutationRate; i++) {
			int index = Parameters.rnd.nextInt(child.pacingStrategy.length);
			child.pacingStrategy[index] = Parameters.rnd.nextInt(Parameters.WOMENS_PACING_STRATEGY_RANGE[1] - Parameters.WOMENS_PACING_STRATEGY_RANGE[0] + 1)  + Parameters.WOMENS_PACING_STRATEGY_RANGE[0];
		}
		
		return child;
	}


	private Individual crossover(Individual parent1, Individual parent2) {
		if(Parameters.rnd.nextDouble() > Parameters.crossoverProbability){
			return parent1;
		}
		Individual child = new Individual();
		
		int crossoverPoint = Parameters.rnd.nextInt(parent1.transitionStrategy.length);
		
		// just copy the pacing strategy from p1 - not evolving in this version
		for(int i = 0; i < parent1.pacingStrategy.length; i++){			
			child.pacingStrategy[i] = parent1.pacingStrategy[i];
		}
		
		
		for(int i = 0; i < crossoverPoint; i++){
			child.transitionStrategy[i] = parent1.transitionStrategy[i];
		}
		for(int i = crossoverPoint; i < parent2.transitionStrategy.length; i++){
			child.transitionStrategy[i] = parent2.transitionStrategy[i];
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
		//System.out.println("best idx: " + idx);
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
