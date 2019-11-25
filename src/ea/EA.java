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

//import com.sun.org.apache.xalan.internal.xsltc.runtime.Parameter;

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
		Parameters.maxIterations = 12000;
		//for(int outerIteration = 0; outerIteration < 25; outerIteration++) {
			initialisePopulation();	
			double bestOverallTime = Double.POSITIVE_INFINITY;
			Individual bestOverall = new Individual();
			iteration = 0;
			int counterPreviousChrom = 0;
			int counterInner = 0;
			ArrayList<Individual> modified = new ArrayList<Individual>(Parameters.applyMoveOperator);
			int index = 0;
			while(iteration < Parameters.maxIterations && Parameters.currentTemp > Parameters.finalTemp){
				modified.clear();
				iteration++;
				
				//int index = getRandomIndividualIndex(); //not using because only 1 chormosome in the population
				
				Individual i1 = population.get(index).copy();
				double previousFitness = population.get(index).getFitness();
				moveOperatorMutate(i1);
				i1.evaluate(teamPursuit);
				replace(i1, index);
				
				if(previousFitness == population.get(index).getFitness()) { //most likely same chromosome as before since the fitness didnt change 
					counterPreviousChrom++;
				}
				else
					counterPreviousChrom = 0;
				if(counterPreviousChrom == 400) { //same chromosome for 100 iterations
					Parameters.currentTemp = Parameters.reheatTemp * 0.6;
					Parameters.reheatTemp = Parameters.currentTemp;
					counterInner++;
					counterPreviousChrom = 0;
				}
				if(population.get(index).getFitness() < bestOverallTime ) {
					bestOverall = population.get(index).copy();
					bestOverallTime = bestOverall.getFitness();
				}
//				population.get(index).evaluate(teamPursuit);
//				for(int i = 0; i < Parameters.applyMoveOperator; i++) {
//					Parameters.mutationRateMax = Parameters.WOMENS_TRANSITION_N / Parameters.applyMoveOperator * (i + 1);//setting mutation rate 4-8-12-16-20
//					Individual i1 = population.get(index).copy();
//					moveOperatorMutate(i1);
//					i1.evaluate(teamPursuit);
//					modified.add(i1);
//				}
//				
//				for(Individual i : modified) {
//					printStatsIndividual(i);
//					replace(i, index);
//				}
				
				
				//Individual best = getBest(population);
				//best.print();
				printStatsPopulation();
		
		}			
			Individual best = getBest(population);
			best.print();
			printBestStatsEnergy(best);		
			System.out.println("times temp up: " + counterInner);
			printStatsIndividual(bestOverall);
//		System.out.println("Success rounds: " + successrounds);
	}

	//Debug functions
	
	private void printStatsIndividual(Individual i) {
		System.out.println("individual: " + "time: " + i.result.getFinishTime() + "\t completed: " + i.result.getProportionCompleted() + "\t");
	}
	
	private void printPacingPop() {
		for(Individual i : population) {
			for(int index = 0; index < i.pacingStrategy.length; index++)
				System.out.print(i.pacingStrategy[index] + " ");
			System.out.print("\n");
		}
		
	}
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
		System.out.println("meanFit: " + getMeanFitness()); //+ "\t stddev: " + getStdDevFitness());
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
	
	private void replace(Individual child, int index) {
		
		if(child.getFitness() < population.get(index).getFitness()){
			System.out.println("replaced"); 
			population.set(index, child);
		}
		else {
			//accept worse solution with some probability
			double diffFitness = child.getFitness() - population.get(index).getFitness();
			double rnd = Parameters.rnd.nextDouble();
			double prob = Math.exp(- Math.abs(diffFitness) / Parameters.currentTemp);
			if(prob > rnd) {
				System.out.println("up");
				population.set(index, child);
			}
		}
		
		Parameters.currentTemp *= Parameters.coolingRate; //reduce geometrically
//		Parameters.currentTemp = Parameters.currentTemp / (1 + (1 - Parameters.coolingRate) * Parameters.currentTemp); //reduce Lundy&Mees 	
		System.out.println("temp: " + Parameters.currentTemp);
	}

	private Individual getRandomIndividual() {
		int index = Parameters.rnd.nextInt(population.size());
		return population.get(index);
	}
	
	private int getRandomIndividualIndex() {
		return Parameters.rnd.nextInt(population.size());
	}

	private Individual mutate(Individual child) {
		
		// choose how many elements to alter
		int mutationN = Parameters.rnd.nextInt(Parameters.mutationRateMax) + 1;
		
		// mutate the transition strategy

		//mutate the transition strategy by flipping boolean value
		for(int i = 0; i < mutationN; i++){
			int index = Parameters.rnd.nextInt(child.transitionStrategy.length);
			child.transitionStrategy[index] = !child.transitionStrategy[index];
		}
			
		//mutate the pacing strategy
		for(int i = 0; i < mutationN; i++) {
			int index = Parameters.rnd.nextInt(child.pacingStrategy.length);
			child.pacingStrategy[index] = Parameters.rnd.nextInt(Parameters.WOMENS_PACING_STRATEGY_RANGE_MUTATION[1] - Parameters.WOMENS_PACING_STRATEGY_RANGE_MUTATION[0] + 1)  + Parameters.WOMENS_PACING_STRATEGY_RANGE_MUTATION[0];
		}
		
		return child;
	}

	/**
	 *move operator. at the moment using the usual mutate operator
	 */
	private Individual moveOperatorMutate(Individual i) {
		return mutate(i);	
	}
	
	//not very good because there is no exploration for the pacing strategy
	private Individual moveOperatorSwap(Individual i) {
		//transition
		int[] indexSwapTransition = new int[Parameters.swapOperator];
		indexSwapTransition[0] = Parameters.rnd.nextInt(i.transitionStrategy.length);
		for(int index = 1; index < Parameters.swapOperator; index++) {
			indexSwapTransition[index] = Parameters.rnd.nextInt(i.transitionStrategy.length);
			boolean temp = i.transitionStrategy[index - 1];
			i.transitionStrategy[index - 1] = i.transitionStrategy[index];
			i.transitionStrategy[index] = temp;
		}
		
				
		//pacing
		int[] indexSwapPacing = new int[Parameters.swapOperator];
		indexSwapPacing[0] = Parameters.rnd.nextInt(i.pacingStrategy.length);
		for(int index = 1; index < Parameters.swapOperator; index++) {
			indexSwapPacing[index] = Parameters.rnd.nextInt(i.pacingStrategy.length);
			int temp = i.pacingStrategy[index - 1];
			i.pacingStrategy[index - 1] = i.pacingStrategy[index];
			i.pacingStrategy[index] = temp;
		}
		
		return i;
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
