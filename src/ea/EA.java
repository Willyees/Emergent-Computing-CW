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
	
	private ArrayList<ArrayList<Individual>> populations = new ArrayList<ArrayList<Individual>>();
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
		for(int i = 0; i < Parameters.islandsN; i++)
			initialisePopulation();
		iteration = 0;
			
		while(iteration < Parameters.maxIterations * Parameters.islandsN){
			for(int island = 0; island < Parameters.islandsN; island++) {
				population = populations.get(island);
				for(int i = 0 ; i < Parameters.swapIterations; i++) {
					iteration++;
					//Individual parent1 = tournamentSelection();
					//Individual parent2 = tournamentSelection();
					Individual parent1 = tournamentSelection();
					Individual parent2 = tournamentSelection();
					ArrayList<Individual> parents = new ArrayList<Individual>();
					parents.add(parent1);parents.add(parent2);
					Individual child = crossoverArithmetic(parents);
					//Individual child = crossoverUniform(parent1, parent2);
					
					child = mutateGaussian(child, 100.0);
					child.evaluate(teamPursuit);
	
					replaceWorst(child);
					//printNumNoFinished();
					//Individual best = getBest(population);
					//best.print();
					//printStatsPopulation();
					//printStdDevPop();
				}	
			}
			MoveIndividualsNextIslandRandom(populations);
				
			//}
			
		}			
			Individual best = getBestIslands(populations);
			best.print();
			printBestStatsEnergy(best);
//		
//		System.out.println("Success rounds: " + successrounds);
	}


	/**
	 * will move individuals to the next island according to the parameters set in the Parameters class.
	 * Individuals are moved 1 by 1 for each island. So for each individual moved the population of the island in current focus will increase of 1 (moved individual can be selected again to be moved into next island)
	 * Movin them randomly
	 * @param islands
	 */
	private void MoveIndividualsNextIslandRandom(ArrayList<ArrayList<Individual>> islands) {
		
		//first island handled outside the loop
		//ArrayList<Individual> swap = new ArrayList<Individual>(Parameters.swapIndividuals);
 		for(int i = 0; i < Parameters.swapIndividuals; i++) {
			//add individuals from last island
			int lastIslandIndex = Parameters.rnd.nextInt(islands.get(islands.size() - 1).size());//get random index from last island
			islands.get(0).add(islands.get(islands.size() - 1).get(lastIslandIndex).copy());
			islands.get(islands.size() - 1).remove(lastIslandIndex);
		}
 		System.out.println("island0 pop:" + islands.get(0).size());
 		
		//other islands
		for(int island_index = 1; island_index < Parameters.islandsN; island_index++) {
			//swap by random individuals
			for(int i = 0; i < Parameters.swapIndividuals; i++) {
				int rndIndex = Parameters.rnd.nextInt(islands.get(island_index - 1).size()); //ramdom individual from previous island
				islands.get(island_index).add(islands.get(island_index - 1).get(rndIndex).copy()); //get previous random island individual
				islands.get(island_index - 1).remove(rndIndex);
			}
	 		System.out.println("island" + island_index + " pop:" + islands.get(island_index).size());

		}
		
		//DEBUG
		for(int i = 0; i < Parameters.islandsN; i++) {
			if(islands.get(i).size() != Parameters.popSize)
				throw new IllegalArgumentException("wrong size island");
		}
//		//last done outside loop
//		for(int i = 0; i < Parameters.swapIndividuals; i++) {
//			int rndIndex = Parameters.rnd.nextInt(islands.get(0).size());
//			swap[i] = rndIndex;
//			
//		}
		
	}
	
	//Debug functions
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
		System.out.println("meanFit: " + getMeanFitness() + "\t stddev: " + getStdDevFitness());
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
	
	//need index of the parents related to the population list.
//	private void replaceRandomParent(ArrayList<Individual> parents, Individual child) {
//		int indexReplace = Parameters.rnd.nextInt(parents.size());
//		population.remove(parents.get(indexReplace))
//	}
	
	private void replaceRandom(Individual child) {
		int indexReplace = Parameters.rnd.nextInt(population.size());
		population.set(indexReplace, child);
	}
	/**
	 * replace oldest and add + 1 to age of all population 
	 */
	private void replaceTournament(Individual child, int tournmentSize) {
		int[] indexes = new int[tournmentSize];
		for(int i = 0; i < tournmentSize; i++)
			indexes[i] = Parameters.rnd.nextInt(population.size());
		int indexReplace = getWorstIndex(population, indexes);
		
		population.set(indexReplace, child);
	}
	
	private void replaceOldest(Individual child) {
//		if(population.size() < Parameters.popSize * 2) {
//			population.add(child);
//			return;
//		}
		int index_oldest = 0;
		int age_oldest = -1;
		for(int i = 0; i < population.size(); i++){
			if(population.get(i).getAge() > age_oldest) {
				age_oldest = population.get(i).getAge();
				index_oldest = i;
			}
		}
		System.out.println("replaced oldest: " + index_oldest + ". Age is: " + population.get(index_oldest).getAge());
		population.set(index_oldest, child);
		
		//add +1 to age
		for(Individual i : population) 
			i.increaseAge();
	}
	
	private void replaceWorst(Individual child) {
		Individual worst = getWorst(population);
		if(child.getFitness() < worst.getFitness()){
			int idx = population.indexOf(worst);
			System.out.println("replaced worst: " + idx); 
			population.set(idx, child);
		}
	}

	//mutate with more chance on neighbour values (gaussian distribution). expecting mean
	private Individual mutateGaussian(Individual child, double stddev, double mean) {
		//pacing
		int mutationN = Parameters.rnd.nextInt(Parameters.mutationRateMax) + 1;
		for(int i = 0; i < mutationN; i++) {
			int index = Parameters.rnd.nextInt();
			int pacing = (int) (Parameters.rnd.nextGaussian() * stddev + mean);
			if(pacing > Parameters.WOMENS_PACING_STRATEGY_RANGE_MUTATION[0] && pacing < Parameters.DEFAULT_WOMENS_PACING_STRATEGY[1])
				child.pacingStrategy[index] = pacing;
		}
		
		//transition
		//mutate the transition strategy by flipping boolean value
		for(int i = 0; i < mutationN; i++){
			int index = Parameters.rnd.nextInt(child.transitionStrategy.length);
			child.transitionStrategy[index] = !child.transitionStrategy[index];
		}
		
		return child;
	}
	
	//overloaded function. not using any mean, instead using the current value and searching around it.
	private Individual mutateGaussian(Individual child, double stddev) {
		//pacing
		int mutationN = Parameters.rnd.nextInt(Parameters.mutationRateMax) + 1;
		for(int i = 0; i < mutationN; i++) {
			int mean = child.pacingStrategy[i];
			int index = Parameters.rnd.nextInt(child.pacingStrategy.length);
			int pacing = (int) (Parameters.rnd.nextGaussian() * stddev + mean);
			if(pacing > Parameters.WOMENS_PACING_STRATEGY_RANGE_MUTATION[0] && pacing < Parameters.DEFAULT_WOMENS_PACING_STRATEGY[1])
				child.pacingStrategy[index] = pacing;
		}
		
		//transition
		//mutate the transition strategy by flipping boolean value
		for(int i = 0; i < mutationN; i++){
			int index = Parameters.rnd.nextInt(child.transitionStrategy.length);
			child.transitionStrategy[index] = !child.transitionStrategy[index];
		}
		
		return child;
	}
	
	
	private Individual mutate(Individual child) {
		if(Parameters.rnd.nextDouble() > Parameters.mutationProbability){
			return child;
		}
		
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

	private Individual crossoverUniform(ArrayList<Individual> parents) {
		Individual child = new Individual();
				
		//pacing
		for(int i = 0; i < parents.get(0).pacingStrategy.length; i++) {
			int rnd = Parameters.rnd.nextInt(parents.size());
			child.pacingStrategy[i] = parents.get(rnd).pacingStrategy[i];

		}
		
		//transition strategy
		for(int i = 0; i < parents.get(0).transitionStrategy.length; i++) {
			int rnd = Parameters.rnd.nextInt(parents.size());
			child.transitionStrategy[i] = parents.get(rnd).transitionStrategy[i];
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
	
	private Individual getBestIslands(ArrayList<ArrayList<Individual>> populations) {
		double bestFitness = Double.MAX_VALUE;
		Individual best = null;
		int idx = 0;
		for(ArrayList<Individual> aPopulation : populations){
			for(Individual individual : aPopulation) {
				if(individual.getFitness() < bestFitness || best == null){
					best = individual;
					bestFitness = best.getFitness();
					idx++;
				}
			}
		}
		//System.out.println("best idx: " + idx);
		return best;
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
	
	//requires parameters of indexes 
	private int getWorstIndex(ArrayList<Individual> aPopulation, int[] rangeIndexes) {
		double worstFitness = Double.MIN_VALUE;
		int idx = 0;
		for(int index : rangeIndexes) {
			if(aPopulation.get(index) != null && aPopulation.get(index).getFitness() > worstFitness) {
				worstFitness = aPopulation.get(index).getFitness();
				idx = index;
			}
		}
		
		//System.out.println("best idx: " + idx);
		return idx;
	}
	
	private int getWorstIndex(ArrayList<Individual> aPopulation) {
		double worstFitness = Double.MIN_VALUE;
		int idx = 0;
		for(int index = 0; index < aPopulation.size(); index++){
			if(aPopulation.get(index) != null && aPopulation.get(index).getFitness() > worstFitness) {
				worstFitness = aPopulation.get(index).getFitness();
				idx = index;
			}
		}
		//System.out.println("best idx: " + idx);
		return idx;
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
		ArrayList<Individual> population = new ArrayList<Individual>();
		while(population.size() < Parameters.popSize){
			Individual individual = new Individual();
			individual.initialise();		
			individual.evaluate(teamPursuit);
			population.add(individual);
		}		
		populations.add(population);
	}	
	/**
	 * 
	 * @param randomP: 0-100 integer to specify the random part of parents from old population to use for new population
	 */
	private void replacePopulation(int randomP) {
		//create pool parents. 60 random and 40 best parents
		int popSize = population.size();
		int randomN = (int) (popSize / 100.0 * randomP);
		int bestN = popSize - randomN;
		ArrayList<Individual> parents = new ArrayList<Individual>(popSize);
		//random part
		for(int i = 0; i < randomN; i++) {
			parents.add(population.get(Parameters.rnd.nextInt(popSize)).copy());
		}
		//best ones
		for(int i = 0; i < bestN; i++) {
			int indexBest = getBestIndex(population);
			parents.add(population.get(indexBest).copy());
			population.remove(indexBest); //remove best index so second best can be found on next iteration
		}
		//
		
		
		//create children
		ArrayList<Individual> children = new ArrayList<Individual>(popSize);
		for(int i = 0; i < popSize; i++) {
			int[] indexParents = new int[Parameters.parentsN];
			for(int j = 0; j < Parameters.parentsN; j++) {
				indexParents[j] = Parameters.rnd.nextInt(parents.size());
			}
			Individual child = crossoverUniform(parents);
			child = mutate(child);
			child.evaluate(teamPursuit);
			children.add(child);
		}
	
		population.clear();
		for(Individual c : children)
			population.add(c.copy());
		System.out.println("Replaced old population with new one, " + randomP + " % parents");
		//replace original population
	}
}
