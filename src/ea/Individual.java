package ea;

import teamPursuit.*;

public class Individual {

	
	boolean[] transitionStrategy = new boolean[22] ;
	int[] pacingStrategy = new int[23];
	
	SimulationResult result = null;	
	
	public Individual() {		
		
	}

	// this code just evolves the transition strategy
	// an individual is initialised with a random strategy that will evolve
	// the pacing strategy is initialised to the default strategy and remains fixed
	
	public void initialise(){
		for(int i = 0; i < transitionStrategy.length; i++){
			transitionStrategy[i] = Parameters.rnd.nextBoolean();
		}
		
		for(int i = 0; i < pacingStrategy.length; i++){
			pacingStrategy[i] = Parameters.rnd.nextInt(Parameters.WOMENS_PACING_STRATEGY_RANGE[1] - Parameters.WOMENS_PACING_STRATEGY_RANGE[0] + 1)  + Parameters.WOMENS_PACING_STRATEGY_RANGE[0];
			
			if(pacingStrategy[i] > Parameters.WOMENS_PACING_STRATEGY_RANGE[1] || pacingStrategy[i] < Parameters.WOMENS_PACING_STRATEGY_RANGE[0])
				throw new IllegalArgumentException("outsiderange " + pacingStrategy[i]);
		}
		
		
	}
	
	// this is just there in case you want to check the default strategies
	public void initialise_default(){
		for(int i = 0; i < transitionStrategy.length; i++){
			transitionStrategy[i] = Parameters.DEFAULT_WOMENS_TRANSITION_STRATEGY[i];
		}
		
		for(int i = 0; i < pacingStrategy.length; i++){
			pacingStrategy[i] = Parameters.DEFAULT_WOMENS_PACING_STRATEGY[i];
		}
		
	}
	
	
	public void evaluate(TeamPursuit teamPursuit){		
		try {
			result = teamPursuit.simulate(transitionStrategy, pacingStrategy);
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	public double getFitness(){
//		double fitness = Parameters.lowFitness;		
//		if (result == null || result.getProportionCompleted() < 0.999){
//			double additional = fitness * (0.999 - result.getProportionCompleted());
//			fitness += additional; 
//		}
//		else{				
//			fitness = result.getFinishTime();
//			
//		}
		double total = 0.0;
		int segment_finished = -1;
		double[] velocities = result.getVelocityProfile();
		for(int race_segment = 0; race_segment < 23; race_segment++) {
			double velocity = velocities[race_segment];
			if(velocity == 0.0) {
				segment_finished = race_segment;
				break;
			}
			if(race_segment == 0 || race_segment == 22) { //if first or last
				total += (187.5 / velocity);
				
			}
			else {
				total += (Parameters.WOMENS_HALF_LAP_DISTANCE / velocity);
			}
		}
		if(segment_finished != -1) {
			//System.out.print("didnt finish");
			double average_spd = 0.0;
			for(int i = 0; i < segment_finished; i++) {
				average_spd += velocities[i];
				 //segment finished is the next segment that was not done
			}
			average_spd /= segment_finished;
			int segments_left = 23 - segment_finished; //how many segments are missing
			//DEBUG
			int count = 0;
			for(double v : velocities) {
				if(v == 0.0)
					count++;
			}
			if(count != segments_left)
				throw new IllegalArgumentException("segments left different from 0.0 velociteis");
			//END DEBUG
			total += (segments_left - 1) * Parameters.WOMENS_HALF_LAP_DISTANCE /  average_spd * 1.2 + (187.5 / average_spd * 1.2); 
				
			//System.out.print("fitness: " + total + "\n");	
		}
//		else
//			System.out.print("finsihed");
//		double[] energies = result.getEnergyRemaining();
//		for(double energy : energies) {
//			fitness += energy;
//		}
		
		return total;
	}
	
	
	
	public Individual copy(){
		Individual individual = new Individual();
		for(int i = 0; i < transitionStrategy.length; i++){
			individual.transitionStrategy[i] = transitionStrategy[i];
		}		
		for(int i = 0; i < pacingStrategy.length; i++){
			individual.pacingStrategy[i] = pacingStrategy[i];
		}		
		individual.evaluate(EA.teamPursuit);	
		return individual;
	}
	
	@Override
	public String toString() {
		String str = "";
		if(result != null){
			str += getFitness();
		}
		return str;
	}

	public void print() {
		for(int i : pacingStrategy){
			System.out.print(i + ",");			
		}
		System.out.println();
		for(boolean b : transitionStrategy){
			if(b){
				System.out.print("true,");
			}else{
				System.out.print("false,");
			}
		}
		System.out.println("\ntime best: " + result.getFinishTime());
		System.out.println("proportion completed: " + result.getProportionCompleted());
		System.out.println("\r\n" + this);
	}
}
