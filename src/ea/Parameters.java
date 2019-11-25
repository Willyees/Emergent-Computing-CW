package ea;

import java.util.Random;

public class Parameters {

	public static Random rnd = new Random(System.currentTimeMillis());
	
	static final boolean [] DEFAULT_WOMENS_TRANSITION_STRATEGY = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
	public static final int [] DEFAULT_WOMENS_PACING_STRATEGY = {300, 300, 300, 300, 300, 300, 300, 350, 350, 300, 300, 350, 350, 350, 350, 300, 300, 350, 350, 350, 350, 300, 300};
	static final int WOMENS_PACING_STRATEGY_RANGE [] = {200,600}; //min: 200; max: 1200
	static final int WOMENS_PACING_STRATEGY_RANGE_MUTATION [] = {200,1200}; //min: 200; max: 1200
	
	public static int popSize = 10;
	public static int tournamentSize = 2;
	public static int parentsN = 2;
	
	public static double scalingFactorRankingSelection = 1.5;
	public static int mutationRateMax = 20;//out of len
	public static double mutationProbability = 0.6;
	public static double crossoverProbability = 1.0;
	
	public static int maxIterations = 100;
	
	//DEBUG
	public static double lowFitness = 1000.0;
	
	
}
