package ea;

import java.util.Random;

public class Parameters {

	public static Random rnd = new Random(System.currentTimeMillis());
	
	/**
	 * used as a seed
	 * 
	 */	
	static final boolean [] DEFAULT_WOMENS_TRANSITION_STRATEGY = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
	public static final int [] DEFAULT_WOMENS_PACING_STRATEGY = {300, 300, 300, 300, 300, 300, 300, 350, 350, 300, 300, 350, 350, 350, 350, 300, 300, 350, 350, 350, 350, 300, 300};
	static final int WOMENS_PACING_STRATEGY_RANGE [] = {200,600}; //min: 200; max: 1200
	static final int WOMENS_PACING_STRATEGY_RANGE_MUTATION [] = {200,1200}; //min: 200; max: 1200
	static final int WOMENS_TRANSITION_N = 21; 
	public static int popSize = 1;
	public static int tournamentSize = 2;
	
	public static int mutationRateMax = 2;//out of len
	//public static double mutationProbability = 0.5; always mutating (move operator)
	public static double crossoverProbability = 1.0;
	
	public static int maxIterations = 100;
	public static int applyMoveOperator = 6;
	public static int swapOperator = 3;
	
	public final static double coolingRate = 0.98;
	public final static double startTemp = 400.0;
	public final static double finalTemp = 0.0;
	public static double currentTemp = startTemp;
	//DEBUG
	public static double lowFitness = 1000.0;
	
	
}
