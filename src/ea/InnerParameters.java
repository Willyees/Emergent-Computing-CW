package ea;

public class InnerParameters {
	public int tournamentSize = 2;
	public int parentsN = 2;
	public int swapIndividuals = 5;
	
	public double scalingFactorRankingSelection = 1.5;
	public int mutationRateMax = 4;//out of len
	public double mutationProbability = 0.5;
	public double crossoverProbability = 1.0;
	
	//not on outer parameters
	public double gaussianStdDev = 50.0;
	public int replaceTournamentN = 2;
	
	public InnerParameters() {
		//TODO copy param from Parameters class xes: this.x = Parameters.x
	}
	
	public void setOutParam() {
//		Parameters.tournamentSize = this.tournamentSize;
//		Parameters.parentsN = this.parentsN;
//		Parameters.swapIndividuals = this.swapIndividuals;
//		Parameters.scalingFactorRankingSelection = this.swapIndividuals;
//		Parameters.mutationProbability = this.mutationProbability;
//		Parameters.crossoverProbability = this.crossoverProbability;
	}
}
