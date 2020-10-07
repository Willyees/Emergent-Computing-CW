package ea;
import java.util.ArrayList;

public class Island {
	InnerParameters parameters = new InnerParameters();
	ArrayList<Individual> population = new ArrayList<Individual>();
	String mutation;
	String replacement;
	String selection;
	String initialization;
	String crossover;
	//int type;
	
	public Island() {
		//to be removed, just quick debugging
	}
	
	public Island(String mutation, String replacement, String selection, String initialization, String crossover) {
		this.mutation = mutation;
		this.replacement = replacement;
		this.selection = selection;
		this.initialization = initialization;
		//this.type = type;
	}
	
}
