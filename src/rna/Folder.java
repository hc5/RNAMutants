package rna;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import static rna.Debug.log;
import static rna.Debug.print;
public class Folder {

	public static void main(String[] args) throws IOException {
		try{
		boolean determ = false;
		int numTries=1000, listLimit=10;
		if(args[0].equals("-d"))
			determ = true;
		else if(args[0].equals("-s")){
			determ = false;
			numTries = Integer.parseInt(args[1]);
			listLimit = Integer.parseInt(args[2]);
		}
		System.out.println("Please enter the sequence");
		BufferedReader br = new BufferedReader(new InputStreamReader((System.in)));
		String sequence = br.readLine();
		if(determ){
			MinimumEnergy me = new MinimumEnergy(sequence);
			System.out.println(me.solve());
		}
		else{
			PartitionFunction pf = new PartitionFunction(sequence);
			pf.solve();
			pf.list(numTries,listLimit);
		}
		}
		catch(Exception e){
			System.out.println("USAGE: [-d|-s] [n] [l]\n-d: deterministic backtracking\n-s: stochastic backtracking with n tries, and list the top l number of configurations");
		}
		//System.out.println(me.runRandom(1));
	}
}
