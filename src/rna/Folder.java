package rna;

import java.util.HashMap;

public class Folder {

	public static void main(String [] args){
		Sequence s = new Sequence("A.....U");
		s.fold();
		HashMap<String,Double> f = (s.sample(10000));
		for(String k : f.keySet()){
			System.out.printf("%s : %2f\n", k, f.get(k));
		}
		System.out.println(s.getFoldedStruct());
		
	}
}
