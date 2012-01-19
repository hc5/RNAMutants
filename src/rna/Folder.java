package rna;

public class Folder {

	public static void main(String [] args){
		Sequence s = new Sequence("AUGCAUUUGUUUGACC");
		s.setBeta(0.5);
		System.out.println(s.getMaxPairings());
		System.out.println(s.getFoldedStruct());
		
	}
}
