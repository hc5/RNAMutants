package rna;

public class Folder {

	public static void main(String [] args){
		Sequence s = new Sequence("AUGC");
		System.out.println(s.getMaxPairings());
		System.out.println(s.getFoldedStruct());
		
	}
}
