
public class Folder {

	public static void main(String [] args){
		Sequence s = new Sequence("GGGAAAUCC");
		s.fold();
		System.out.println(s.getMaxPairings());
		
	}
}
