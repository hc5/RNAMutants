package rna;

import junit.framework.TestCase;

public class FolderTest extends TestCase {
	public void testFolding(){
		SequenceOld s = new SequenceOld("GGGCCCC");
		assertTrue(s.getMaxPairings()==3);
		assertTrue(s.getFoldedStruct().equals("(((.)))"));
		s.setBasepairs("AUGC");
		assertTrue(s.getMaxPairings()==2);
		assertTrue(s.getFoldedStruct().equals("()()"));
		s.setBasepairs("AAAAAAAA");
		assertTrue(s.getMaxPairings()==0);
		assertTrue(s.getFoldedStruct().equals("........"));
	}
}
