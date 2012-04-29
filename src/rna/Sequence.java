package rna;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static rna.Debug.print;
import static rna.Debug.log;
public class Sequence {
	protected double RT = 0.6156;
	protected String sequence;
	public int len;
	public HashMap<String, Double> loopEnergies2x2;
	public HashMap<String, Double> loopEnergies2x1;
	public HashMap<String, Double> loopEnergies1x1;

	public Sequence(String seq) {
		this.sequence = seq;
		len = sequence.length();
		this.loopEnergies2x2 = new HashMap<String, Double>();
		try {
			loadLoopEnergy2x2("INTERNAL-LOOPS2x2.dat");
			loadLoopEnergy2x1("INTERNAL-LOOPS2x1.dat");
			loadLoopEnergy1x1("INTERNAL-LOOPS1x1.dat");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/***
	 * 5' x1 --- x2 3' 
	 * 3' y1 --- y2 5'
	 * 
	 */
	protected double internalLoop(int x1, int x2, int y1, int y2) {
		double energy = 0;
		String seq = sequence.substring(x1, x2 + 1)
				+ rev(sequence.substring(y2, y1 + 1));
		int topLen = x2 - x1 - 1;
		int botLen = y1 - y2 - 1;
		if (topLen == 1 && botLen == 1) {// 1x1 loop
			return loopEnergies1x1.get(seq);
		}
		if (topLen == 2 && botLen == 2) {// 2x2 loop

			return loopEnergies2x2.get(seq);
		}
		if ((topLen == 1 && botLen == 2)) {// 1x2 loop
			return loopEnergies2x1.get(seq);
		}
		if (topLen == 2 && botLen == 1) {// 2x1 loop
			seq = rev(sequence.substring(y2, y1 + 1))
					+ (sequence.substring(x1, x2 + 1));
			try{
			return loopEnergies2x1.get(seq);
			}
			catch(Exception e){
			log(seq);
			log(topLen+","+botLen);
			log(sequence.substring(x1, x2 + 1)+"|"+rev(sequence.substring(y2, y1 + 1)));
			log(rev(sequence.substring(y2, y1 + 1))+"|"+(sequence.substring(x1, x2 + 1)));
				return loopEnergies2x1.get(seq);
			}
		}
		int n = topLen + botLen;
		int mismatch = Math.abs(topLen - botLen);
		switch (n) { // initiation
		case 4:
			energy += 1.1;
			break;
		case 5:
			energy += 2.0;
			break;
		case 6:
			energy += 2.0;
			break;
		default:
			energy += 2.0 + 1.08 * Math.log(n / 6.0);
			break;
		}

		energy += mismatch * 0.6; // asymmetry bonus
		
		
		//bad code needs refactor
		if ((topLen == 2 && botLen == 3) || (topLen == 3 && botLen == 2)) {// 2x3
																			// terminal
																			// mismatch
			if (isY(y2) && (at(y2 + 1) == 'A' && at(x2 - 1) == 'G'))
				energy -= 0.5;
			if (isY(y2) && (at(y2 + 1) == 'G' && at(x2 - 1) == 'A'))
				energy -= 1.1;
			if (!isY(y2) && (at(y2 + 1) == 'G' && at(x2 - 1) == 'A'))//not a typo
				energy -= 1.2;
			if (at(y2 + 1) == 'U' && at(x2 - 1) == 'U')
				energy -= 0.4;
			if (at(y2 + 1) == 'G' && at(x2 - 1) == 'G')
				energy -= 0.8;

			if (isY(x1) && (at(x1 + 1) == 'A' && at(y1 - 1) == 'G'))
				energy -= 0.5;
			if (isY(x1) && (at(x1 + 1) == 'G' && at(y1 - 1) == 'A'))
				energy -= 1.1;
			if (!isY(x1) && (at(x1 + 1) == 'G' && at(y1 - 1) == 'A'))//not a typo
				energy -= 1.2;
			if (at(x1 + 1) == 'U' && at(y1 - 1) == 'U')
				energy -= 0.4;
			if (at(x1 + 1) == 'G' && at(y1 - 1) == 'G')
				energy -= 0.8;
		} else if (topLen > 1 && botLen > 1) {
			if (at(x1 + 1) == 'A' && at(y1 - 1) == 'G')
				energy -= 0.8;
			if (at(x1 + 1) == 'G' && at(y1 - 1) == 'A')
				energy -= 1.0;
			if (at(x1 + 1) == 'G' && at(y1 - 1) == 'G')
				energy -= 1.2;
			if (at(x1 + 1) == 'U' && at(y1 - 1) == 'U')
				energy -= 1.2;
			
			if (at(y2 + 1) == 'A' && at(x2 - 1) == 'G')
				energy -= 0.8;
			if (at(y2 + 1) == 'G' && at(x2 - 1) == 'A')
				energy -= 1.0;
			if (at(y2 + 1) == 'G' && at(x2 - 1) == 'G')
				energy -= 1.2;
			if (at(y2 + 1) == 'U' && at(x2 - 1) == 'U')
				energy -= 1.2;
		}
		
		//AU/GU closure
		if((at(x1)=='A'&&at(y1)=='U')||(at(x1)=='U'&&at(y1)=='A')||(at(x1)=='U'&&at(y1)=='G')||(at(x1)=='G'&&at(y1)=='U'))
			energy += 0.7;
		if((at(x2)=='A'&&at(y2)=='U')||(at(x2)=='U'&&at(y2)=='A')||(at(x2)=='U'&&at(y2)=='G')||(at(x2)=='G'&&at(y2)=='U'))
			energy += 0.7;
		return energy;
	}

	protected char at(int x1) {
		return sequence.charAt(x1);
	}

	protected boolean isR(int n) {
		return isR(at(n));
	}

	protected boolean isR(char c) {
		return c == 'A' || c == 'G';
	}

	protected boolean isY(char c) {
		return !isR(c);
	}

	protected boolean isY(int n) {
		return !isR(n);
	}

	private String rev(String s) {
		return new StringBuffer(s).reverse().toString();
	}

	private void loadLoopEnergy1x1(String f) throws IOException {
		loopEnergies1x1 = new HashMap<String, Double>();
		String[] boundingPairs = new String[] { "AU", "CG", "GC", "UA", "GU",
				"UG" };
		char[] baseCoord = new char[] { 'A', 'C', 'G', 'U' };
		BufferedReader br = new BufferedReader(new FileReader(f));
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 10; j++)
				br.readLine();
			for (int k = 0; k < 4; k++) {
				String[] lineArray = (br.readLine()).trim().split(" +");
				for (int j = 0; j < lineArray.length; j++) {
					String first = boundingPairs[i];
					String second = "" + baseCoord[k] + baseCoord[j % 4];
					String third = boundingPairs[j / 4];
					String seq = "" + first.charAt(0) + second.charAt(0)
							+ third.charAt(0) + first.charAt(1)
							+ second.charAt(1) + third.charAt(1);
					loopEnergies1x1.put(seq, Double.parseDouble(lineArray[j]));
				}
			}
		}
	}

	private void loadLoopEnergy2x1(String f) throws IOException {
		loopEnergies2x1 = new HashMap<String, Double>();
		String[] boundingPairs = new String[] { "AU", "CG", "GC", "UA", "GU",
				"UG" };
		char[] baseCoord = new char[] { 'A', 'C', 'G', 'U' };
		BufferedReader br = new BufferedReader(new FileReader(f));
		for (int i = 0; i < 24; i++) {
			char base = baseCoord[i % 4];
			for (int k = 0; k < 10; k++)
				br.readLine();
			for (int k = 0; k < 4; k++) {
				String[] lineArray = (br.readLine()).trim().split(" +");
				for (int j = 0; j < lineArray.length; j++) {
					String seq = boundingPairs[i / 4].charAt(0) + "";
					seq += baseCoord[k];
					seq += boundingPairs[j / 4].charAt(0);
					seq += boundingPairs[i / 4].charAt(1) + "";
					seq += baseCoord[j % 4];
					seq += baseCoord[i % 4];
					seq += boundingPairs[j / 4].charAt(1);
					loopEnergies2x1.put(seq, Double.parseDouble(lineArray[j]));
				}
			}
		}

	}

	private void loadLoopEnergy2x2(String f) throws IOException {
		String[] coordinates = new String[] { "AA", "AC", "AG", "AU", "CA",
				"CC", "CG", "CU", "GA", "GC", "GG", "GU", "UA", "UC", "UG",
				"UU" };
		String[] boundingPairs = new String[] { "AU", "CG", "GC", "GU", "UA",
				"UG" };
		BufferedReader br = new BufferedReader(new FileReader(f));
		for (int i = 0; i < 36; i++) {
			String firstPair = boundingPairs[i / 6];
			String secondPair, thirdPair;
			String fourthPair = boundingPairs[i % 6];
			for (int j = 0; j < 10; j++)
				br.readLine();
			for (int j = 0; j < coordinates.length; j++) {
				String[] lineArray = br.readLine().trim().split(" +");
				for (int k = 0; k < lineArray.length; k++) {
					secondPair = coordinates[j];
					thirdPair = coordinates[k];
					double energy = Double.parseDouble(lineArray[k]);
					loopEnergies2x2.put(
							"" + firstPair.charAt(0) + secondPair.charAt(0)
									+ thirdPair.charAt(0)
									+ fourthPair.charAt(0)
									+ firstPair.charAt(1)
									+ secondPair.charAt(1)
									+ thirdPair.charAt(1)
									+ fourthPair.charAt(1), energy);
				}
			}
		}
	}

	protected boolean isValidPair(int x, int y) {
		char b1 = sequence.charAt(x);
		char b2 = sequence.charAt(y);
		switch (b1) {
		case 'U':
			return b2 == 'G' || b2 == 'A';
		case 'A':
			return b2 == 'U';
		case 'G':
			return b2 == 'U' || b2 == 'C';
		case 'C':
			return b2 == 'G';
		}
		return false;
	}

	protected double hairpin(int x, int y) {
		int length = y - x - 1;
		if(length<3)
			return Double.MAX_VALUE;
		double energy;
		
		double[] energies = new double[] { Double.MAX_VALUE,
				Double.MAX_VALUE,Double.MAX_VALUE, 5.4, 5.6, 5.7, 5.4, 6.0, 5.5, 6.4};
		if(length<=9)
			energy = energies[length];
		else
			energy = energies[9] + 1.75 * RT * Math.log(length/9.0);
		HashMap<String, Double> specialCases = new HashMap<String, Double>() {
			{
				put("CAACG", 6.8);
				put("GUUAC", 6.9);
				put("CUACGG", 2.8);
				put("CUCCGG", 2.7);
				put("CUUCGG", 3.7);
				put("CUUUGG", 3.7);
				put("CCAAGG", 3.3);
				put("CCCAGG", 3.4);
				put("CCGAGG", 3.5);
				put("CCUAGG", 3.7);
				put("CCACGG", 3.7);
				put("CCGCGG", 3.6);
				put("CCUCGG", 2.5);
				put("CUAAGG", 3.6);
				put("CUCAGG", 3.7);
				put("CUUAGG", 3.5);
				put("CUGCGG", 2.8);
				put("CAACGG", 5.5);
				put("ACAGUGC", 2.9);
				put("ACAGUGAU", 3.6);
				put("ACAGUGUU", 1.8);
				put("ACAGUACU", 2.8);
			}
		};
		String seq = sequence.substring(x, y + 1);
		String loopSeq = sequence.substring(x + 1, y);
		if (specialCases.containsKey(seq))
			energy = specialCases.get(seq);
		if ((sequence.charAt(x) == 'G' && sequence.charAt(y) == 'U')
				|| (sequence.charAt(x) == 'U' && sequence.charAt(y) == 'G'))
			energy -= 2.2;
		if (length > 3) {
			if (at(x) == 'G' && sequence.charAt(y - 1) == 'G')
				energy -= 0.8;
			if ((at(x) == 'U' && sequence.charAt(y - 1) == 'U')
					|| (at(x) == 'A' && sequence.charAt(y - 1) == 'G')
					|| (at(x) == 'G' && sequence.charAt(y - 1) == 'A'))
				energy -= 0.9;
			if (allC(loopSeq)) {
				energy += 0.3 * length + 1.6;
			}
		}
		if (loopSeq.equals("CCC"))
			energy += 1.5;

		return energy;
	}

	private boolean allC(String s) {
		boolean allC = true;
		for (char c : s.toCharArray())
			allC &= c == 'C';
		return allC;
	}

	/***
	 * 5' x1 x2 3'
	 * 3' y1 y2 5'
	 * 
	 */

	protected double stack(int x1, int x2, int y1, int y2) {
		String bases = "" + at(x1) + at(x2) + at(y1) + at(y2);
		HashMap<String, Double> energies = new HashMap<String, Double>() {
			{
				put("AAUU", -0.93);
				put("AUUA", -1.10);
				put("AUUG", -1.4);
				put("UUAA", -0.93);
				put("UUAG", -1.3);
				put("UGAU", -1.0);
				put("UAGU", -1.0);
				put("UCGG", -1.5);
				put("UGGC", -1.4);
				put("UGGU", 0.3);
				put("UUGA", -0.6);
				put("UUGG", -0.5);
				put("UAAU", -1.33);
				put("CUGA", -2.08);
				put("CUGG", -2.1);
				put("AGUC", -2.08);
				put("AGUU", -0.6);
				put("CAGU", -2.11);
				put("UGAC", -2.11);
				put("GUCA", -2.24);
				put("GUCG", -2.5);
				put("ACUG", -2.24);
				put("GACU", -2.35);
				put("UCAG", -2.35);
				put("CGGC", -2.36);
				put("CGGU", -1.4);
				put("GGCC", -3.26);
				put("GGUU", -0.5);
				put("GGUC", -2.1);
				put("GUUG", 1.3);
				put("GAUU", -1.3);
				put("GUUA", -1.4);
				put("GCUG", -2.5);
				put("GGCU", -3.3);
				put("CCGG", -3.26);
				put("GCCG", -3.42);
			}
		};
		return energies.get(bases);
	}
	
	
	/***
	 * 5' x1 --- x2 3' 
	 * 3' y1 --- y2 5'
	 * 
	 */
	protected double bulge(int x1, int x2, int y1, int y2) {
		double[] initiation = new double[] { 3.81, 2.8, 3.2,3.6,4.0,4.4 };
		double energy = 0.0;
		int bulgeLength = Math.max(x2-x1-1, y1 - y2-1);
		if(bulgeLength<=6)
			energy += initiation[bulgeLength-1];
		else
			energy += 4.4 + 1.75*RT * Math.log(bulgeLength/6.0);
		if(bulgeLength == 1){
			char e1,e2,middle;
			if(x2-x1 == 2){
				e1 = at(x1);
				e2 = at(x2);
				middle = at(x1 +1);
			}
			else{
				e1 = at(y2);
				e2 = at(y1);
				middle = at(y2 +1);
			}
			if(middle == 'C' && (e1 == 'C' || e2 == 'C'))
				energy -= 0.9;
		}
		else{
			//AU/GU closure
			if((at(x1)=='A'&&at(y1)=='U')||(at(x1)=='U'&&at(y1)=='A')||(at(x1)=='U'&&at(y1)=='G')||(at(x1)=='G'&&at(y1)=='U'))
				energy += 0.45;
			if((at(x2)=='A'&&at(y2)=='U')||(at(x2)=='U'&&at(y2)=='A')||(at(x2)=='U'&&at(y2)=='G')||(at(x2)=='G'&&at(y2)=='U'))
				energy += 0.45;
		}
		return energy;
	}
}
