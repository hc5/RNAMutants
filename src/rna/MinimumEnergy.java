package rna;

import static rna.Debug.print;
import static rna.Debug.log;

public class MinimumEnergy extends EnergyFunction{
	

	public MinimumEnergy(String seq) {
		super(seq);
		double defaultVal = 0;
		try{
		for (int i = 0; i < len; i++) {
			setF(i, i, defaultVal);
			setF(i, i + 1, defaultVal);
			setC(i, i, Double.MAX_VALUE);
			setC(i, i + 1, Double.MAX_VALUE);
			setM(i, i, Double.MAX_VALUE);
			setM(i, i + 1, Double.MAX_VALUE);
			setM1(i, i, Double.MAX_VALUE);
			setM1(i, i + 1, Double.MAX_VALUE);
		}
		}
		catch(ArrayIndexOutOfBoundsException e){
		}
		
	}

	

	protected String solve() {
		getF(0, sequence.length() - 1);
		return backtrackF(0,len-1);
	}

	protected String backtrackF(int i, int j) {
		if(j<i)
			return "";
		String seq = ".";
		if(j == i)
			return seq;
		double min = F[i][j];
		if (min == F[i][j-1]){
			return  backtrackF(i , j-1)+".";
		}
		for (int k = i; k < j; k++) {
			if (!isValidPair(k, j))
				continue;
			double term = getC(k, j) + getF(i, k-1);
			if (min == term) {
				return  backtrackF(i, k-1)+backtrackC(k, j);
			}
		}
		return seq;
	}

	protected String backtrackC(int i, int j) {
		
		String seq = "(" + repeat(".", j - i - 1) + ")";
		double min = C[i][j];
		if (min == hairpin(i, j)){
			//log("hairpin");
			return seq;
		}
		for (int k = i + 1; k < j - 1; k++) {
			for (int l = k + 1; l < j; l++) {
				if (isValidPair(k, l)) {
					double term = C[k][l];
					if (k == i + 1 && l == j - 1) {
						//log("stack + %f",stack(i, k, j, l));
						term += stack(i, k, j, l);
					}
					else if (k - i >= 2 && j-l >= 2){
						//log("internalLoop + %f",internalLoop(i, k, j,l));
						term += internalLoop(i, k, j,l);
					}
					else if (k - i >= 2 || j-l >= 2){
						term += bulge(i, k, j,l);
					}
					if (min == term) {
						return "(" + repeat(".", k - i - 1)
								+ backtrackC(k, l)
								+ repeat(".", j - l - 1) + ")";
					}
				}
			}
		}
		for (int k = i + 2; k < j - 1; k++) {
			double term = getM(i + 1, k) + getM1(k + 1, j - 1);
			
			if (min == term){
				return "("+backtrackM(i + 1, k) + backtrackM1(k + 1, j - 1)+")";
			}
		}
		return seq;
	}

	protected String backtrackM(int i, int j) {
		String seq = "";
		double min = M[i][j];
		
		for (int k = i + 1; k < j; k++) {
			double term = getM1(k, j);
			if (term == min){
				
				return repeat(".", k - i - 1) + backtrackM1(k, j);
			}
		}
		for (int k = i + 2; k < j; k++) {
			double term = getM(i, k - 1) + getM1(k, j);
			if (term == min){
				return backtrackM(i, k - 1) + backtrackM1(k, j);
			}
		}
		return seq;
	}

	protected String backtrackM1(int i, int j) {
		String seq = "";
		double min = M1[i][j];
		for (int k = i + 1; k < j; k++) {
			if (isValidPair(i, k)) {
				double term = getC(i, k);
				if (term == min){
					return backtrackC(i, k) + repeat(".", j - k - 1);
				}
			}
		}
		return seq;
	}

	protected double getC(int i, int j) {
		if (calculatedC[i][j])
			return C[i][j];
		calculatedC[i][j] = true;
		double min = Double.MAX_VALUE;
		if (hairpin(i, j) < min)
			min = hairpin(i, j);
		for (int k = i + 1; k < j - 1; k++) {
			for (int l = k + 1; l < j; l++) {
				if (isValidPair(k, l)) {
					double term = getC(k, l);
					if (k == i + 1 && l == j - 1) {
						term += stack(i, k, j, l);
					}
					else if (k - i >= 2 && j-l >= 2){
						term += internalLoop(i, k, j,l);
					}
					else if (k - i >= 2 || j-l >= 2){
						term += bulge(i, k, j,l);
					}
					if (term <= min)
						min = term;
				}
			}
		}
		for (int k = i + 2; k < j - 1; k++) {
			double term = getM(i + 1, k) + getM1(k + 1, j - 1);
			if (term < min)
				min = term;
		}
		C[i][j] = min;
		return min;
	}

	protected double getF(int i, int j) {
		if(j<i)
			return 0;
		if (calculatedF[i][j])
			return F[i][j];
		calculatedF[i][j] = true;
		double min = getF(i , j-1);
		for (int k = i; k < j; k++) {
			if (!isValidPair(k, j))
				continue;
			double term = getC(k, j) + getF(i, k-1);
			if (term < min)
				min = term;
		}
		F[i][j] = min;
		return min;
	}

	protected double getM(int i, int j) {
		if (calculatedM[i][j])
			return M[i][j];
		calculatedM[i][j] = true;
		double min = Double.MAX_VALUE;
		for (int k = i + 1; k < j; k++) {
			double term = getM1(k, j);
			if (term < min)
				min = term;
		}
		for (int k = i + 2; k < j; k++) {
			double term = getM(i, k - 1) + getM1(k, j);
			if (term < min)
				min = term;
		}
		M[i][j] = min;
		return min;
	}

	protected double getM1(int i, int j) {
		if (calculatedM1[i][j]){
			return M1[i][j];
		}
		calculatedM1[i][j] = true;
		double min = Double.MAX_VALUE;
		for (int k = i + 1; k < j; k++) {
			if (isValidPair(i, k)) {
				double term = getC(i, k);
				if (term < min)
					min = term;
			}
		}
		M1[i][j] = min;
		return min;
	}

	
}
