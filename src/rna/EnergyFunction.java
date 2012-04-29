package rna;

public abstract class EnergyFunction extends Sequence{
	protected double F[][];
	protected double C[][];
	protected double M[][];
	protected double M1[][];
	protected boolean calculatedF[][];
	protected boolean calculatedC[][];
	protected boolean calculatedM[][];
	protected boolean calculatedM1[][];
	public EnergyFunction(String seq) {
		super(seq);
		int l = seq.length();
		F = new double[l][l];
		C = new double[l][l];
		M = new double[l][l];
		M1 = new double[l][l];
		calculatedF = new boolean[l][l];
		calculatedC = new boolean[l][l];
		calculatedM = new boolean[l][l];
		calculatedM1 = new boolean[l][l];
	}
	protected String repeat(String s,int n){
		return new String(new char[n]).replace("\0", s);
	}
	protected abstract String solve();
	protected abstract String backtrackF(int i, int j);
	protected abstract String backtrackC(int i, int j) ;
	protected abstract String backtrackM(int i, int j) ;
	protected abstract String backtrackM1(int i, int j);
	protected abstract double getC(int i, int j);
	protected abstract double getF(int i, int j);
	protected abstract double getM(int i, int j) ;
	protected abstract double getM1(int i, int j);
	protected void set(int i,int j, double m,double [][] d, boolean [][] mset){
		d[i][j]=m;
		mset[i][j]= true;
	}
	protected void setF(int i, int j, double m) {
		F[i][j] = m;
		calculatedF[i][j] = true;
	}

	protected void setC(int i, int j, double m) {
		C[i][j] = m;
		calculatedC[i][j] = true;
	}

	protected void setM(int i, int j, double m) {
		M[i][j] = m;
		calculatedM[i][j] = true;
	}

	protected void setM1(int i, int j, double m) {
		M1[i][j] = m;
		calculatedM1[i][j] = true;
	}
	protected String folded;
}
