package rna;

import java.util.HashMap;

public class SequenceOld {
	private final boolean DEBUG = true;
	private String basepairs;
	private String foldedStruct;
	private int maxPairings;
	private int [][] maxBasepairs = null;
	private double beta = 1/(8.3144621*310.2);// 1/RT, T = 37 C
	private double Z [][] = null;
	private double Zb [][] = null;
	private double ZM1 [][] = null;
	private double ZM [][]= null;
	public double getBeta() {
		return beta;
	}

	public void setBeta(double beta) {
		this.beta = beta;
	}
	public String getBasepairs() {
		return basepairs;
	}
	public void setBasepairs(String basepairs) {
		this.basepairs = basepairs.toLowerCase();
		this.foldedStruct = "";
		this.maxPairings = -1;
		this.maxBasepairs = null;
		fold();
	}

	public SequenceOld(String basepairs) {
		this.basepairs = basepairs.toLowerCase();
	}
	public String getFoldedStruct(){
		if(this.maxBasepairs==null)
			fold();
		return foldedStruct;
	}
	public int getMaxPairings(){
		if(this.maxBasepairs==null)
			fold();
		return maxPairings;
	}
	public void fold(){
		getPairingMatrix();
		maxPairings = maxBasepairs[0][maxBasepairs.length-1];
		//print(partitionMatrix);
	}

	//currently it's just the max number of bp between i and j
	public double getFreeEnergy(int i, int j){
		if(maxBasepairs[i][j]==0)
			return Double.MAX_VALUE;
		return 1.0 * maxBasepairs[i][j];
	}

	public HashMap<String,Double> sample(int n){
		calculatePartitionMatrix();
		HashMap<String,Double> freq = new HashMap<String,Double>();
		for(int i=0;i<n;i++){
			String seq = stochBacktrack(0,length()-1);
			if(freq.containsKey(seq)){
				freq.put(seq, freq.get(seq).doubleValue()+1.0/n);
			}
			else
				freq.put(seq, 1.0/n);
		}
		return freq;
	}

	public void calculatePartitionMatrix(){
		//initialize matricies
		int l = length();
		Z = new double[l][l];
		Zb = new double[l][l];
		ZM1 = new double[l][l];
		ZM = new double[l][l];
		for(int i =0;i<l;i++){
			Z[i][i] = 1;			
		}
		for(int end=1;end<l;end++){
			for(int start =end-1;start>-1;start--){
				double total = 0;
				total +=Z[start+1][end];
				for(int k = start;k<=end;k++){
					total += Zb[start][k]*Z[start+1][end];					
				}
			}
		}
	}



	public double energyFunction(double E){
		return Math.exp(this.beta * -1*E * -1);
	}


	public int length(){
		return this.basepairs.length();		
	}
	public String determBacktrack(int start,int end){
		String seq="";
		int max = this.maxBasepairs[start][end];
		if(start == end)
			return ".";
		if(start > end)
			return "";
		if(isValidPair(start,end) && this.maxBasepairs[start+1][end-1]+1 == max)
			return "("+determBacktrack(start+1,end-1)+")";
		if(this.maxBasepairs[start][end-1]==max)
			return determBacktrack(start,end-1)+".";
		if(this.maxBasepairs[start+1][end]==max)
			return "."+determBacktrack(start+1,end);
		for(int k = start +1;k<end;k++){
			if(this.maxBasepairs[start][k]+this.maxBasepairs[k+1][end]==max)
				return determBacktrack(start,k)+determBacktrack(k+1,end);
		}
		return seq;
	}
	private double E(String s){
		double e = 0;
		return e;
	}
	public String stochBacktrack(int start, int end){
		String seq="";
		double randCutoff = randChoose(Z[start][end]);
		if(start == end)
			return ".";
		if(start > end)
			return "";
		double base = 0;
		if(isValidPair(start,end)){
			base += Z[start+1][end-1] * energyFunction(1) ;
			if(base >= randCutoff){
				return "("+stochBacktrack(start+1,end-1)+")";
			}
		}
		base += Z[start][end-1];
		if(base >= randCutoff){
			return stochBacktrack(start,end-1)+".";
		}
		for(int k =start+1;k<end;k++){
			if(isValidPair(k,end)){
				double partial =Z[start][k-1] * Z[k+1][end] * energyFunction(1);
				base += partial;
				if(base  >= randCutoff){
					return stochBacktrack(start,k-1)+"("+stochBacktrack(k+1,end)+")";
				}
			}
		}
		return seq;
	}

	public double randChoose(double total){
		return Math.random() * total;		
	}

	private void getPairingMatrix(){
		int [][] m = new int [basepairs.length()][basepairs.length()];
		//initialization
		for(int i = 0;i<m.length;i++){
			m[i][i]=0;
			if(i>0){
				m[i][i-1]=0;
			}
			if(i<m.length-1){
				if(isValidPair(i,i+1)){
					m[i][i+1]=1;
				}
				else{
					m[i][i+1] = 0;
				}
			}

		}
		for(int end=1;end<m.length;end++){
			for(int start =end-1;start>-1;start--){
				int max = 0;

				if(isValidPair(start,end)){
					if(m[start+1][end-1]+1>max)
						max = m[start+1][end-1]+1;
				}

				if(m[start+1][end]>max){
					max = m[start+1][end];
				}
				if(m[start][end-1]>max){
					max = m[start][end-1];
				}
				for(int i =start+1;i<end;i++){
					int curBpCount = m[start][i]+m[i+1][end];
					if(curBpCount>max){
						max = curBpCount;
					}
				}

				m[start][end]=max;
			}
		}
		this.maxBasepairs=m;
		this.foldedStruct=determBacktrack(0,length()-1);
	}
	private String repeat(String c,int n){
		return new String(new char[n]).replace("\0", c);

	}
	public double [][] getPartitionMatrix(){
		if(Z == null){
			this.fold();
		}
		return this.Z;
	}

	private boolean isValidPair(int x,int y){
		char b1 = basepairs.charAt(x);
		char b2 = basepairs.charAt(y);
		switch(b1){
		case 'u':
			return b2=='g'||b2=='a';
		case 'a':
			return b2=='u';
		case 'g':
			return b2=='u'||b2=='c';
		case 'c':
			return b2=='g';
		}
		return false;
	}
	@Override
	public String toString() {
		return "Sequence [basepairs=" + basepairs + "]";
	}



	/**** 
	 * 
	 * 
	 * all print array functions 
	 *  
	 * 
	 * ***/



}
