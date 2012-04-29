package rna;

import static rna.Debug.log;
import static rna.Debug.print;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;


public class PartitionFunction extends EnergyFunction {
	// the fields/methods starting with 'd' are deterministic algorithms - should be equivalent to MinimumEnergy
	// the ones not starting with 'd' are stochastic ones
	private double[][] dF;
	private double[][] dC;
	private double[][] dM;
	private double[][] dM1;
	private boolean[][] dcalculatedF;
	private boolean[][] dcalculatedC;
	private boolean[][] dcalculatedM;
	private boolean[][] dcalculatedM1;
	public PartitionFunction(String seq) {
		super(seq);
		int l = len;
		dF = new double[l][l];
		dC = new double[l][l];
		dM = new double[l][l];
		dM1 = new double[l][l];
		dcalculatedF = new boolean[l][l];
		dcalculatedC = new boolean[l][l];
		dcalculatedM = new boolean[l][l];
		dcalculatedM1 = new boolean[l][l];
		double defaultVal = 0;
		try{
			for (int i = 0; i < len; i++) {
				set(i,i,1,dF,dcalculatedF);
				set(i,i+1,1,dF,dcalculatedF);
				set(i,i,0,dC,dcalculatedC);
				set(i,i+1,0,dC,dcalculatedC);
				set(i,i,0,dM,dcalculatedM);
				set(i,i+1,0,dM,dcalculatedM);
				set(i,i,0,dM1,dcalculatedM1);
				set(i,i+1,0,dM1,dcalculatedM1);
				setF(i, i, 1);
				setF(i, i + 1, 1);
				setC(i, i, 0);
				setC(i, i + 1, 0);
				setM(i, i, defaultVal);
				setM(i, i + 1, defaultVal);
				setM1(i, i, defaultVal);
				setM1(i, i + 1, defaultVal);
			}
		}
		catch(ArrayIndexOutOfBoundsException e){
		}
	}
	private double e(double d){
		return Math.exp(-d*1000/(1.9858775*310.2));
	}
	private double rand(double d){
		return Math.random() * d;
	}

	public double[][] determSolve(){
		dgetF(0,len-1);
		return dF;
	}
	private void testRand(){//test the correctness of rand
		int dist [] = new int[]{0,0,0,0,0,0,0,0,0,0};
		for(int i = 0;i<100000000;i++){
			double d = rand(1000);
			int di = (int)d;
			dist[di/100]++;
		}
		for(int i = 0;i<10;i++){
			System.out.printf("%d - %d: %d\n", i,(i+1)*100,dist[i]);
		}
	}
	@Override
	public String solve() {
		getF(0, sequence.length() - 1);
		return backtrackF(0,len-1);
	}

	@Override
	protected String backtrackF(int i, int j) {
		if(j<i)
			return "";
		String seq = ".";
		if(j == i)
			return seq;
		double val = F[i][j];
		double r = rand(val);
		double base = F[i][j-1];
		log("F-1:"+base+"\t"+base+"\t"+r+"\t"+val);
		if (r<=base){
			return  backtrackF(i , j-1)+".";
		}
		for (int k = i; k < j; k++) {
			if (!isValidPair(k, j))
				continue;
			double term =getC(k, j) * getF(i, k-1); 
			log("is valid pair:"+term+"\t"+base+"\t"+r+"\t"+val);
			base += term;
			if (r<=base) {
				return  backtrackF(i, k-1)+backtrackC(k, j);
			}
		}

		return seq;
	}

	@Override
	protected String backtrackC(int i, int j) {
		String seq = "(" + repeat(".", j - i - 1) + ")";
		double val = C[i][j];
		double r = rand(val);
		double base = e(hairpin(i, j));
		if(r<=base){
			log("hairpin = "+base+"/"+r+"/"+val);
			return seq;
		}
		for (int k = i + 1; k < j - 1; k++) {
			for (int l = k + 1; l < j; l++) {
				if (isValidPair(k, l) && l-k-1>=3) {
					double term = getC(k, l);
					if (k == i + 1 && l == j - 1) {
						term *= e(stack(i, k, j, l));
						log("stack");
					}
					else if (k - i >= 2 && j-l >= 2){
						log("internalLoop");
						term *= e(internalLoop(i, k, j,l));
					}
					else if (k - i >= 2 || j-l >= 2){
						log("bulge");
						term *= e(bulge(i, k, j,l));
					}

					base += (term);

					log(term+"\t"+(term)+"\t"+r+"\t"+val);
					if(r<=base)
						return "(" + repeat(".", k - i - 1)
								+ backtrackC(k, l)
								+ repeat(".", j - l - 1) + ")";
				}
			}
		}
		for (int k = i + 2; k < j - 1; k++) {
			base += getM(i + 1, k) * getM1(k + 1, j - 1);
			if(r<=base)
				return "("+backtrackM(i + 1, k) + backtrackM1(k + 1, j - 1)+")";
		}
		return seq;
	}

	@Override
	protected String backtrackM(int i, int j) {
		String seq = "";
		double val = M[i][j];
		double r = rand(val);
		double base = 0;
		for (int k = i + 1; k < j; k++) {
			base += getM1(k, j);
			if (r<=base){

				return repeat(".", k - i - 1) + backtrackM1(k, j);
			}
		}
		for (int k = i + 2; k < j; k++) {
			base += getM(i, k - 1) * getM1(k, j);
			if (r<=base){
				return backtrackM(i, k - 1) + backtrackM1(k, j);
			}
		}
		log("M "+base+"/"+r+"/"+val);
		return seq;
	}

	@Override
	protected String backtrackM1(int i, int j) {
		String seq = "";
		double val = M1[i][j];
		double r = rand(val);
		double base = 0;
		for (int k = i + 1; k < j; k++) {
			if (isValidPair(i, k)) {
				base += getC(i, k);
				if (r<base){
					return backtrackC(i, k) + repeat(".", j - k - 1);
				}
			}
		}
		log("M1 "+base+"/"+r+"/"+val);
		return seq;
	}

	protected String dbacktrackF(int i, int j) {
		if(j<i)
			return "";
		String seq = ".";
		if(j == i)
			return seq;
		double max = dF[i][j];
		double v = dF[i][j-1];
		if(max==v){
			return  dbacktrackF(i , j-1)+".";
		}
		for (int k = i; k < j; k++) {
			if (!isValidPair(k, j))
				continue;

			double term =dgetC(k, j) * dgetF(i, k-1); 
			if(max == term){
				return  dbacktrackF(i, k-1)+dbacktrackC(k, j);
			}
		}
		return seq;
	}

	protected String dbacktrackC(int i, int j) {
		String seq = "(" + repeat(".", j - i - 1) + ")";
		double max = dC[i][j];
		if(max == e(hairpin(i,j)))
			return seq;

		for (int k = i + 1; k < j - 1; k++) {
			for (int l = k + 1; l < j; l++) {
				if (isValidPair(k, l) && l-k-1>=3) {
					double term = dgetC(k, l);

					if (k == i + 1 && l == j - 1) {
						term += stack(i, k, j, l);
					}
					else if (k - i >= 2 && j-l >= 2){
						term += internalLoop(i, k, j,l);
					}
					else if (k - i >= 2 || j-l >= 2){
						term += bulge(i, k, j,l);
					}
					if(term == max)
						return "(" + repeat(".", k - i - 1)
								+ dbacktrackC(k, l)
								+ repeat(".", j - l - 1) + ")";
				}
			}
		}
		for (int k = i + 2; k < j - 1; k++) {
			double term = dgetM(i + 1, k) * dgetM1(k + 1, j - 1);
			if(term==max)
				return "("+dbacktrackM(i + 1, k) + dbacktrackM1(k + 1, j - 1)+")";
		}
		return seq;
	}

	protected String dbacktrackM(int i, int j) {
		String seq = "";
		double max = dgetM(i,j);

		for (int k = i + 1; k < j; k++) {
			double term = dgetM1(k, j);
			if(term==max)
				return repeat(".", k - i - 1) + dbacktrackM1(k, j);
		}
		for (int k = i + 2; k < j; k++) {
			double term =  dgetM(i, k - 1) * dgetM1(k, j);
			if(term==max)
				return dbacktrackM(i, k - 1) + dbacktrackM1(k, j);
		}
		return seq;
	}

	protected String dbacktrackM1(int i, int j) {
		String seq = "";
		double max = dgetM1(i,j);
		for (int k = i + 1; k < j; k++) {
			if (isValidPair(i, k)) {
				double term = dgetC(i, k);
				if (term==max){
					return dbacktrackC(i, k) + repeat(".", j - k - 1);
				}
			}
		}
		return seq;
	}

	protected double dgetC(int i, int j) {
		if (dcalculatedC[i][j])
			return dC[i][j];
		dcalculatedC[i][j] = true;
		double max = e(hairpin(i, j));


		for (int k = i + 1; k < j - 1; k++) {
			for (int l = k + 1; l < j; l++) {
				if (isValidPair(k, l) && l-k-1>=3) {
					double term = dgetC(k, l);
					if (k == i + 1 && l == j - 1) {
						term += stack(i, k, j, l);
					}
					else if (k - i >= 2 && j-l >= 2){
						term += internalLoop(i, k, j,l);
					}
					else if (k - i >= 2 || j-l >= 2){
						term += bulge(i, k, j,l);
					}

					if(term>max)
						max = term;

				}
			}
		}
		for (int k = i + 2; k < j - 1; k++) {
			double term = dgetM(i + 1, k) * dgetM1(k + 1, j - 1);
			if(term>max)
				max = term;
		}
		dC[i][j] = max;
		return max;
	}

	protected double dgetF(int i, int j) {
		if(j<i)
			return 1;
		if (dcalculatedF[i][j])
			return dF[i][j];
		dcalculatedF[i][j] = true;
		double max = 0;
		double v = dgetF(i,j-1);
		if (v>max){
			max = v;
		}
		for (int k = i; k < j; k++) {
			if (!isValidPair(k, j))
				continue;
			double term =dgetC(k, j) * dgetF(i, k-1); 

			if(term>max)
				max = term;
		}
		dF[i][j] = max;
		return max;
	}

	protected double dgetM(int i, int j) {
		if (dcalculatedM[i][j])
			return dM[i][j];
		dcalculatedM[i][j] = true;
		double max = 0;
		for (int k = i + 1; k < j; k++) {
			double term = dgetM1(k, j);
			if (term>max){
				max = term;

			}
		}
		for (int k = i + 2; k < j; k++) {
			double term =  dgetM(i, k - 1) * dgetM1(k, j);
			if (term>max){
				max = term;
			}
		}
		dM[i][j] = max;
		return max;
	}

	protected double dgetM1(int i, int j) {

		if (dcalculatedM1[i][j]){
			return dM1[i][j];
		}
		dcalculatedM1[i][j] = true;
		double max = 0;
		for (int k = i + 1; k < j; k++) {
			if (isValidPair(i, k)) {
				double term = dgetC(i, k);
				if (term>max){
					max = term;
				}
			}
		}
		dM1[i][j] = max;
		return max;
	}


	@Override
	protected double getC(int i, int j) {
		if (calculatedC[i][j])
			return C[i][j];
		calculatedC[i][j] = true;
		double val = e(hairpin(i, j));
		
		for (int k = i + 1; k < j - 1; k++) {
			for (int l = k + 1; l < j; l++) {
				if (isValidPair(k, l)&& l-k-1>=3) {
					double term = getC(k, l);
					if (k == i + 1 && l == j - 1) {
						term *= e(stack(i, k, j, l));
						//	log("stack=%f",stack(i,k,j,l));
					}
					else if (k - i >= 2 && j-l >= 2){
						term *= e(internalLoop(i, k, j,l));
						//	log("internalLoop=%f",internalLoop(i,k,j,l));
					}
					else if (k - i >= 2 || j-l >= 2){
						term *= e(bulge(i, k, j,l));
						//	log("bulge=%f",bulge(i,k,j,l));
					}
					//log("term = "+e(term));
					val += (term);
				}
			}
		}
		for (int k = i + 2; k < j - 1; k++) {
			double term = getM(i + 1, k) * getM1(k + 1, j - 1);
			val += term;
		}
		C[i][j] = val;
		return val;
	}

	@Override
	protected double getF(int i, int j) {
		if(j<i)
			return 1;
		if (calculatedF[i][j])
			return F[i][j];
		calculatedF[i][j] = true;
		double val = getF(i , j-1);
		//log("%d - %d",i,j);
		for (int k = i; k < j; k++) {
			if (!isValidPair(k, j))
				continue;
			val += getC(k, j) * getF(i, k-1);

			//	log("%f,%f,%f",getC(k,j),getF(i,k-1),getC(k, j) * getF(i, k-1));
		}
		F[i][j] = val;
		return val;
	}

	@Override
	protected double getM(int i, int j) {
		if (calculatedM[i][j])
			return M[i][j];
		calculatedM[i][j] = true;
		double val = 0;
		for (int k = i + 1; k < j; k++) {
			val += getM1(k, j);
		}
		for (int k = i + 2; k < j; k++) {
			val += getM(i, k - 1) * getM1(k, j);
		}
		M[i][j] = val;
		return val;
	}

	@Override
	protected double getM1(int i, int j) {
		if (calculatedM1[i][j]){
			return M1[i][j];
		}
		calculatedM1[i][j] = true;
		double val = 0;
		for (int k = i + 1; k < j; k++) {
			if (isValidPair(i, k)) {
				val += getC(i, k);
			}
		}
		M1[i][j] = val;
		return val;
	}
	private <K,V extends Comparable<? super V>> SortedSet<Map.Entry<K,V>> sortByFreq(Map<K,V> map) {
		SortedSet<Map.Entry<K,V>> sortedEntries = new TreeSet<Map.Entry<K,V>>(
				new Comparator<Map.Entry<K,V>>() {
					@Override 
					public int compare(Map.Entry<K,V> e1, Map.Entry<K,V> e2) {
						int res = e2.getValue().compareTo(e1.getValue());
						return res != 0 ? res : 1;
					}
				}
				);
		sortedEntries.addAll(map.entrySet());
		return sortedEntries;
	}
	public SortedSet<Map.Entry<String,Integer>> runRandom(final int n){
		final Map<String,Integer> results = new HashMap<String,Integer>();
		final int threadCount = 10;
		final int threadLoad = n/threadCount;
		if(n<threadCount){
			for(int i = 0;i<n;i++){
				String seq = backtrackF(0, len - 1);
				Integer freq = results.get(seq);
				results.put(seq, freq == null?(1):freq+(1));
			}
			return sortByFreq(results);
			
		}
		class workerThread extends Thread{
			
			public void run(){
				
					for(int i=0;i<threadLoad;i++){
						String seq = backtrackF(0, len - 1);
						Integer freq = r.get(seq);
						r.put(seq, freq == null?(1):freq+(1));
					}
				
			}
			public HashMap<String,Integer> r = new HashMap<String,Integer>();
		
	}
		workerThread ts[] = new workerThread[threadCount];
		
		for(int i = 0 ;i <ts.length;i++){
			ts[i] = new workerThread();
			ts[i].start();
		}
		
		for(int i =0;i<threadCount;i++){
			try {
				ts[i].join();
				for(Map.Entry<String,Integer> entry:ts[i].r.entrySet()){
					Integer count = entry.getValue();
					Integer existing = results.get(entry.getKey());
					results.put(entry.getKey(), existing==null?count:existing+count);
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		return sortByFreq(results);
	}
	public void list(int n){
		for(Map.Entry entry : runRandom(n)){
			System.out.println(entry.getKey()+":"+entry.getValue());
		}
	}
	public void list(int n, int limit){
		int i = 0;
		for(Map.Entry entry : runRandom(n)){
			if(i++>limit)
				break;
			System.out.println(entry.getKey()+":"+entry.getValue());
		}
	}
	
}
