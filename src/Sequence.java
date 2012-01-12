
public class Sequence {
	private final boolean DEBUG = true;
	private String basepairs;
	private String foldedStruct;
	private int maxPairings;
	private int [][] completedMatrix = null;
	private String [][] completeStructMatrix = null;
	public String getBasepairs() {
		return basepairs;
	}

	public void setBasepairs(String basepairs) {
		this.basepairs = basepairs;
		this.foldedStruct = "";
		this.maxPairings = -1;
		this.completedMatrix = null;
		this.completeStructMatrix = null;
	}

	public Sequence(String basepairs) {
		this.basepairs = basepairs.toLowerCase();
	}
	public String getFoldedStruct(){
		if(this.completedMatrix==null)
			fold();
		return foldedStruct;
	}
	public int getMaxPairings(){
		if(this.completedMatrix==null)
			fold();
		return maxPairings;
	}
	public void fold(){
		getPairingMatrix();
		maxPairings = completedMatrix[0][completedMatrix.length-1];
		foldedStruct = completeStructMatrix[0][completeStructMatrix.length-1];
	}
	private void print(int[][]a){
		if(!DEBUG)return;
		System.out.println("~~~~~");
		for(int i =0;i<a.length;i++){
			for(int j =0;j<a[0].length;j++){
				//System.out.print("start:"+i+","+"end:"+j+"-"+a[i][j]+"  ");
				System.out.print(a[i][j]+" ");
			}
			System.out.println();
		}
	}

	private void getPairingMatrix(){
		int [][] m = new int [basepairs.length()][basepairs.length()];
		String [][] optimalStructs = new String [basepairs.length()][basepairs.length()];
		//initialization
		for(int i = 0;i<m.length;i++){
			m[i][i]=0;
			optimalStructs[i][i]=".";
			if(i>0){
				m[i][i-1]=0;
				optimalStructs[i][i-1]="";
			}
			if(i<m.length-1){
				if(isValidPair(i,i+1)){
					m[i][i+1]=1;
					optimalStructs[i][i+1]="()";
				}
			}
		}
		for(int end=1;end<m.length;end++){
			for(int start =end-1;start>-1;start--){
				int max = 0;
				
				//worst case is all '.' for the entire sequence, no basepairing
				String bestStruct = new String(new char[end-start+1]).replace("\0", ".");;
				if(isValidPair(start,end)){
					if(m[start+1][end-1]+1>max)
						max = m[start+1][end-1]+1;
						bestStruct = "("+optimalStructs[start+1][end-1]+")";
				}

				if(m[start+1][end]>max){
					max = m[start+1][end];
					bestStruct = "."+optimalStructs[start+1][end];
				}
				if(m[start][end-1]>max){
					max = m[start][end-1];
					bestStruct = optimalStructs[start][end-1]+".";
				}
				for(int i =start+1;i<end;i++){
					int curBpCount = m[start][i]+m[i+1][end];
					if(curBpCount>max){
						max = curBpCount;
						bestStruct = optimalStructs[start][i]+optimalStructs[i+1][end];
					}
				}

				m[start][end]=max;
				optimalStructs[start][end] = bestStruct;
			}
		}
		this.completedMatrix=m;
		this.completeStructMatrix=optimalStructs;
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

}
