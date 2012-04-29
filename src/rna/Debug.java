package rna;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Debug {
	private static boolean DEBUG = false;
	private static boolean toFile = false;
	private static PrintWriter pr ;
	public static void log(String a){
		if(DEBUG){
			if(toFile){
				try {
					pr = new PrintWriter(new FileWriter("LOG"));
					pr.println(a);
					pr.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			System.out.println(a);
		}
	}
	public static void log(double a){
		log(a+"");
	}
	public static void log(int a){
		log(a+"");
	}
	public static void log(boolean a){
		log(a+"");
	}
	public static void log(Object a){
		log(a.toString());
	}
	public static void log(String format, Object ...objects ){
		log(String.format(format, objects));
	}
	public static void print(int[][]a){
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
	public static void print(String[][] a) {
		if(!DEBUG)return;
		System.out.println("~~~~~");
		for(int i =0;i<a.length;i++){
			for(int j =0;j<a[0].length;j++){
				//System.out.print("start:"+i+","+"end:"+j+"-"+a[i][j]+"  ");
				System.out.print(a[i][j]+"\t\t");
			}
			System.out.println();
		}
	}
	public static void print(double[][] a) {
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
}
