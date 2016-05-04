package main;

import java.io.File;

public class main {

	public static void main(String[] args) {
		
		/**
		 * 文件名 不带.jpg
		 * 最小圆半径 9
		 * 腐蚀半径  2
		 * 灰度行程矩阵大小  96
		 * 种子区域生长时距离值 1.25
		 * 种子点半径阈值  20
		 */
		String picDir="H:/image/2/";
		String outDir="H:/ob1";
		File outFile =new File(outDir);
		if(!outFile.exists()){
			outFile.mkdir();
		}
		
		File picFile =new File(picDir);
		File [] pic=picFile.listFiles();
		
		for (int i=0;i<pic.length;i++){
			String picName=picDir+pic[i].getName().substring(0,pic[i].getName().length()-4);
			String outName=outDir+picDir+pic[i].getName().substring(0,pic[i].getName().length()-4);
			String []in={picName,"3","2","96","1.25","20"};
			Runlength.mainRun(in);
		}
		
		
		
		
	}

}
