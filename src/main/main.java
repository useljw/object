package main;

import java.io.File;

public class main {

	public static void main(String[] args) {
		
		/**
		 * �ļ��� ����.jpg
		 * ��СԲ�뾶 9
		 * ��ʴ�뾶  2
		 * �Ҷ��г̾����С  96
		 * ������������ʱ����ֵ 1.25
		 * ���ӵ�뾶��ֵ  20
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
