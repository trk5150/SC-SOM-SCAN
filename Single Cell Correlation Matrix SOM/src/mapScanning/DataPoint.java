package mapScanning;
public class DataPoint 
{
	public String look;
	public MiniNode myMini;
	public boolean count;
	public int x,y;
	public double mag;
	public double[] g;
	public String name;
	public String chrome;
	public int minLocus, maxLocus;
	
	public DataPoint()
	{
		count = true;
	}
	public DataPoint(String chr, int min, int max, MiniNode n)
	{
		chrome = chr;
		minLocus = min;
		maxLocus = max;
		myMini = n;
		count = true;
	}
	public DataPoint(String gene, MiniNode n)
	{
		name = gene;
		myMini = n;
	}
	public DataPoint(String s)
	{
		name = s;
		look = "";
	}
	public double updateMag()
	{
		double magD = 0;
		for(int i = 0; i<g.length; i++)
		{
			magD += (g[i])*(g[i]);
		}
		magD = Math.sqrt(magD);
		return magD;
	}
	public void gInit(int length) 
	{
		g = new double[length];
	}
	public void removeDuds(int newSize, int[] locs)
	{
		double[] gg = new double[newSize];
		int index = 0;
		for(int i = 0; i< g.length; i++)
		{
			for(int j = 0; j< locs.length; j ++)
			{
				gg[index] = g[i];
				if(i == locs[j])
					index--;
			}
			index++;
		}
		g = gg;
	}
	public void setName(String s)
	{
		name = s;
	}
}