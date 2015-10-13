package pb.sim;

import java.util.*;

public class Point implements Comparable <Point> {

	public double x;
	public double y;

	public Point(double x, double y)
	{
		this.x = x;
		this.y = y;
	}

	public Point()
	{
		x = y = 0.0;
	}

	// length of vector
	public double magnitude()
	{
		return Math.sqrt(x * x + y * y);
	}

	// direction of vector
	public double direction()
	{
		return Math.atan2(y, x);
	}

	// Euclidean distance between two points
	public static double distance(Point p1, Point p2)
	{
		double dx = p1.x - p2.x;
		double dy = p1.y - p2.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	// find which (sets of) overlaps (null circles are ignored)
	public static int[][] overlaps(Point[] center, double[] radius)
	{
		if (center.length != radius.length)
			throw new IllegalArgumentException();
		// create a graph of overlapping circles (ignore if isolated)
		HashMap <Integer, HashSet <Integer>> E = null;
		for (int i = 0 ; i != center.length ; ++i) {
			if (center[i] == null) continue;
			if (Double.isNaN(center[i].x))
				throw new IllegalArgumentException("Center x is not a number");
			if (Double.isInfinite(center[i].x))
				throw new IllegalArgumentException("Center x is infinite");
			if (Double.isNaN(center[i].y))
				throw new IllegalArgumentException("Center y is not a number");
			if (Double.isInfinite(center[i].y))
				throw new IllegalArgumentException("Center y is infinite");
			if (Double.isNaN(radius[i]))
				throw new IllegalArgumentException("Radius is not a number");
			if (Double.isInfinite(radius[i]))
				throw new IllegalArgumentException("Radius is infinite");
			for (int j = 0 ; j != i ; ++j) {
				if (center[j] == null) continue;
				double R = radius[i] + radius[j];
				double dx = center[i].x - center[j].x;
				double dy = center[i].y - center[j].y;
				if (dx * dx + dy * dy < R * R) {
					if (E == null)
						E = new HashMap <Integer, HashSet <Integer>>();
					if (E.get(i) == null)
						E.put(i, new HashSet <Integer> ());
					if (E.get(j) == null)
						E.put(j, new HashSet <Integer> ());
					E.get(i).add(j);
					E.get(j).add(i);
				}
			}
		}
		if (E == null) return none;
		// split graph in connected components
		ArrayList <int[]> Cs = new ArrayList <int[]> ();
		HashSet <Integer> C = new HashSet <Integer> ();
		Stack <Integer> S = new Stack <Integer> ();
		do {
			int i0 = E.keySet().iterator().next();
			S.add(i0);
			C.add(i0);
			do {
				int i = S.pop();
				for (int j : E.get(i))
					if (C.add(j)) S.push(j);
				E.remove(i);
			} while (!S.empty());
			Cs.add(toInt(C));
			C.clear();
		} while (!E.isEmpty());
		return Cs.toArray(none);
	}

	// keep a zero array to avoid reallocating it
	private static final int[][] none = new int[0][];

	// array of integers from set (internal function)
	private static int[] toInt(Collection <Integer> set)
	{
		int[] array = new int [set.size()];
		int i = 0;
		for (int value : set)
			array[i++] = value;
		return array;
	}

	// point equality using bit representation
	public boolean equals(Object o)
	{
		if (o instanceof Point) {
			Point p = (Point) o;
			return Double.doubleToLongBits(x) ==
			       Double.doubleToLongBits(p.x) &&
			       Double.doubleToLongBits(y) ==
			       Double.doubleToLongBits(p.y);
		}
		return false;
	}

	// print point
	public String toString()
	{
		return "(" + x + ", " + y + ")";
	}

	// hash points
	public int hashCode()
	{
		double[] xy = new double [2];
		xy[0] = x;
		xy[1] = y;
		return Arrays.hashCode(xy);
	}

	// order points lexicographically
	public int compareTo(Point p)
	{
		if (x < p.x) return -1;
		if (x > p.x) return  1;
		if (y < p.y) return -1;
		if (y > p.y) return  1;
 		return 0;
	}
}
