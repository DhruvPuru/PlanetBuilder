package pb.g7;

public class asteroid_index implements Comparable<asteroid_index> {

    public int index;
    public double r;
    
    public asteroid_index(int index, double r){
	this.index = index;
	this.r = r;
    }

    public int compareTo(asteroid_index other) {
	/*	if (Math.max(a,b) < Math.max(other.a, other.b)) {return -1;}
		else if(Math.max(a,b) > Math.max(other.a, other.b)) {return 1;}*/
	if (r > other.r) {return 1;}
	else if (r < other.r) {return -1;}
	else {return 0;}
    }

    public double getRadius() {
	return r;
    }
}
