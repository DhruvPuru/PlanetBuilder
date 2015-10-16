package pb.g2;

public class ComparableAsteroid implements Comparable<ComparableAsteroid> {
    public int index;
    public double radius;
    public double mass;

    public ComparableAsteroid(int index, double radius, double mass) {
        this.index = index;
        this.radius = radius;
        this.mass = mass;
    }

    private double getScore() {
        return radius*radius*mass;
    }

    public int compareTo(ComparableAsteroid other) {
        double thisScore = this.getScore();
        double otherScore = other.getScore();
        if (thisScore > otherScore) {
            return -1;
        } else if (thisScore < otherScore) {
            return 1;
        } else {
            return 0;
        }
    }
}
