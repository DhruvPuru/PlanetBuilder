package pb.g6;

import pb.sim.Point;

public class Ellipse {
    public double a, b, A, h, k;
    public Point[] foci;

    public Ellipse() {
        a = 0.0;
        b = 0.0;
        A = 0.0;
        h = 0.0;
        k = 0.0;
        foci = new Point[] {};
    }

    public Ellipse(double a, double b, double A, double h, double k, Point[] f) {
        this.a = a;
        this.b = b;
        this.A = A;
        this.h = h;
        this.k = k;
        this.foci = f;
    }
}
