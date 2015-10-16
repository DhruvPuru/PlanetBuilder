package pb.g6;

import pb.sim.Point;
//import pb.sim.Orbit;
import pb.sim.Asteroid;
//import pb.sim.InvalidOrbitException;

import java.util.ArrayList;

/**
 * Example usage:
 * Collision collision = new Collision(asteroids);
 * for allKindsOfPush:
 * collision.updateAsteroid(0, pushedAsteroid);
 * nearestTime = collision.findCollisionTime(asteroidIndex1, asteroidIndex2);
 */
public class Collision {
    public final long TIME_LIMIT = 1000 * 24 * 60 * 60;
    public Asteroid[] asteroids;
    public int curPos;  // position of the current push
    public Asteroid oriAsteroid;
    public final double DELTA = 10;

    public Collision(Asteroid[] asteroids) {
        this.asteroids = asteroids;  // is shallow copy good enough?
        oriAsteroid = null;
    }

    public void updateAsteroid(int pushedAsteroidPos, Asteroid pushedAsteroid) {
        if (oriAsteroid != null) {
            // reverse the previous push
            asteroids[curPos] = oriAsteroid;
        }
        assert pushedAsteroidPos < asteroids.length : "Failure message";
        curPos = pushedAsteroidPos;
        oriAsteroid = asteroids[curPos];
        asteroids[curPos] = pushedAsteroid;
    }

    public Point[] getFoci(double a, double b, double A, double h, double k) {
        double c = 0.0;
        if (a == b) {
            return new Point[]{new Point(c, c)};
        } else if (a > b) {
            c = Math.sqrt(Math.pow(a, 2) - Math.pow(b, 2));
            Point f1 = new Point(c, 0.0);
            Point f2 = new Point(-c, 0.0);

            f1 = new Point(f1.x * Math.cos(A) - f1.y * Math.sin(A) + h, f1.x * Math.sin(A) + f1.y * Math.cos(A) + k);
            f2 = new Point(f2.x * Math.cos(A) - f2.y * Math.sin(A) + h, f2.x * Math.sin(A) + f2.y * Math.cos(A) + k);

            return new Point[]{f1, f2};
        } else {
            c = Math.sqrt(Math.pow(b, 2) - Math.pow(a, 2));
            Point f1 = new Point(0.0, c);
            Point f2 = new Point(0.0, -c);

            f1 = new Point(f1.x * Math.cos(A) - f1.y * Math.sin(A) + h, f1.x * Math.sin(A) + f1.y * Math.cos(A) + k);
            f2 = new Point(f2.x * Math.cos(A) - f2.y * Math.sin(A) + h, f2.x * Math.sin(A) + f2.y * Math.cos(A) + k);

            return new Point[]{f1, f2};
        }
    }

    public Point findPoint(double a, double b, double A, double h, double k, double t) {
        t = (t * Math.PI) / 180;
        return new Point(a * Math.cos(t) * Math.cos(A) - b * Math.sin(t) * Math.sin(A) + h, a * Math.cos(t) * Math.sin(A) + b * Math.sin(t) * Math.cos(A) + k);
    }

    public ArrayList<ArrayList<Point>> findIntersection(Asteroid a, int ignore_idx) {
        ArrayList<ArrayList<Point>> intersect = new ArrayList<ArrayList<Point>>(asteroids.length);

        Ellipse[] ellipses = new Ellipse[asteroids.length];

        double a1 = a.orbit.a;
        double b1 = a.orbit.b;
        double A1 = a.orbit.A;
        Point c1 = a.orbit.center();

        for (int i = 0; i < asteroids.length; i++) {
            ArrayList<Point> inner = new ArrayList<Point>();
            intersect.add(inner);

            double a2 = asteroids[i].orbit.a;
            double b2 = asteroids[i].orbit.b;
            double A2 = asteroids[i].orbit.A;
            Point c2 = asteroids[i].orbit.center();
            Point[] f = getFoci(a2, b2, A2, c2.x, c2.y);

            ellipses[i] = new Ellipse(a2, b2, A2, c2.x, c2.y, f);
        }


        Point p = null;
        for (double angle = 0; angle < 360; angle += 0.2d) {
            p = findPoint(a1, b1, A1, c1.x, c1.y, angle);
            for (int i = 0; i < asteroids.length; i++) {
                if (i == ignore_idx) {
                    continue;
                }

                double sum = 0;
                for (Point f : ellipses[i].foci) {
                    sum += Point.distance(f, p);
                }

                double radius = asteroids[i].radius();
                if (ellipses[i].foci.length == 0) {
                    throw new Error();
                } else if (ellipses[i].foci.length == 1) {
                    if (sum > ellipses[i].a - radius && sum < ellipses[i].a + radius) {
                        intersect.get(i).add(p);
                    }
                } else {
                    if (sum > 2 * ellipses[i].a - radius && sum < 2 * ellipses[i].a + radius) {
                        intersect.get(i).add(p);
                    }
                }
            }
        }

        return intersect;
    }

    /**
     * Find time of asteroid at position p
     */
    public long timeAt(Asteroid a, Point p) {
        long res;
        long period = a.orbit.period();
        Point curPos = new Point();
        a.orbit.positionAt(0, curPos);
        double curAngle = Math.atan2(curPos.y, curPos.x);
        double lastAngle;
        double pAngle = Math.atan2(p.y, p.x);

        long dt = period / 2;
        long t = 0;
        System.out.println("==================================================");
        while (Point.distance(curPos, p) > a.radius()) {
            // System.out.println("pAngle " + pAngle);
            lastAngle = curAngle;
            System.out.println("try " + dt);
            System.out.println("dist " + Point.distance(curPos, p));
            // jump
            t += dt;
            a.orbit.positionAt(t, curPos);
            curAngle = Math.atan2(curPos.y, curPos.x);
            // System.out.println("try " + curPos);

            if ((lastAngle*curAngle >= 0 && (curAngle-pAngle)*(lastAngle-pAngle) <= 0)
                    || lastAngle*curAngle < 0 && (curAngle-pAngle)*(lastAngle-pAngle) >= 0) {
                if (dt == 1 || dt == 0) {
                    System.err.println("timeAt error");
                }
                // if pass by the point, get the last time
                t -= dt; // come back to before jumping
                dt /= 2;
            }
        }
        // System.err.println("timeAt found");
        return t;
    }

    /**
     * Find collision time when push asteroids[ignore_idx] to the state of asteroid a
     * @return list of all collision time under TIME_LIMIT.
     */
    public ArrayList<Long> findCollisionTime(Asteroid a, int ignore_idx) {
        ArrayList<ArrayList<Point>> intersect = findIntersection(a, ignore_idx);
        ArrayList<Long> collision = new ArrayList<Long>();
        long t1, t2, k1, k2, period1, period2, colTime;

        for (int i = 0; i < intersect.size(); i++) {
            if (i == ignore_idx) {
                continue;
            }

            // find all collisions between a and asteroids[i]
            period1 = a.orbit.period();
            period2 = asteroids[i].orbit.period();
            for (int j = 0; j < intersect.get(i).size(); j++) {
                Point p = intersect.get(i).get(j);
                t1 = timeAt(a, p);
                t2 = timeAt(asteroids[i], p);
                System.out.println(t1 + " " + t2);
                k2 = 1;
                // find int values of k1 and k2 such that t1 + k1*period1 = t2 + k2*period2
                do {
                    k1 = (k2*period2 + t2 - t1) / (long)period1;
                    ++k2;
                } while (k2 < (TIME_LIMIT/period2) && Math.round(k1) != k1);

                colTime = t1 + k1*period1;
                System.err.println("Collision between " + ignore_idx + " & " + i);
                System.err.println("  Year: " + (1 + colTime / 365));
                System.err.println("  Day: " + (1 + colTime % 365));
                collision.add(colTime);
            }
        }

        return collision;
    }
}

