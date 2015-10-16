package pb.sim;

import java.util.Iterator;

public class Orbit implements Iterable <Point> {

	// real time per quantum time unit (default is 1 Earth day)
	static double dt = 24 * 60 * 60;

	// standard gravitational parameter and mass of the Sun
	public static final double GM = 1.32712440018e20;
	public static final double  M = 1.98855e30;

	// the paramters of the ellipse
	public final double a;
	public final double b;
	public final double A;

	// the initial mean anomaly
	public final double Mo;

	// store extra parameters to save time
	private final double e;
	private final double c_x;
	private final double c_y;
	private final long T_dt;

	// avoid instance creation
	public Orbit()
	{
		throw new UnsupportedOperationException();
	}

	// generate circular orbit from current position
	public Orbit(Point r)
	{
		if (!r.finite()) throw new IllegalArgumentException("r not finite");
		// the radius is the same
		a = b = Math.sqrt(r.x * r.x + r.y * r.y);
		// rotation is undefined
		A = Double.NaN;
		// find the initial anomaly
		Mo = Math.atan2(r.y, r.x);
		// center and eccentricity is zero
		c_x = c_y = e = 0.0;
		// compute the period
		double T = (Math.PI + Math.PI) * Math.sqrt(a / GM) * a;
		T_dt = (long) Math.ceil(T / dt);
	}

	// generate elliptic orbit from current position and velocity
	public Orbit(Point r, Point v)
	{
		if (!r.finite()) throw new IllegalArgumentException("r not finite");
		if (!v.finite()) throw new IllegalArgumentException("v not finite");
		// compute the eccentricity e using r and v
		double r2 = r.x * r.x + r.y * r.y;
		double r1 = Math.sqrt(r2);
		double v2 = v.x * v.x + v.y * v.y;
		double rv = r.x * v.x + r.y * v.y;
		double f = v2 - GM / r1;
		double e_x = (f * r.x - rv * v.x) / GM;
		double e_y = (f * r.y - rv * v.y) / GM;
		double e2 = e_x * e_x + e_y * e_y;
		// rule out fly-by orbits
		if (e2 == 1.0) throw new InvalidOrbitException("Parabolic orbit");
		if (e2 > 1.0) throw new InvalidOrbitException("Hyperbolic orbit");
		// set very small eccentricities to zero
		if (e2 < 1.0e-12) e2 = 0.0;
		e = Math.sqrt(e2);
		// compute the major radius using r and υ
		a = GM / (2.0 * GM / r1 - v2);
		// compute the minor radius from a and e
		b = a * Math.sqrt(1.0 - e2);
		if (a == b) {
			// circular orbit
			if (e != 0.0) throw new ArithmeticException("e > 0 & a = b");
			// undefined rotation
			A = Double.NaN;
			// mean = true = eccentric anomaly
			Mo = Math.atan2(r.y, r.x);
			// center is zero
			c_x = c_y = 0.0;
		} else {
			if (e == 0.0) throw new ArithmeticException("e = 0 & a > b");
			// compute the initial true anomaly θ
			double cos_theta = (e_x * r.x + e_y * r.y) / (e * r1);
			double theta = acos(cos_theta);
			if (rv < 0.0) theta = -theta;
			// compute the ellipse rotation using θ and r
			A = Math.atan2(r.y, r.x) - theta;
			// compute the initial eccentric anomaly Eo using θ and e
			double Eo = acos((e + cos_theta) / (1.0 + e * cos_theta));
			if (Double.isNaN(Eo)) throw new ArithmeticException("Eo is NaN");
			if (rv < 0.0) Eo = -Eo;
			// compute the initial mean anomaly Mo using Eo and e
			Mo = Eo - e * Math.sin(Eo);
			if (Double.isNaN(Mo)) throw new ArithmeticException("Mo is NaN");
			// compute the center
			double c1 = a * e;
			double cD = A + Math.PI;
			c_x = c1 * Math.cos(cD);
			c_y = c1 * Math.sin(cD);
		}
		// compute the period
		double T = (Math.PI + Math.PI) * Math.sqrt(a / GM) * a;
		T_dt = (long) Math.ceil(T / dt);
	}

	// return the turn time
	public static double dt() { return dt; }

	// return the period of the orbit (quantized by dt)
	public long period() { return T_dt; }

	// return the center of the orbit
	public void center(Point c)
	{
		c.x = c_x;
		c.y = c_y;
	}

	// find the position at time t
	public void positionAt(long t, Point r)
	{
		// avoid modulo if possible
		if (t < 0) t %= T_dt;
		else if (t < T_dt) ;
		else if (t < T_dt + T_dt) t -= T_dt;
		else t %= T_dt;
		// mean anomaly
		double M = Mo + (Math.PI + Math.PI) * t / T_dt;
		// circular orbit case
		if (a == b) {
			r.x = a * Math.cos(M);
			r.y = a * Math.sin(M);
			return;
		}
		// compute the eccentric anomaly using Newton-Raphson
		double E = M;
		for (;;) {
			// f(x) = x - e sin(x) - M
			double fE = E - e * Math.sin(E) - M;
			if (Math.abs(fE) < 1.0e-12) break;
			// f'(x) = 1 - e cos(x)
			E -= fE / (1.0 - e * Math.cos(E));
			if (Double.isNaN(E)) throw new ArithmeticException("E is NaN");
		}
		// center based coordinate system with foci on xx'
		r.x = a * Math.cos(E);
		r.y = b * Math.sin(E);
		// rotate the axes by A and move the center of the axes
		double r1 = Math.sqrt(r.x * r.x + r.y * r.y);
		double rD = Math.atan2(r.y, r.x) + A;
		r.x = c_x + r1 * Math.cos(rD);
		r.y = c_y + r1 * Math.sin(rD);
	}

	// find the velocity at time t
	public void velocityAt(long t, Point v)
	{
		// find the position in orbit (use same Point object)
		Point r = v;
		positionAt(t, r);
		// compute the angle between the foci using the law of cosines
		double r1 = Math.sqrt(r.x * r.x + r.y * r.y);
		double r_a = r1 / a;
		double phi = acos(cos(e + e, r_a, 2.0 - r_a));
		// compute the velocity direction assuming the foci are at xx'
		double vD = Math.atan2(r.y, r.x) + (Math.PI - phi) * 0.5;
		// compute the velocity magnitude using the vis-visa equation
		double v1 = Math.sqrt((GM + GM) / r1 - GM / a);
		// compute the velocity vector
		v.x = v1 * Math.cos(vD);
		v.y = v1 * Math.sin(vD);
	}

	// allocate Point and call center()
	public Point center()
	{
		Point p = new Point();
		center(p);
		return p;
	}

	// allocate Point and call positionAt()
	public Point positionAt(long t)
	{
		Point p = new Point();
		positionAt(t, p);
		return p;
	}

	// allocate Point and call velocityAt()
	public Point velocityAt(long t)
	{
		Point p = new Point();
		velocityAt(t, p);
		return p;
	}

	// iterator over orbital positions
	public Iterator <Point> iterator()
	{
		return new Iterator <Point> () {

			private long t = -1;

			public boolean hasNext()
			{
				return t + 1 != T_dt;
			}

			public Point next()
			{
				if (t + 1 == T_dt) return null;
				Point r = new Point();
				positionAt(++t, r);
				return r;
			}
		};
	}

	// compute arc cosine and include some error bound
	private static double acos(double cos)
	{
		if (Math.abs(cos) - 1.0 > 1.0e-12)
			throw new ArithmeticException("|cos(x)| > 1");
		if (cos > +1.0) cos = +1.0;
		if (cos < -1.0) cos = -1.0;
		return Math.acos(cos);
	}

	// the law of cosines for angle opposite to a
	private static double cos(double a, double b, double c)
	{
		double bc = b * c;
		return (b * b - a * a + c * c) / (bc + bc);
	}
}
