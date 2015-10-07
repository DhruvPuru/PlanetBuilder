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

	// avoid instance creation
	public Orbit()
	{
		throw new UnsupportedOperationException();
	}

	// generate circular orbit from current position
	public Orbit(Point r)
	{
		// the radius is the same
		a = b = Math.sqrt(r.x * r.x + r.y * r.y);
		// rotation is undefined
		A = Double.NaN;
		// find the initial anomaly
		Mo = Math.atan2(r.y, r.x);
	}

	// generate elliptic orbit from current position and velocity
	public Orbit(Point r, Point v)
	{
		// find the eccentricity vector
		double r2 = r.x * r.x + r.y * r.y;
		double r1 = Math.sqrt(r2);
		double v2 = v.x * v.x + v.y * v.y;
		double rv = r.x * v.x + r.y * v.y;
		double f = v2 - GM / r1;
		double e_x = (f * r.x - rv * v.x) / GM;
		double e_y = (f * r.y - rv * v.y) / GM;
		double e2 = e_x * e_x + e_y * e_y;
		double e = Math.sqrt(e2);
		if (e == 1) throw new InvalidOrbitException("Parabolic orbit");
		if (e > 1) throw new InvalidOrbitException("Hyperbolic orbit");
		// compute the major radius from r and υ
		a = GM / (2.0 * GM / r1 - v2);
		// compute the minor radius from a and e
		b = a * Math.sqrt(1.0 - e2);
		// find the initial true anomaly θ
		double cos_theta = (e_x * r.x + e_y * r.y) / (e * r1);
		double theta = acos(cos_theta);
		if (rv < 0.0) theta = -theta;
		// set the rotation using θ
		A = Math.atan2(r.y, r.x) - theta;
		// find the initial eccentric anomaly
		double Eo = acos((e + cos_theta) / (1.0 + e * cos_theta));
		if (rv < 0.0) Eo = -Eo;
		// set the initial mean anomaly
		Mo = Eo - e * Math.sin(Eo);
	}

	// return the turn time
	public static double dt()
	{
		return dt;
	}

	// return the period of the orbit (quantized by dt)
	public long period()
	{
		double T = 2.0 * Math.PI * Math.sqrt(a) * a / Math.sqrt(GM);
		return (long) Math.ceil(T / dt);
	}

	// return the center of the orbit
	public void center(Point c)
	{
		// resolve circular case
		if (a == b) {
			c.x = c.y = 0.0;
			return;
		}
		// generic ellipse case
		double b_a = b / a;
		double e = Math.sqrt(1.0 - b_a * b_a);
		double len = a * e;
		double dir = A + Math.PI;
		c.x = len * Math.cos(dir);
		c.y = len * Math.sin(dir);
	}

	// find the position at time t
	public void positionAt(long t, Point r)
	{
		// compute the period and the time modulo
		long T_dt = period();
		t %= T_dt;
		// find the mean anomaly
		double M = Mo + 2.0 * Math.PI * t / T_dt;
		// resolve circular case
		if (a == b) {
			r.x = a * Math.cos(M);
			r.y = a * Math.sin(M);
			return;
		}
		// compute the eccentricity
		double b_a = b / a;
		double e = Math.sqrt(1.0 - b_a * b_a);
		// compute the eccentric anomaly using Newton-Raphson
		double E = M;
		for (;;) {
			// f(x) = x - e sin(x) - M
			double fE = E - e * Math.sin(E) - M;
			if (Math.abs(fE) < 1.0e-12) break;
			// f'(x) = 1 - e cos(x)
			E -= fE / (1.0 - e * Math.cos(E));
		}
		// center based coordinate system with foci on xx'
		r.x = a * Math.cos(E);
		r.y = b * Math.sin(E);
		// rotate the axes by A
		double len = Math.sqrt(r.x * r.x + r.y * r.y);
		double dir = Math.atan2(r.y, r.x) + A;
		r.x = len * Math.cos(dir);
		r.y = len * Math.sin(dir);
		// move the center of the axes
		len = a * e;
		dir = A + Math.PI;
		r.x += len * Math.cos(dir);
		r.y += len * Math.sin(dir);
	}

	// find the velocity at time t
	public void velocityAt(long t, Point v)
	{
		// find the position in orbit
		Point r = v;
		positionAt(t, r);
		// compute the eccentricity
		double b_a = b / a;
		double e = Math.sqrt(1.0 - b_a * b_a);
		// compute the angle between the foci using the law of cosines
		double r1 = Math.sqrt(r.x * r.x + r.y * r.y);
		double r_a = r1 / a;
		double phi = acos(cos(e + e, r_a, 2.0 - r_a));
		// compute the velocity direction assuming the foci are at xx'
		double dir = Math.atan2(r.y, r.x) + (Math.PI - phi) * 0.5;
		// compute the velocity magnitude using the vis-visa equation
		double len = Math.sqrt(GM * (2.0 / r1 - 1.0 / a));
		// compute the velocity vector
		v.x = len * Math.cos(dir);
		v.y = len * Math.sin(dir);
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
				return t + 1 != period();
			}

			public Point next()
			{
				if (t + 1 == period())
					return null;
				Point r = new Point();
				positionAt(++t, r);
				return r;
			}
		};
	}

	// compute arc cosine and include some error bound
	private static double acos(double cos)
	{
		if (Math.abs(cos) - 1.0 > 1.0e-9)
			throw new NumberFormatException("Invalid cosine");
		if (cos > +1.0) cos = +1.0;
		if (cos < -1.0) cos = -1.0;
		return Math.acos(cos);
	}

	// the law of cosines for angle opposite to a
	private static double cos(double a, double b, double c)
	{
		return (b * b - a * a + c * c) / (2.0 * b * c);
	}
}
