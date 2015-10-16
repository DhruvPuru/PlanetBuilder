package pb.sim;

import java.io.*;
import java.net.*;
import java.util.*;
import javax.tools.*;
import java.awt.Desktop;
import java.util.concurrent.*;

class Simulator {

	// root directory
	private static final String root = "pb";

	// main function (change default parameters here)
	public static void main(String[] args)
	{
		// default number of asteroids without state file
		int default_asteroids = 5;
		int random_asteroids = -1;
		// default game_time limit (in days) is one millenium
		long game_time_limit = 365 * 1000;
		// GUI parameters
		boolean gui = false;
		boolean gui_fast_forward = false;
		boolean gui_planets = true;
		int     gui_refresh_rate = 20;
		// state file path
		String state_file = null;
		// min and max orbit
		double min_orbit = 250.0e9;
		double max_orbit = 700.0e9;
		// min and max mass
		double min_mass = 0.8 * Orbit.M;
		double max_mass = 1.2 * Orbit.M;
		// initial state
		Point[] asteroid_position = null;
		double[] asteroid_mass = null;
		// CPU timeouts in ms
		long init_timeout = 0;
		long play_timeout = 0;
		long total_timeout = 0;
		// the player
		Class <Player> player = null;
		String group = "g0";
		try {
			for (int a = 0 ; a != args.length ; ++a)
				if (args[a].equals("-g") || args[a].equals("--group")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing group name");
					group = args[a];
				} else if (args[a].equals("-s") || args[a].equals("--state")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing state file");
					state_file = args[a];
				} else if (args[a].equals("-a") || args[a].equals("--asteroids")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing state file");
					random_asteroids = Integer.parseInt(args[a]);
					if (random_asteroids < 2)
						throw new IllegalArgumentException("Cannot have less than 2 asteroids");
				} else if (args[a].equals("-t") || args[a].equals("--game_time")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing game_time limit parameter");
					game_time_limit = Long.parseLong(args[a]);
					if (game_time_limit < 0)
						throw new IllegalArgumentException("Invalid game_time limit");
				} else if (args[a].equals("--gui-fps")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing the FPS parameter");
					double fps = Double.parseDouble(args[a]);
					if (fps < 0.0)
						throw new IllegalArgumentException("Invalid FPS parameter");
					if (fps == 0.0)
						gui_refresh_rate = -1;
					else
						gui_refresh_rate = (int) Math.round(1000 / fps);
					gui = true;
				} else if (args[a].equals("--init-timeout")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing the init() timeout parameter");
					init_timeout = Long.parseLong(args[a]);
				} else if (args[a].equals("--play-timeout")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing the play() timeout parameter");
					play_timeout = Long.parseLong(args[a]);
				} else if (args[a].equals("--total-timeout")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing the total timeout parameter");
					total_timeout = Long.parseLong(args[a]);
				} else if (args[a].equals("--orbit-range")) {
					if (a + 2 >= args.length)
						throw new IllegalArgumentException("Missing orbit range parameters");
					min_orbit = Integer.parseInt(args[++a]) * 1.0e9;
					max_orbit = Integer.parseInt(args[++a]) * 1.0e9;
					if (min_orbit <= 0 || max_orbit <= 0 || min_orbit > max_orbit)
						throw new IllegalArgumentException("Invalid orbit range parameters");
				} else if (args[a].equals("--mass-range")) {
					if (a + 2 >= args.length)
						throw new IllegalArgumentException("Missing orbit range parameters");
					min_mass = Double.parseDouble(args[++a]) * Orbit.M;
					max_mass = Double.parseDouble(args[++a]) * Orbit.M;
					if (min_mass <= 0 || max_mass <= 0 || min_mass > max_mass)
						throw new IllegalArgumentException("Invalid orbit range");
				} else if (args[a].equals("--density")) {
					if (++a == args.length)
						throw new IllegalArgumentException("Missing density parameter");
					Asteroid.density *= Double.parseDouble(args[a]);
					if (Asteroid.density <= 0.0)
						throw new IllegalArgumentException("Invalid density parameter");
				} else if (args[a].equals("--gui")) {
					gui = true;
				} else if (args[a].equals("--gui-no-planets")) {
					gui_planets = false;
					gui = true;
				} else if (args[a].equals("--gui-fast-forward")) {
					gui_fast_forward = true;
					gui = true;
				} else throw new IllegalArgumentException("Unknown argument: " + args[a]);
			// load player
			player = load(group);
			// figure out where to get input from
			if (random_asteroids < 0 && state_file != null) {
				// load state from file
				double[][] asteroid_mass_addr = new double [1][];
				asteroid_position = load(state_file, asteroid_mass_addr);
				asteroid_mass = asteroid_mass_addr[0];
				int asteroids = asteroid_position.length;
				System.err.println("Loaded asteroids from state file");
			} else {
				// generate random orbits
				int asteroids = random_asteroids < 0 ? default_asteroids
				                                     :  random_asteroids;
				asteroid_mass = random_masses(asteroids, min_mass, max_mass);
				asteroid_position = random_orbits(asteroid_mass, min_orbit,
				                                                 max_orbit);
				System.err.println("Generated random asteroids");
				// store to file
				if (state_file != null) {
					store(state_file, asteroid_position, asteroid_mass);
					System.err.println("Stored to state file");
				}
			}
		} catch (Exception e) {
			System.err.println("Error during setup: " + e.getMessage());
			e.printStackTrace();
			System.err.println("Exiting the simulator ...");
			System.exit(1);
		}
		// print other parameters
		System.err.println("Asteroids: " + asteroid_position.length);
		System.err.println("Group: " + group);
		System.err.println("Game time limit (in dt units): " +
		                   (game_time_limit < 0 ? "+oo" : game_time_limit));
		System.err.println("Timeout for init() in ms: " +
		                   (init_timeout > 0 ? init_timeout : "+oo"));
		System.err.println("Timeout for play() in ms: " +
		                   (play_timeout > 0 ? play_timeout : "+oo"));
		System.err.println("Total timeout in ms: " +
		                   (total_timeout > 0 ? total_timeout : "+oo"));
		// print the ranges
		min_orbit = max_orbit = asteroid_position[0].magnitude();
		min_mass = max_mass = asteroid_mass[0];
		for (int i = 1 ; i != asteroid_position.length ; ++i) {
			double mass = asteroid_mass[i];
			double orbit = asteroid_position[i].magnitude();
			if (min_mass > mass) min_mass = mass;
			if (max_mass < mass) max_mass = mass;
			if (min_orbit > orbit) min_orbit = orbit;
			if (max_orbit < orbit) max_orbit = orbit;
		}
		System.err.println("Mass  range:  [" + min_mass + ", " +
		                                       max_mass + "]");
		System.err.println("Orbit range:  [" + min_orbit + ", " +
		                                       max_orbit + "]");
		// print GUI parameters
		if (!gui) System.err.println("GUI: disabled");
		else {
			if (gui_refresh_rate < 0)
				System.err.println("GUI: 0 FPS  [reload manually]");
			else if (gui_refresh_rate == 0)
				System.err.println("GUI: max FPS");
			else {
				double fps = (100 * 1000 / gui_refresh_rate) / 100.0;
				System.err.println("GUI: up to " + (int) fps + " FPS");
			}
			System.err.println("GUI fast forward: " +
			                   (gui_fast_forward ? "yes" : "no"));
			System.err.println("GUI planets: " +
			                   (gui_planets ? "yes" : "no"));
		}
		Info i = null;
		try {
			i = game(group, player, asteroid_position, asteroid_mass,
			         gui, gui_fast_forward, gui_planets, gui_refresh_rate,
			         game_time_limit, init_timeout, play_timeout, total_timeout);
		} catch (Exception e) {
			System.err.println("Error during play: " + e.getMessage());
			e.printStackTrace();
			System.err.println("Exiting the simulator ...");
			System.exit(1);
		}
		if (i == null)
			System.err.println("A game error occured during the simulation");
		else {
			System.err.println("Game game_time: " + i.game_time + " \"days\"");
			System.err.println("CPU game_time: " + i.cpu_time / 1.0e9 + " seconds");
			if (i.max_mass > i.total_mass * 0.5)
				System.err.println("Planet built with " + i.sum_energy + " Joules");
			else {
				double m_p = i.max_mass * 100.0 / i.total_mass;
				System.err.println("Failed to build the planet in game_time");
				System.err.println("Largest asteroid has " + m_p + "% of mass");
			}
		}
		System.exit(0);
	}

	// game result
	private static class Info {

		public final double sum_energy;
		public final double max_mass;
		public final double total_mass;
		public final long game_time;
		public final long cpu_time;

		public Info(double s_e, double m_m, double t_m, long g_t, long c_t)
		{
			sum_energy = s_e;
			max_mass = m_m;
			total_mass = t_m;
			game_time = g_t;
			cpu_time = c_t;
		}
	}

	// play the game
	private static Info game(String player_name,
	                         Class <Player> player_class,
	                         Point[]  asteroid_position,
	                         double[] asteroid_mass,
	                         boolean gui,
	                         boolean gui_fast_forward,
	                         boolean gui_planets,
	                         int gui_refresh_rate,
	                         long game_time_limit,
	                         long init_timeout,
	                         long play_timeout,
	                         long total_timeout) throws Exception
	{
		// generate the circular orbits of asteroids
		int n_asteroids = asteroid_position.length;
		if (n_asteroids != asteroid_mass.length)
			throw new IllegalArgumentException();
		Asteroid[] asteroids = new Asteroid [n_asteroids];
		double total_mass = 0.0;
		for (int i = 0 ; i != n_asteroids ; ++i) {
			Orbit orbit = new Orbit(asteroid_position[i]);
			double mass = asteroid_mass[i];
			asteroids[i] = new Asteroid(orbit, mass, 0);
			total_mass += mass;
		}
		// start the timer
		Timer timer = new Timer();
		timer.start();
		// initialize the player
		Asteroid[] asteroids_copy = Arrays.copyOf(asteroids,
		                                          asteroids.length);
		final Asteroid[] a0_final = asteroids_copy;
		final Class <Player> c_final = player_class;
		final long t_final = game_time_limit;
		final Player player = timer.call(new Callable <Player> () {

			public Player call() throws Exception
			{
				Player p = c_final.newInstance();
				p.init(a0_final, t_final);
				return p;
			}
		}, init_timeout);
		// push information
		double sum_energy = 0.0;
		List <Push> pushes = new ArrayList <Push> ();
		// initialize the GUI
		Orbit[] planets = null;
		HTTPServer server = null;
		long gui_time = 0;
		if (gui) {
			// initialize server and print the port
			server = new HTTPServer();
			int port = server.port();
			System.err.println("HTTP port: " + port);
			// try to open web browser automatically
			if (!Desktop.isDesktopSupported())
				System.err.println("Desktop operations not supported");
			else {
				Desktop desktop = Desktop.getDesktop();
				if (!desktop.isSupported(Desktop.Action.BROWSE))
					System.err.println("Desktop browsing not supported");
				else {
					URI uri = new URI("http://localhost:" + port);
					desktop.browse(uri);
				}
			}
			// get orbits of known planets and show their periods
			planets = gui_planets ? planet_orbits() : new Orbit[0];
			// send initial state until successful
			int refresh = game_time_limit == 0 ? -1 : gui_refresh_rate;
			long cpu_time = timer.cpu_ns();
			String content = state(player_name, planets, asteroids,
			                       pushes, sum_energy, 0, cpu_time, refresh);
			gui(server, content);
			gui_time = System.currentTimeMillis();
		}
		// start playing the game
		double[] energy = new double [n_asteroids];
		double[] direction = new double [n_asteroids];
		boolean game_over = false;
		long game_time_of_last_event = 0;
		long game_time = -1;
		// check game time and game termination conditions
		while (!game_over && ++game_time != game_time_limit) {
			// check total CPU time timeout
			if (total_timeout > 0) {
				long cpu_ms = timer.cpu_ns() / 1000000;
				if (cpu_ms > total_timeout) break;
			}
			// reset the energy and direction and copy the asteroids
			for (int i = 0 ; i != asteroids.length ; ++i) {
				asteroids_copy[i] = asteroids[i];
				energy[i] = 0.0;
				direction[i] = Double.NaN;
			}
			// call the play() method of the player
			final Asteroid[] a_final = asteroids_copy;
			final double[] e_final = energy;
			final double[] d_final = direction;
			timer.call(new Callable <Object> () {

				public Object call()
				{
					player.play(a_final, e_final, d_final);
					return null;
				}
			}, play_timeout);
			// check for pushes from the player
			for (int i = 0 ; i != asteroids.length ; ++i) {
				// validate energy
				if (Double.isNaN(energy[i]) ||
				    Double.isInfinite(energy[i]))
					throw new IllegalArgumentException("Undefined energy");
				if (energy[i] < 0)
					throw new IllegalArgumentException("Negative energy");
				// skip if no energy
				if (energy[i] == 0) continue;
				// validate direction
				if (Double.isNaN(direction[i]) ||
				    Double.isInfinite(direction[i])) throw new
					IllegalArgumentException("Undefined direction");
				// attempt to push the asteroid
				try {
					asteroids[i] = Asteroid.push(asteroids[i], game_time,
					                             energy[i], direction[i]);
					sum_energy += energy[i];
					game_time_of_last_event = game_time;
					pushes.add(new Push(energy[i], game_time));
				} catch (InvalidOrbitException e) {
					System.err.println("Push failed: " + e.getMessage());
				}
			}
			// check for collisions of asteroids
			int[][] index = Asteroid.test_collision(asteroids, game_time);
			// merge colliding asteroids
			for (int i = 0 ; i != index.length ; ++i) {
				// build array of asteroids that collide
				Asteroid[] as = new Asteroid [index[i].length];
				for (int j = 0 ; j != index[i].length ; ++j)
					as[j] = asteroids[index[i][j]];
				// merge the asteroids in a single asteroid
				System.err.println(game_time + ": Collision of "
				                   + as.length + " asteroids");
				Asteroid a = null;
				try {
					a = Asteroid.force_collision(as, game_time);
				} catch (InvalidOrbitException e) {
					System.err.println("Invalid collision orbit: "
					                   + e.getMessage());
					return null;
				}
				// stop the game if the asteroid mass is more than half
				if (a.mass > total_mass * 0.5) game_over = true;
				// store merged asteroids
				asteroids[index[i][0]] = a;
				for (int j = 1 ; j != index[i].length ; ++j)
					asteroids[index[i][j]] = null;
			}
			// compact asteroids after merging
			if (index.length != 0) {
				asteroids = compact(asteroids);
				game_time_of_last_event = game_time;
				System.err.println("Remaining: " + asteroids.length);
				energy    = new double [asteroids.length];
				direction = new double [asteroids.length];
				asteroids_copy = new Asteroid [asteroids.length];
			}
			// process GUI refresh
			long turns = game_time - game_time_of_last_event;
			if (game_over || game_time + 1 == game_time_limit) turns = -1;
			else if (!gui_fast_forward) turns = 0;
			else {
				long t = System.currentTimeMillis();
				if (t - gui_time > 100) turns = 0;
			}
			if (gui && turns < 50) {
				int refresh = turns < 0 ? -1 : gui_refresh_rate;
				String content = state(player_name, planets, asteroids,
				                       pushes, sum_energy, game_time,
				                       timer.cpu_ns(), refresh);
				gui(server, content);
				gui_time = System.currentTimeMillis();
			}
		}
		// find the mass of the largest asteroid
		double max_mass = 0.0;
		for (int i = 0 ; i != asteroids.length ; ++i)
			if (max_mass < asteroids[i].mass)
				max_mass = asteroids[i].mass;
		// return all info
		return new Info(sum_energy, max_mass, total_mass,
		                game_time, timer.cpu_ns());
	}

	// remove null Asteroid objects and compact
	private static Asteroid[] compact(Asteroid[] a)
	{
		int j;
		for (int i = j = 0 ; i != a.length ; ++i)
			if (a[i] != null) j++;
		Asteroid[] b = new Asteroid [j];
		for (int i = j = 0 ; i != a.length ; ++i)
			if (a[i] != null) b[j++] = a[i];
		return b;
	}

	// generate random masses in range
	private static double[] random_masses(int n, double min_mass,
	                                             double max_mass)
	{
		if (n < 0) throw new IllegalArgumentException();
		double[] masses = new double [n];
		for(int i = 0 ; i != n ; ++i)
			masses[i] = Math.random() * (max_mass - min_mass) + min_mass;
		return masses;
	}

	// create random circular orbits with non-overlapping distances
	private static Point[] random_orbits(double[] mass, double min_orbit,
	                                                    double max_orbit)
	{
		int n = mass.length;
		double[] radius = new double [n];
		for (int i = 0 ; i != n ; ++i) {
			double volume = mass[i] / Asteroid.density;
			radius[i] = Math.cbrt(0.75 * volume / Math.PI);
		}
		// random distances from the sun
		double[] orbit = new double [n];
		double orbit_range = max_orbit - min_orbit;
		for (int i = 0 ; i != n ; ++i) {
			orbit[i] = Math.random() * orbit_range + min_orbit;
			// check that the radii do not overlap
			for (int j = 0 ; j != i ; ++j) {
				double distance = Math.abs(orbit[i] - orbit[j]);
				if (distance < radius[i] + radius[j]) {
					i--;
					break;
				}
			}
		}
		// random initial angle
		Point[] orbit_pos = new Point [n];
		for (int i = 0 ; i != n ; ++i) {
			double len = orbit[i];
			double dir = Math.random() * 2 * Math.PI;
			orbit_pos[i] = new Point(len * Math.cos(dir),
			                         len * Math.sin(dir));
		}
		return orbit_pos;
	}

	// store points and masses
	private static void store(String path, Point[] positions,
	                          double[] masses) throws IOException
	{
		if (positions.length != masses.length)
			throw new IllegalArgumentException();
		// overwrite all contents of file
		FileOutputStream file = new FileOutputStream(path, false);
		for (int i = 0 ; i != masses.length ; ++i) {
			Point r = positions[i];
			String line = "" + r.x + ", " + r.y + ", " + masses[i] + "\n";
			file.write(line.getBytes());
		}
		file.close();
	}

	// load points from file
	private static Point[] load(String path, double[][] masses_addr)
	                            throws IOException
	{
		BufferedReader file = new BufferedReader(new FileReader(path));
		List <Point> point_list = new ArrayList <Point> ();
		List <Double> mass_list = new ArrayList <Double> ();
		String line;
		while ((line = file.readLine()) != null) {
			String [] xym = line.split(",");
			if (xym.length != 3)
				throw new IllegalStateException("Invalid state file");
			double x = Double.parseDouble(xym[0]);
			double y = Double.parseDouble(xym[1]);
			double m = Double.parseDouble(xym[2]);
			if (m <= 0.0)
				throw new IllegalStateException("Negative mass");
			point_list.add(new Point(x, y));
			mass_list.add(m);
		}
		file.close();
		if (point_list.size() < 2)
			throw new IllegalStateException("Less than 2 asteroids");
		double[] masses = new double [mass_list.size()];
		masses_addr[0] = masses;
		int i = 0;
		for (Double m : mass_list)
			masses[i++] = m;
		return point_list.toArray(new Point [0]);
	}

	// information for a push
	private static class Push {

		public final double energy;
		public final long game_time;

		public Push(double energy, long game_time)
		{
			this.energy = energy;
			this.game_time = game_time;
		}
	}

	// return the state of the game
	private static String state(String player_name,
	                            Orbit[] planets,
	                            Asteroid[] asteroids,
	                            List <Push> pushes,
	                            double sum_energy,
	                            long game_time,
	                            long cpu_time,
	                            int refresh)
	{
		StringBuffer buf = new StringBuffer();
		// compute mass ratio
		double sum_mass = 0.0;
		double max_mass = 0.0;
		for (Asteroid a : asteroids) {
			sum_mass += a.mass;
			if (max_mass < a.mass) max_mass = a.mass;
		}
		// header
		buf.append(game_time + ", " +
		           refresh + ", " +
		           asteroids.length + ", " +
		           planets.length + ", " +
		           pushes.size() + ", " +
		           human_power(sum_energy, 2) + ", " +
		           human_no_power(max_mass * 100.0 / sum_mass, 2) + ", " +
		           human_no_power(cpu_time / 1.0e9, 2) + ", " +
		           player_name);
		// solar system planets
		Point p = new Point();
		Point c = new Point();
		for (int i = 0 ; i != planets.length ; ++i) {
			Orbit o = planets[i];
			o.positionAt(game_time, p);
			o.center(c);
			double p_r = Math.sqrt(planet_radius[i] / planet_radius[0]) * 2;
			buf.append("\n" + planet_names[i]);
			buf.append(", " + planet_colors[i]);
			buf.append(", " + p.x + ", " + p.y + ", " + p_r);
			buf.append(", " + c.x + ", " + c.y);
			buf.append(", " + o.a + ", " + o.b);
		}
		// asteroids
		for (int i = 0 ; i != asteroids.length ; ++i) {
			Orbit o = asteroids[i].orbit;
			o.positionAt(game_time - asteroids[i].epoch, p);
			o.center(c);
			double p_r = Math.cbrt(asteroids[i].mass / Orbit.M) * 2;
			buf.append("\n, black");
			buf.append(", " + p.x + ", " + p.y + ", " + p_r);
			buf.append(", " + c.x + ", " + c.y);
			buf.append(", " + o.a + ", " + o.b);
		}
		// latest pushes
		int i = pushes.size() - 50;
		if (i < 0) i = 0;
		while (i != pushes.size()) {
			Push u = pushes.get(i++);
			buf.append("\n" + human_power(u.energy, 2) + ", " + u.game_time);
		}
		return buf.toString();
	}

	// orbits of solar system planets
	private static String[] planet_names =
	{"Mercury", "Venus", "Earth", "Mars", "Jupiter"};

	// colors for solar system planets
	private static String[] planet_colors =
	{"gray", "gold", "blue", "red", "brown"};

	// distances of solar system planets
	private static double[] planet_distances =
	{46.0e9, 107.5e9, 147.1e9, 206.7e9, 740.6e9};

	// periods of solar system planets
	private static double[] planet_periods =
	{87.969, 224.701, 365.256, 686.971, 4332.59};

	// radii of solar system planets
	private static int[] planet_radius =
	{2440, 6052, 6378, 3397, 71492};

	// compute orbits of known planets
	private static Orbit[] planet_orbits()
	{
		Orbit[] orbits = new Orbit [planet_names.length];
		for (int i = 0 ; i != orbits.length ; ++i) {
			double r1 = planet_distances[i];
			double T = planet_periods[i] * Orbit.dt;
			double w = 2 * Math.PI / T;
			double a = Math.cbrt(Orbit.GM / (w * w));
			double v1 = Math.sqrt(Orbit.GM * (2 / r1 - 1 / a));
			double r_dir = Math.random() * 2 * Math.PI;
			double v_dir = r_dir + Math.PI * 0.5;
			Point r = new Point(r1 * Math.cos(r_dir), r1 * Math.sin(r_dir));
			Point v = new Point(v1 * Math.cos(v_dir), v1 * Math.sin(v_dir));
			orbits[i] = new Orbit(r, v);
		}
		return orbits;
	}

	// serve static files and return dynamic file
	private static void gui(HTTPServer server, String content)
	                        throws UnknownServiceException
	{
		String path = null;
		for (;;) {
			// get request
			for (;;)
				try {
					path = server.request();
					break;
				} catch (IOException e) {
					System.err.println("HTTP request error: " + e.getMessage());
				}
			// dynamic content
			if (path.equals("data.txt"))
				// send dynamic content
				try {
					server.reply(content);
					return;
				} catch (IOException e) {
					System.err.println("HTTP dynamic reply error: " + e.getMessage());
					continue;
				}
			// static content
			if (path.equals("")) path = "webpage.html";
			else if (!path.equals("favicon.ico") &&
			         !path.equals("apple-touch-icon.png") &&
			         !path.equals("script.js")) break;
			// send file
			File file = new File(root + File.separator + "sim"
			                          + File.separator + path);
			try {
				server.reply(file);
			} catch (IOException e) {
				System.err.println("HTTP static reply error: " + e.getMessage());
			}
		}
		if (path == null)
			throw new UnknownServiceException("Unknown HTTP request (null path)");
		else
			throw new UnknownServiceException("Unknown HTTP request: \"" + path + "\"");
	}

	// recursive directory scan for files with given extension
	private static Set <File> directory(String path, String extension)
	{
		Set <File> files = new HashSet <File> ();
		Set <File> prev_dirs = new HashSet <File> ();
		prev_dirs.add(new File(path));
		do {
			Set <File> next_dirs = new HashSet <File> ();
			for (File dir : prev_dirs)
				for (File file : dir.listFiles())
					if (!file.canRead()) ;
					else if (file.isDirectory())
						next_dirs.add(file);
					else if (file.getPath().endsWith(extension))
						files.add(file);
			prev_dirs = next_dirs;
		} while (!prev_dirs.isEmpty());
		return files;
	}

	// last modified
	private static long last_modified(Iterable <File> files)
	{
		long last_date = 0;
		for (File file : files) {
			long date = file.lastModified();
			if (last_date < date)
				last_date = date;
		}
		return last_date;
	}

	// compile and load
	private static Class <Player> load(String group) throws IOException,
	                                       ReflectiveOperationException
	{
		String sep = File.separator;
		Set <File> source_files = directory(root + sep + group, ".java");
		File class_file  = new File(root + sep + group + sep + "Player.class");
		long player_modified = last_modified(source_files);
		long simulator_modified = last_modified(directory(root + sep + "sim", ".java"));
		if (!class_file.exists() || simulator_modified >= player_modified
		                         || player_modified >= class_file.lastModified()) {
			JavaCompiler compiler = ToolProvider.getSystemJavaCompiler();
			if (compiler == null)
				throw new IOException("Cannot find Java compiler");
			StandardJavaFileManager manager = compiler.
			                        getStandardFileManager(null, null, null);
			System.err.print("Compiling " + source_files.size() + " files ... ");
			if (!compiler.getTask(null, manager, null, null, null,
			     manager.getJavaFileObjectsFromFiles(source_files)).call())
				throw new IOException("Compilation failed");
			System.err.println("done!");
			class_file = new File(root + sep + group + sep + "Player.class");
			if (!class_file.exists())
				throw new FileNotFoundException("Missing class file");
		}
		ClassLoader loader = ToolProvider.getSystemToolClassLoader();
		if (loader == null)
			throw new IOException("Cannot find Java class loader");
		@SuppressWarnings("rawtypes")
		Class raw_class = loader.loadClass(root + "." + group + ".Player");
		@SuppressWarnings("unchecked")
		Class <Player> player_class = raw_class;
		return player_class;
	}

	// parse a percentage with 1 decimal point
	private static String human_ratio(double x)
	{
		if (x < 0.0) throw new IllegalArgumentException();
		if (x >= 1.0) return "100";
		int i = (int) (x / 0.001);
		int a =  i % 10;
		int b = (i / 10) % 10;
		int c = (i / 100);
		return c == 0 ? ("" + b + "." + a) : ("" + c + b + "." + a);
	}

	// parse a large real number and present in exponent mode
	private static String human_power(double x, int d)
	{
		if (x == 0.0) return "0";
		if (d <= 0) throw new IllegalArgumentException();
		StringBuffer buf = new StringBuffer();
		if (x < 0.0) {
			buf.append("-");
			x = -x;
		}
		int e = 0;
		double b = 1.0;
		while (b <= x) {
			b *= 10.0;
			e++;
		}
		b *= 0.1;
		int i = (int) (x / b);
		x -= b * i;
		buf.append(i + ".");
		do {
			b *= 0.1;
			i = (int) (x / b);
			x -= b * i;
			buf.append(i);
		} while (--d != 0);
		buf.append(" x 10^" + (e - 1));
		return buf.toString();
	}

	// parse a real number and cut the number of decimals
	private static String human_no_power(double x, int d)
	{
		if (x == 0.0) return "0";
		if (d < 0) throw new IllegalArgumentException();
		int e = 1;
		double b = 10.0;
		while (b <= x) {
			b *= 10.0;
			e++;
		}
		StringBuffer buf = new StringBuffer();
		do {
			b *= 0.1;
			int i = (int) (x / b);
			x -= b * i;
			if (e == 0) buf.append(".");
			buf.append(i);
		} while (--e != -d);
		return buf.toString();
	}
}
