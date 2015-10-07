package pb.sim;

public interface Player {

	// the initial asteroids and their orbits and
	// the time limit of the game (negative if none)
	public void init(Asteroid[] asteroids,
	                 long time_limit);

	// asteroids -> the current asteroids and their orbits
	// energy    -> energy used to push asteroid (leave as 0 if no push)
	// direction -> direction of the push (ignored if 0 energy)
	public void play(Asteroid[] asteroids,
	                 double[] energy,
	                 double[] direction);
}
