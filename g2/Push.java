package pb.g2;

import pb.sim.Asteroid;

public class Push {
    Asteroid asteroid;
    int index;
    double energy;
    double direction;
    long expected_collision_time;

    Push(Asteroid asteroid, int index, double energy, double direction, long expected_collision_time) {
        this.asteroid = asteroid;
        this.index = index;
        this.energy = energy;
        this.direction = direction;
        this.expected_collision_time = expected_collision_time;
    }
}
