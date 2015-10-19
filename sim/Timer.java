package pb.sim;

import java.util.concurrent.*;
import java.lang.management.*;

class Timer extends Thread {

	private boolean start = false;
	private boolean finished = false;
	private Callable <?> task = null;
	private Exception exception = null;
	private Object result = null;

	public long cpu_ns()
	{
		return ManagementFactory.getThreadMXBean().getThreadCpuTime(getId());
	}

	public <T> T call(Callable <T> task, long timeout) throws Exception
	{
		if (!isAlive()) throw new IllegalStateException();
		if (task == null) throw new IllegalArgumentException();
		this.task = task;
		synchronized (this) {
			start = true;
			notify();
		}
		if (timeout <= 0)
			synchronized (this) {
				if (finished == false)
					try {
						wait();
					} catch (InterruptedException e) {}
			}
		else
			do {
				long cpu_time = cpu_ns();
				synchronized (this) {
					if (finished) break;
					try {
						wait(timeout);
					} catch (InterruptedException e) {}
				}
				cpu_time = cpu_ns() - cpu_time;
				timeout -= cpu_time / 1000000;
			} while (timeout > 0);
		if (finished == false)
			throw new TimeoutException();
		finished = false;
		if (exception != null) throw exception;
		@SuppressWarnings("unchecked")
		T result_T = (T) result;
		return result_T;
	}

	public void run()
	{
		for (;;) {
			synchronized (this) {
				if (start == false)
					try {
						wait();
					} catch (InterruptedException e) {}
			}
			start = false;
			exception = null;
			try {
				result = task.call();
			} catch (Exception e) {
				exception = e;
			}
			synchronized (this) {
				finished = true;
				notify();
			}
		}
	}
}
