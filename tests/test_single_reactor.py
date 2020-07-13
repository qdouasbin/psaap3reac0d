import matplotlib.pyplot as plt
import toml

from psaap3reac0d import single_reactor

if __name__ == "__main__":
    # load TOML input file as a dictionary
    input_file = '../tests/BKD_LP4.toml'
    params = toml.load(input_file)

    # Conduct simulation
    result = single_reactor.single_reactor_simulation(params)

    # Plot results
    plt.figure()
    plt.plot(result['time'], result['T'])
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [K]')

    plt.figure()
    plt.plot(result['time'], 1e-5 * result['P'])
    plt.xlabel('Time [s]')
    plt.ylabel('Pressure [bar]')

    plt.show()
