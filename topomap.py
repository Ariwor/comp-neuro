from brian2 import *
import numpy as np
import matplotlib.pyplot as plt

# Network parameters
Vr = 10 * mV
theta = 20 * mV
tau = 20 * ms
delta = 1.5 * ms
taurefr = 2 * ms

sparseness = 0.1
Je = .1 * mV
Ji = -.3 * mV

muext = 40 * mV
sigmaext = 2.1 * mV

# Neuron model
eqs = '''
dV/dt = (-V + muext + sigmaext * sqrt(tau) * xi) / tau : volt (unless refractory)
x : 1
y : 1
'''

n_cells = 12500  # Total number of cells
n_exc = int(0.8 * n_cells)  # 4:1 ratio for exc/inh
size = 1.  # Size of the network
simtime = 100 * ms  # Simulation time
sim_step = 1 * ms  # Display snapshots every sim_step ms
epsilon = 0.02  # Probability density
s_lat = 0.2  # Spread of the lateral connections
velocity = 0.3 * mm/ms  # velocity
max_distance = size * mm/np.sqrt(2)  # Since this is a torus
max_delay = max_distance / velocity  # Needed for the connectors

all_cells = NeuronGroup(n_cells, eqs, threshold='V > theta', reset='V = Vr', refractory=taurefr, method='euler')
all_cells.V = Vr
all_cells.x = size * np.random.rand(n_cells)
all_cells.y = size * np.random.rand(n_cells)

exc_cells = all_cells[:n_exc]
inh_cells = all_cells[n_exc:]

# Function returning the probabilities of connections as a function of distances
def probas(xi, xj):
    d1 = abs(xi - xj)
    min_d = np.minimum(d1, size - d1)
    distance = np.sqrt(np.sum(min_d**2, axis=1))
    return epsilon * np.exp(-distance**2 / (2 * s_lat**2))

# Function returning linear delays as function of distances
def delays(xi, xj):
    d1 = abs(xi - xj)
    min_d = np.minimum(d1, size - d1)
    distance = np.sqrt(np.sum(min_d**2, axis=1))
    return 0.1 * ms + (distance * mm) / velocity

# Define the network and add the groups and connections
net = Network(all_cells)

Ce1 = Synapses(exc_cells, exc_cells, model='''w : volt''',
               on_pre='V_post += w')
Ce1.connect(p='epsilon * exp(-sqrt((abs(x_pre - x_post) % size)**2 + (abs(y_pre - y_post) % size)**2) ** 2 / (2*s_lat**2))')
Ce1.w = Je
Ce1.delay = '0.1*ms + sqrt((abs(x_pre - x_post) % size)**2 + (abs(y_pre - y_post) % size)**2) * mm / velocity'
net.add(Ce1)

Ce2 = Synapses(exc_cells, inh_cells, model='''w : volt''',
               on_pre='V_post += w')
Ce2.connect(p='epsilon * exp(-sqrt((abs(x_pre - x_post) % size)**2 + (abs(y_pre - y_post) % size)**2) ** 2 / (2*s_lat**2))')
Ce2.w = Je
Ce2.delay = '0.1*ms + sqrt((abs(x_pre - x_post) % size)**2 + (abs(y_pre - y_post) % size)**2) * mm / velocity'
net.add(Ce2)

Ci1 = Synapses(inh_cells, inh_cells, model='''w : volt''',
               on_pre='V_post += w')
Ci1.connect(p='epsilon * exp(-sqrt((abs(x_pre - x_post) % size)**2 + (abs(y_pre - y_post) % size)**2) ** 2 / (2*s_lat**2))')
Ci1.w = Ji
Ci1.delay = '0.1*ms + sqrt((abs(x_pre - x_post) % size)**2 + (abs(y_pre - y_post) % size)**2) * mm / velocity'
net.add(Ci1)

Ci2 = Synapses(inh_cells, exc_cells, model='''w : volt''',
               on_pre='V_post += w')
Ci2.connect(p='epsilon * exp(-sqrt((abs(x_pre - x_post) % size)**2 + (abs(y_pre - y_post) % size)**2) ** 2 / (2*s_lat**2))')
Ci2.w = Ji
Ci2.delay = '0.1*ms + sqrt((abs(x_pre - x_post) % size)**2 + (abs(y_pre - y_post) % size)**2) * mm / velocity'
net.add(Ci2)

# Run simulation
net.run(simtime)

# Plotting results
plt.ion()  # To enter the interactive mode
print("Initializing the plots...")

fig, axs = plt.subplots(2, 1, figsize=(10, 8))

sc1 = axs[0].scatter(all_cells.x[:n_exc], all_cells.y[:n_exc], c='k', marker='.')
axs[0].set_title("Spikes")
plt.colorbar(sc1, ax=axs[0])

sc2 = axs[1].scatter(all_cells.x[:n_exc], all_cells.y[:n_exc], c='k', marker='.')
axs[1].set_title("Voltage")
plt.colorbar(sc2, ax=axs[1])

plt.show(block=False)  # Show the plot without blocking

print("Running network ...")
for t in range(int((simtime / sim_step) / ms)):
    net.run(sim_step)
    if plt.fignum_exists(fig.number):  # Check if the figure window is still open
        axs[1].cla()
        sc2 = axs[1].scatter(all_cells.x[:n_exc], all_cells.y[:n_exc], c=all_cells.V[:n_exc] / mV)
        axs[1].set_xlim(0, size)
        axs[1].set_ylim(0, size)
        axs[1].set_title(f"Voltage at t = {t * sim_step}")
        plt.draw()
        plt.pause(0.01)
    else:
        print("Figure window closed by user.")
        break

plt.ioff()  # To leave the interactive mode
if plt.fignum_exists(fig.number):  # Only show the final plot if the figure is still open
    plt.show()  # Display the final result

# Close the figure if it's still open
if plt.fignum_exists(fig.number):
    plt.close(fig)