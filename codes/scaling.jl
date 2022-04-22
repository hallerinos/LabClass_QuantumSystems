using PyPlot

Ns = 1:16
dim_wavefunction = 2 .^ Ns
dim_hamiltonian = dim_wavefunction .^ 2
gb = 8e+9
size_f64 = 64

plt.figure()
plt.plot(Ns, dim_wavefunction, label="wavefunction")
plt.plot(Ns, dim_hamiltonian, label="Hamiltonian")
plt.legend()
plt.yscale("log")
plt.xscale("linear")
plt.xlabel("N")
plt.ylabel("number of elements")
plt.savefig("scaling_analysis_dim.png", dpi=300)
plt.close()

memory_wavefunction, memory_hamiltonian = [dim*size_f64/gb for dim in [dim_wavefunction, dim_hamiltonian]]
plt.figure()
plt.plot(Ns, memory_wavefunction, label="wavefunction")
plt.plot(Ns, memory_hamiltonian, label="Hamiltonian")
plt.legend()
plt.yscale("log")
plt.xscale("linear")
plt.xlabel("N")
plt.ylabel("memory [GB]")
plt.savefig("scaling_analysis_memory.png", dpi=300)
plt.close()