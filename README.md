# Welcome to the Course: [Introduction to the Theory of Quantum Systems](https://wwwfr.uni.lu/recherche/fstm/dphyms/research/theory_of_mesoscopic_quantum_systems)!

During this two-day lab class, you will learn basic analytic and computational methods we use in the [Theory of Mesoscopic Quantum Systems](https://wwwfr.uni.lu/recherche/fstm/dphyms/research/theory_of_mesoscopic_quantum_systems) group.
On both days, you will receive lectures in the morning, followed by programming sessions in the afternoon.
We focus on two important techniques, exact diagonalization (ED) and matrix product states (MPS).

You will need to develop your own codes in order to investigate an assignment, and you may choose one of three subjects, with increasing complexity
1. The Ising model: order and disorder in quantum magnets (easy)
2. The Bose-Hubbard model: superfluids and Mott-insulators (intermediate)
3. Topological quantum matter: exotic bound modes and non-local order parameters (advanced)

You need to tell me which topic you like after the first day, based on which I prepare a lecture and a specific assignment for the second day.

# Exact Diagonalization
This technique is a brute-force approach to numerical simulations of condensed matter systems, because it constructs the Hamiltonian operator in a suitable vector space, followed by a straightforward diagonalization in order to access the spectrum.

Mathematically, exact diagonalization corresponds to finding solutions of the stationary Schrödinger equation in Dirac notation
$$
\hat H \ket{\psi_n} = E_n\ket{\psi_n}
$$
in which $\hat H$ corresponds to the Hamiltonian operator and $\ket{\psi_n}\in\mathcal H$ is the corresponding eigenstate of $\hat H$ with eigenvalue $E_n$.
Note that the Hamiltonian can be represented as a matrix in an arbitrary basis $\{\phi_n\}$, with matrix elements $H_{mn}^\phi = \braket{\phi_m | \hat H|\phi_n}$.
The Hamiltonian in this notation is simply 
$$
\hat H = H_{mn}^\phi\ket{\phi_m}\bra{\phi_n}.
$$
Please note that I use the sum convention.
The matrix elements depend of course on the choice of the basis, and $H_{mn}^\psi = E_n\delta_{mn}$ is diagonal.
In a typical situation, the Hamiltonian is written in a canonical basis (e.g. in a real-space basis where $\ket{\phi_n}$ is localized at position $n$), such that the matrix elements $H_{mn}^\phi$ are known, but the eigenvalues $E_n$ are not.

> Convince yourself that the eigendecomposition of $H^\phi_{mn}$ leads to the spectral decomposition of $\hat H$.

A diagonalization of the Hamiltonian matrix $H^\phi = U^\phi D U^{\phi\dag}$ will always lead to a unitary transformation of the form
$$
\hat H = (U^\phi D U^{\phi\dag})_{mn}\ket{\phi_m}\bra{\phi_n}
= U^\phi_{mk} D_k U^{\phi*}_{nk}\ket{\phi_m}\bra{\phi_n}
= D_k\ket{\psi_k'}\bra{\psi_k'}
$$
A comparison with the Schrödinger equation then leads to $D_k=E_k$ and $\ket{\psi_k}=\ket{\psi'_k}$ up to an irrelevant arbitrary phase.

Now, if we want to measure eigenstate observables in the canonical basis, they can be expressed in terms of the unitary transformation which diagonalizes $H^\phi$.

> Find a formula for the expectation value of a general Hermitian operator $\braket{\psi_i|\hat O|\psi_j}$. Double-check your result by computing $\braket{\psi_i|\hat H|\psi_j}$.

The expectation value of any operator, initially expressed in the canonical basis, can then be expressed in terms of the unitary rotation of $H$ to the diagonal basis, i.e.
$$
\braket{\psi_i|\hat O|\psi_j} = O^\phi_{mn}\braket{\psi_j|\phi_m}\braket{\phi_n|\psi_j} = (U^{\phi\dag} O^\phi U^\phi)_{ij}
$$

That's about all we need to get started.
So let's move to a particular example and figure out the missing steps along the way.

## The transverse-field Ising model (TFIM)
The TFIM is the quantum analog of the classical (Ernst) Ising model, motivated from the study of ferromagnetism in statistical mechanics.

The Hamiltonian reads
$$
\hat H = J\sum_{\braket{i,j}}\hat S^z_{\bm r_i}\hat S^z_{\bm r_j} + h \sum_{i}\hat S^x_{\bm r_i}
$$
in which $\bm r_i$ denotes some position in a Bravais lattice, and $\hat S^k_{\bm r_i}$ is the local spin-1/2 operator along the $k$ direction.
J denotes the binding energy between adjacent spins ($\braket{i,j}$ denotes a constrained sum over nearest neighbors only), and $h$ is the magnetic field amplitude associated with the Zeeman coupling.
The spin operators satisfy the same algebraic commutation relations as orbital angular momentum, i.e.
$$
[\hat S^i,\hat S^j] = i\hbar\varepsilon_{ijk}\hat S^l
$$
and they can be represented in terms of the Pauli matrices
$$
\hat S^i \simeq \frac{i\hbar}2\sigma^i,
$$

$$
\sigma^0 =
\begin{pmatrix}
       1 & 0 \\
       0 & 1
\end{pmatrix}
,\
\sigma^x =
\begin{pmatrix}
       0 & 1 \\
       1 & 0
\end{pmatrix}
,\
\sigma^y =
\begin{pmatrix}
       0 & -i \\
       i & 0
\end{pmatrix}
,\
\sigma^z =
\begin{pmatrix}
       1 & 0 \\
       0 & -1
\end{pmatrix}
.
$$
Keep in mind that the Pauli matrices are a good basis for any $2\times2$ matrix, it might be imortant for some of the assignments.

In the following, we set $\hbar=1$, and remember that we thereby changed our unit of energy to a reduced Planck quantum.
Furthermore, to make our discussion easier for now, we consider a one-dimensional chain such that $\bm r_i = ia$, and set the lattice constant $a=1$.
We don't need to be careful about that because $a$ does not enter in $\hat H$, although in an actual experiment, $J$ will of course decrease if the distance between the spins is increased.

Remember that a general many-body state, composed by $N$ instances of a local Hilbert space $\ket{i_n}\in\mathcal H_n$, reads
$$
\ket\Psi = \sum_{\{i_n\}}C_{i_1,i_2,\dots,i_n}\ket{i_1,i_2,\dots,i_n}
$$
where
$$
\ket{i_1,i_2,\dots,i_n} = \ket{i_1}\otimes\ket{i_2}\otimes\dots\otimes\ket{i_n}
$$
denotes the Kronecker product.
In order to construct the local spin operators of $\hat H$, we apply the same strategy:
$$
\hat S^k_n \simeq \left(\bigotimes_{j<n}\sigma^0\right)\sigma^k\left(\bigotimes_{j>n}\sigma^0\right)
$$

> Construct the operators $\hat S^k_{n}$ in a `julia` script and demonstrate that $[\hat S^k_{m}, \hat S^l_{n}]=0\ \forall m\neq n$.
> Use the operators to construct the Hamiltonian.

If we now consider a canonical basis on every instance to be $i_n\in\{\uparrow,\downarrow\}$ in the $\hat S^z$ basis, we see that the Hilbert space has dimension $2^N$.
> Consider $N=8,16,32,64$ instances and assume that each coefficient of $C_{\{i_n\}}$ is a `float64`. How much memory is needed to store the wavefunction? How much memory do you need for a naive matrix expression of $H$? What is the maximum $N$ that you can simulate on your laptop? What if the Hamiltonian is perturbed by just a single term $\hat S^y_{\bm r_i}$?

## The Lanczos algorithm