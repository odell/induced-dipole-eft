# Atomic Units

I'm *still* trying to work out $X^+-H$ scattering in SI units. In atomic units, the Hamiltonian is simple
$$
H = -\frac{1}{2\tilde{\mu}}\frac{d^2}{d\tilde{r}^2} - \frac{\tilde{\alpha}}{2\tilde{r}^4}~,
$$
where a tilde indicates that the quantity is "dimensionless". Explicitly,
$$
\begin{split}
\mu & = \tilde{\mu} ~m_e ~,\\
r   & = \tilde{r} ~a_0 ~,\\
\alpha & =  ~\tilde{\alpha} ~a_0^3 ~,
\end{split}
$$
where $a_0$ is the Bohr radius and $m_e$ is the electron mass. Both are set to 1 in atomic units. As are $\hbar$, $e$, and $k_e=1/(4\pi\epsilon_0)$. ($a_0$ is a derived quantity that becomes 1 when the others are chosen to be 1.)



To make sense of the way I have written the interaction, $-C_4/r^4$, it would be helpful to have this Hamiltonian in SI units. So far, I have
$$
H = -\frac{\hbar^2}{2\mu}\frac{d^2}{dr^2} - \frac{\alpha}{4\pi\epsilon_0}\frac{1}{2r^4}~,
$$
where the polarizability $\alpha/(4\pi\epsilon_0) = 6.67\times 10^{-31}{\rm m}^3 = (9/2)a_0^3$. But this is clearly wrong, because the kinetic energy term has units of energy, $J$, and the potential has units of ${\rm m}^{-1}$... so where's my missing factor of $J\cdot m$?

## What do I know?

$$
\begin{split}
\hbar & = 1 \\
m_e   & = 1 \\
e     & = 1 \\
k_e=1/(4\pi\epsilon_0) & = 1
\end{split}
$$

$$
a_0 = \frac{(4\pi\epsilon_0) \hbar^2}{m_e e^2}
$$

$$
E_h = \frac{\hbar^2}{m_e a_0^2}
$$

Curly brackets indicate "the units of" the variable contained inside.
$$
\begin{split}
\{\alpha\} & = \frac{C\cdot {\rm m}^2}{V} = \frac{C^2 \cdot {\rm m}^2}{J} \\
\{\epsilon_0\} & = \frac{C^2}{J \cdot {\rm m}}
\end{split}
$$
And...
$$
\left\{\frac{\hbar^2}{\mu}\frac{d^2}{dr^2}\right\} = \{V(r)\} = J
$$