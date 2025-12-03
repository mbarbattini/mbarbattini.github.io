---
title: Numerically simulating the dc-SQUID equations of motion
date: 2025-12-3 11:00:00
---
Here we numerically simulate the dc-SQUID and plot the voltage across the junction and compare it to data presented in research papers.

The full code can be found on my Github repository [link](https://github.com/mbarbattini/computational_physics/tree/main/1_dc_squid_numerical)

### References

- _Quantum noise theory for the dc-SQUID_, Koch. 

    https://pubs.aip.org/aip/apl/article/38/5/380/525250/Quantum-noise-theory-for-the-dc-SQUID 

- _The SQUID Handbook_, John Clarke

- _Principles and Applications of SQUIDs_, John Clarke

## Introduction

A SQUID is the most sensitive detector of magnetic flux in the world. It has applications from biomagnetism where it can measure the electrical signals of brain cells, to superconducting qubits in quantum computing, and to measuring the magnetic field due to ancient rock samples deep in the Earth's crust.

## Theory

The story begins with a superconducting wire wound to form a closed loop. The area the loop makes defines a magnetic flux when an external magnetic field penetrates the loop. This is the same idea from classical electromagnetism. 

![](/images/dc_squid/superconducting_loop.png)

From DOI: 10.13140/RG.2.2.19768.29446

But quantum mechanics enters the story because the loop cannot allow any value of magnetic flux to penetrate the loop, it has to be quantized in units of the magnetic flux quantum \(\Phi_0=h/2e = 2.067\times 10^{-15}\ \ \mathrm{Wb}\). The reason is that in the superconducting state, the macroscopic wavefunction exists around the entire loop, and since the wavefunction must have a single value everywhere in the loop, this leads to the periodic condition

$$
\Phi=n\frac{h}{q}=n\Phi_0
$$
where \(q=2e\) is the charge of the particle, in this case a Cooper pair (two electrons with charge \(e$\)), and \(n\) is an integer. Therefore the magnetic flux quantum is \(\Phi_0=\frac{h}{2e}\)

### Josephson Junctions
If you add a small, non-superconducting barrier in-between a location in the superconducting loop, you create a Josephson junction. Here, the macroscopic wavefunction will tunnel through the small barrier. 

Brian Josephson discovered two equations, the DC equation which describes the value of the super-current across the junction

$$
I=I_c\sin(\varphi)
$$

where \(I_c\) is the critical current of the junction, and \(\varphi\) is the phase of the wavefunction across the junction,

and the AC equation

$$
\frac{d}{dt}(\Delta\varphi)=\frac{2eV}{\hbar}
$$

which says the time rate of change of the phase difference \(\Delta \varphi = \varphi_2-\varphi_1\) is related to the voltage across the junction.


### DC SQUID

If you add another JJ in the superconducting loop, there are now two elements that are disrupting the macroscopic wavefunction. Because of the sinusoidal nature of the current, the currents can exhibit constructive and destructive interference. This gives the SQUID it's name: superconducting quantum _interference_ device.



In the early days, researchers found that you needed to add a shunt resistance and capacitance on each branch in order to reduce hysteresis in the IV curve, so a practical SQUID uses the RCSJ model, which stands for "resistively and capacitively shunted junction". The circuit diagram is shown below.

![](/images/dc_squid/dc_squid_handbook.png)

From SQUID Handbook, John Clarke.


### Time-Domain Simulations


When simulating the dc-SQUID in the time domain, we need to work on a time scale that is determined by the characteristic frequency of the system \(\omega_c\). To get there, we start with the Josephson frequency \(\omega_J\) which is derived from substituting the AC Josephson equation into the DC Josephson equation

$$
I_s=I_c\sin(\frac{2eV}{\hbar} t)
$$

so the angular frequency term is \(\omega_J\)

$$
\omega_J=\frac{2eV}{\hbar}=\frac{2\pi}{\Phi_0}V
$$

When you apply a voltage \(V\) across a JJ you get a oscillating super-current at frequency \(\omega_J\). For example, a voltage of 1 \(\mu V\) gives a frequency of 

$$
\frac{\omega_J}{2\pi}=\frac{1 \times 10^{-6} \ \ \mathrm{V}}{2.067\times 10^{-15} \ \ \mathrm{Wb}} = 483.7 \ \ \mathrm{MHz}
$$

In this simulation, we can use a characteristic voltage \(V_c=I_0R\) where \(I_0\) is the critical current, which gives the _characteristic frequency_ \(\omega_c\) 

$$
\frac{\omega_c}{2\pi}= I_0 R / \Phi_0
$$

For typical devices, \(I_0 \approx 10 \ \mu \mathrm{A}\) and \(R\approx 20 \ \Omega\) which gives a voltage of 200 \(\mu \mathrm{V}\) so the typical frequency is \(\approx 96\ \mathrm{GHz}\). 

This explains that in a time-domain simulation we need to step forward on the order of \(\frac{1}{96 \ \mathrm{GHz}} \approx 2 \times 10^{-11} \ \mathrm{s}\) to capture how the system changes.


Since the SQUID can also be considered a parallel RC circuit, there is another characteristic frequency

$$
\omega_{RC} = \frac{1}{RC}
$$

and the relation for all the frequencies is captured in this equation
$$
\omega_p^2=\omega_c \omega_{RC}
$$

where \(\omega_P\) is the plasma frequency.

### Dimensionless parameters

When simulating SQUID equations, the values often get really large or small, so it is easier to work with dimensionless parameters. 

Some dimensionless SQUID parameters are:

- The hysteresis parameter, also known as the dimensionless inductance
$$
\beta = \frac{2LI_0}{\Phi_0}
$$

- The Stewart-McCumber parameter
$$
\beta_c = \bigg(\frac{\omega_c}{\omega_p}\bigg)^2 = \frac{2\pi}{\Phi_0} I_0 R^2 C
$$

### Differential equations 

The equations we will be simulating are presented in the Koch paper, which are a set of 4 coupled differential equations

$$
\frac{d}{dt} \delta_1 = \dot \delta_1
$$

$$
\frac{d}{dt} \dot \delta_1 = \frac{2\pi}{\Phi_0}\frac{1}{C}\bigg[\frac{I}{2} + \Gamma I_0 - I\sin(\delta_1)\bigg] - J 
$$

$$

\frac{d}{dt} \delta_2 = \dot \delta_2
$$

$$
\frac{d}{dt} \dot \delta_2 = \frac{2\pi}{\Phi_0}\frac{1}{C}\bigg[\frac{I}{2} + \Gamma I_0 - I\sin(\delta_2)\bigg] + J
$$

where \(\delta_1\) is the phase across the first Josephson junction, and \(\delta_2\) is the phase across the second Josephson Junction. 

The first and third equations are simple: the time rate of change of \(\delta_1\) is \(\dot \delta_1\), and then \(\dot \delta_1\) is determined by the second equation which actually has the physics in it. Then the same goes for \(\delta_2\). There are several terms to discuss:

- The noise parameter \(\Gamma\) is the ratio of the thermal energy to the Josephson energy
$$
\Gamma = \frac{E_{th}}{E_J}=\frac{k_B T}{I_0 \Phi_0 / 2\pi}=\frac{2\pi}{\Phi_0} \frac{k_B T}{I_0}
$$

So if the thermal energy is greater than the Josephson energy, then the noise parameter will be greater than 1. This would correspond to a resistive SQUID, i.e. one that is above its \(T_c\). For a SQUID operating at a typical temperature of ~ 4K, \(\Gamma < 0.1\).

- The DC Josephson equation shows up as the \(-I\sin(\delta_1)\) term, which makes the differential equations non-linear.

- The circulating current 
$$
J=\frac{I_0}{\pi\beta}\bigg(  \delta_1 - \delta_2 - \frac{2\pi \Phi_{ext}}{\Phi_0}\bigg)
$$

is the current that is generated as a consequence of the quantization of magnetic flux. When an external flux \(\Phi_{ext}\) is applied to the SQUID, the circulating current will be generated in order to keep the flux an integer multiple of \(\Phi_0\), akin to Lenz's law. The direction of the current is either clockwise or counterclockwise around the whole SQUID loop, depending on what part of the \(V-\Phi\) curve we are on. So this means that for the first junction the circulating current _adds_ to the current through the SQUID \(I\), but on the other junction it _subtracts_ from \(I\). 

To write the equations more compactly, we can put them in a vector known as the state vector. Then the time evolution of the state vector looks like this

$$
\frac{d}{dt}
\begin{pmatrix}
\delta_1 \\ 
\dot \delta_1 \\ 
\delta_2 \\ 
\dot \delta_2 \\
\end{pmatrix}
=
\begin{pmatrix}
\dot \delta_1\\
f_1 \\
\dot \delta_2 \\
f_2 \\
\end{pmatrix}
$$

where \(f_1\) and \(f_2\) are the right hand sides of the equations 2 and 4 above.

### Dependent Variables

The AC Josephson equation allows us to simulate something that has been measured by experimentalists: the voltage across the junction \(V\) (unlike the phase of a wavefunction which is not an observable). Since this equation is for a single junction, we need to add them together for 2 junctions

$$
\frac{d}{dt}\delta_1 + \frac{d}{dt}\delta_2 = \frac{2eV}{\hbar} + \frac{2eV}{\hbar}
$$

$$
\dot \delta_1 + \dot \delta_2 = \frac{4eV}{\hbar}
$$

$$
V = \frac{\hbar}{4e}(\dot \delta_1 + \dot \delta_2) 
$$

$$
V = \frac{h}{2\pi 4e}(\dot \delta_1 + \dot \delta_2) 
$$

and since the magnetic flux quantum is \(\Phi_0=\frac{h}{2e}\) we get

$$
V = \frac{\Phi_0}{4\pi}(\dot \delta_1 + \dot \delta_2) 
$$

So the voltage across the dc-SQUID is the sum of how the phases change in time, or the phase "velocities" if you use a kinematics analogy

## Numerical Simulations

This project is very easily implemented using `scipy.integrate.solve_ivp()` which uses the Runge-Kutta 4(5) method by default. It can easily accept a system of coupled differential equations by bunching the equations into a numpy array. The code for the equations of motion looks like this

```python
def dcSquidEquationOfMotion(self, t, y):
        
    delta1    = y[0]
    delta1dot = y[1]
    delta2    = y[2]
    delta2dot = y[3]

    # circulating current
    J = 1/(self.inductance*self.capacitance)*(delta1 - delta2 - 2*pi*self.flux)
    
    return np.array((
        delta1dot,
        (2*pi)/(phi0*self.capacitance)*(self.current/2 - self.criticalCurrent*sin(delta1) + self.gamma1*self.criticalCurrent) - J - 1/(self.resistance*self.capacitance)*delta1dot,
        delta2dot,
        (2*pi)/(phi0*self.capacitance)*(self.current/2 - self.criticalCurrent*sin(delta2) + self.gamma2*self.criticalCurrent) + J - 1/(self.resistance*self.capacitance)*delta2dot,
    ))
```

And when passed into the solver 
```python
sol = solve_ivp(
    self.dcSquidEquationOfMotion,
    [0, self.tMax/self.omegaC],
    state,
    method='RK45',
    max_step=1/self.omegaC*self.tStep)
```

the initial condition `state` for all simulations is an arbitrary state vector
$$
y(0)=
\begin{pmatrix}
\pi/4\\
0\\
\pi/6\\
0\\
\end{pmatrix}
$$



### IV Curve Simulation

> **Note**: the simulations in the next 2 sections both set the temperature of the SQUID to 0K, to observe the ideal behavior in the absense of thermal noise

The first thing I wanted to simulate was an IV curve. If everything works correctly, you should be able to simulate the superconducting state, i.e. there is zero voltage even when the current is non-zero. Then at some point the voltage comes back, which is the superconducting transition.

The function looks like this: 

At each time step 
1) Set the current through the SQUID \(I\) to a specific value
2) Simulate how the state vector changes in time with `solve_ivp()`
3) Use the final state as the initial condition to the next \(I\) value in the sweep

Here is the simulated IV curve when sweeping the current from -50 \(\mu \mathrm{A}\) to  +50 \(\mu \mathrm{A}\). The SQUID was set to the following parameters: \(I_0=10\ \mu \mathrm{A}\), \(L=120\ \mathrm{pH}\), \(R=20 \ \Omega\), \(C=0.08 \ \mathrm{pF}\)

![ivcurve](/images/dc_squid/iv_curve.png)

This confirms that the code is working correctly because the SQUID enters the resistive state when the bias current is 20 \(\mu A\) which is twice the critical current \(2I_0\), as expected because the current splits equally on each branch of the SQUID which each have a critical current of \(I_0\).

### V-\(\Phi\) curve

This function sweeps the external magnetic flux through the SQUID loop. The most notable feature is that the voltage across the SQUID is _periodic_ in multiples of the magnetic flux quantum \(\Phi_0\). This is the principles that allows the SQUID to be used to sense magnetic flux, i.e. it is a flux-to-voltage transformer. 

The function performs the same iterative procedure as above, but this time the bias current is set to a constant value and the external flux is swept.

Here is the plot for \(I_b=2.4I_0\) which is within the resisitive regime of the dc-SQUID

![](/images/dc_squid/v_phi_2p4.png)

The plot clearly shows a periodic oscillation with a period of 1 \(\Phi_0\). The maximums appear at \(n\Phi_0\) and the minimums appear at \((n+\frac{1}{2})\Phi_0\). 

The transfer function of the dc-SQUID is the slope of the \(V-\Phi\) curve 
$$
H=\bigg|\frac{\partial V}{\partial \Phi}\bigg|_{I=I_b}
$$

and the goal for using the device in an experiment is to have the largest transfer coefficient, which is at the locations \((n+\frac{1}{4})\Phi_0\). This gives the largest voltage output reading for some magnetic flux penetrating the loop.

Here is another plot for several different bias points. 

![](/images/dc_squid/v_phi_many.png)

At lower bias currents, there are large portions of the plot where the voltage is 0, which would not be good because then the magnetic flux would not be detected. It looks like the \(2.25I_0\) curve has the largest transfer coefficient, which is good, but there is little difference in the output voltage for other areas, which would not allow you to distinguish between certain values of magnetic flux.

Therefore the best bias point would be the \(2.5I_0\) curve because it gives a nice sinusoidal shape with relatively large transfer coefficient. This supports the idea that the SQUID should be biased in the resistive regime \(I_b > 2I_0\). 

Here is a plot from a research paper showing the same thing from measured data

![](/images/dc_squid/measured_v_phi_bolometer.png)

From [Dobbs et. al.](https://arxiv.org/pdf/1112.4215)




