---
title: Calculating Earth's Magnetic Field using Spherical Harmonic Interpolation
date: 2025-07-20 11:00:00
---


At the end of 2024 I found this [NASA White Paper](https://ntrs.nasa.gov/api/citations/20140007516/downloads/20140007516.pdf) which explained how to model the Earth's magnetic field using spherical harmonic interpolation. At the end of the paper the authors provided a challenge: use magnetic field data acquired by the MAGSAT spacecraft and calculate the spherical harmonic coefficients of the Earth's magnetic field up to 3rd order. To compare my solution, they also provided a dataset that contained their calculated coefficients.

I was familiar with magnetism from my undergraduate degree in physics, so this project seemed great because it was familiar but also challenging. I chose to use my strongest language, Python, to complete the project. I planned to use `pandas` to read and manipulate the data files, `matplotlib` for plotting 2D contour plots and 3D vector plots, and `scipy` for writing an optimization routine.

I enjoy projects where it involves making colorful plots of any sort of data and finding hidden detail in them. So I thought it would be cool to model Earth's magnetic field and perhaps be able to identify geographic information within it, such as the continents or metal deposits in the crust.

### Background Physics

The Ampere-Maxwell equation
$$
\nabla\times\mathbf{B} = \mu_0\mathbf{J} + \mu_0\epsilon_0\frac{\partial \mathbf{E}}{\partial t}
$$

shows that the magnetic field $\mathbf{B}$ is generated due to contributions from a current density $\mathbf{J}$ or a changing electric field $\mathbf{E}$. In this project, we will ignore Earth's electric field because its contribution to Earth's total magnetic field is much less than than the current density created by moving magma inside Earth's core. So the equation simplifies to

$$
\nabla \times \mathbf{B} = \mu_0 \mathbf{J}
$$

The goal is to solve for $\mathbf{B}$ at every point in some volume. This equation is hard to solve in its current form for multiple reasons. An improvement is to assume that you are in a source-free volume, which means that $\mathbf{J}=0$. The physical meaning of this assumption is that it is only valid to solve for the magnetic field _outside_ of the surface of the earth, where there are no currents of magma. The equation simplifies to

$$
\nabla \times \mathbf{B} = 0
$$

This equation is much easier to solve because it shows that the curl of the magnetic field is 0, which is the definition of a conservative vector field, i.e. one that only depends on the position rather than the path taken. Then you can conveniently say that the magnetic field can be calculated as the gradient of a scalar potential $V$

$$
\mathbf{B} = -\mu_0 \nabla V
$$

Using the other Maxwell equation for magnetic fields (which says that magnetic fields have no sources or sinks)

$$
\nabla \cdot \mathbf{B} = 0
$$

we can substitute the expression to show that the scalar potential must satisfy Laplace's equation

$$
\nabla \cdot (-\mu_0 \nabla V) = 0 \\
\nabla^2 V = 0 \\
$$

Solving Laplace's equation in spherical coordinates is a popular exercise that takes a bit of mathematical rigor, so I'm going to skip it in this writeup. The main result is that the scalar potential $V$ can be calculated as a series of spherical harmonic functions

$$
V(r,\theta,\phi)=\frac{a}{\mu_0}\sum_{l=1}^{N}\big(\frac{a}{r}\big)^{l+1}\sum_{m=0}^{l}(g_{lm}\cos m\phi + h_{lm}\sin m\phi) P_{l}^m(\cos\theta)
$$

where $P_{l}^m(\cos\theta)$ are the associated Legendre polynomials. We sum over two indices: $l$ and $m$. The coefficients we need to solve for are $g_{lm}$ and $h_{lm}$, up to degree $N=3$, which means there are 15 total coefficients to be determined. Increasing the degree means that the solution is more accurate, but adds many more coeffients in the expansion.


| $l$ | $m$ | $g_{l,m}$ | $h_{l,m}$ |
| --- | --- | --------- | --------- |
| 1 | 0 | $g_{10}$ | Not used |
| 1 | 1 | $g_{11}$ | $h_{11}$ |
| 2 | 0 | $g_{20}$ | Not used |
| 2 | 1 | $g_{21}$ | $h_{21}$ |
| 2 | 2 | $g_{22}$ | $h_{22}$ |
| 3 | 0 | $g_{30}$ | Not used |
| 3 | 1 | $g_{31}$ | $h_{31}$ |
| 3 | 2 | $g_{32}$ | $h_{32}$ |
| 3 | 3 | $g_{33}$ | $h_{33}$ |



The $l$ terms correspond to the multipole expansion. $l=1$ is the dipole term, $l=2$ is the quadrupole term, $l=3$ is the octupole term, and so on. There are no $l=0$ terms included in the expansion because this would correspond to a magnetic monopole, which is not allowed by Gauss's law for magnetism $\nabla \cdot \mathbf{B} = 0$.

So to find the magnetic field, we need to take the gradient of this expression for each component separatley, keeping in mind that we are in spherical coordiantes. I also learned that in geomagnetism there are certain convetions for the directions of the magnetic field, which is not something I had dealt with before.

$$
X = -B_\theta = \mu_0 \frac{1}{r} \frac{\partial V}{\partial \theta} \ \ \ \ \ \ (\mathrm{northward \ \ component}) \\
$$

$$
Y = +B_\phi = -\mu_0 \frac{1}{r \sin \theta} \frac{\partial V}{\partial \phi} \ \ \ \ \ \ (\mathrm{eastward \ \ component}) \\
$$

$$
Z = -B_r = \mu_0 \frac{\partial V}{\partial \theta} \ \ \ \ \ \ (\mathrm{downward \ \ component} ) \\
$$

Here is an example of the first few expressions for the X components of the field. 

$$
X = -\mu_0\frac{1}{r}\frac{\partial V}{\partial \theta} = -\frac{a^3}{r^3}\sin(\theta) g_1^0+\frac{a^3}{r^3}\cos(\phi)\sin(\theta)g_1^1+\frac{a^3}{r^3}\sin(\phi)\cos(\theta)h_1^1 - \frac{3}{2}\frac{a^4}{r^4}\sin(2\theta)g_2^0 + ...
$$

This was manageable to do up to $N=3$. However this gets tedious quickly when trying to caluclate the higher order terms, because you have to take derivatives of the Legendre polynomials which become quite complicated like

$$
P_4^2(\cos(\theta)) = -\frac{\sqrt 5}{4}\bigg[-7\sin^2(\theta) \cos^2(\theta)+\sin^2(\theta)\bigg]
$$

The matrix equation for this problem is

$$
\mathbf{Ag}=\mathbf{b}
$$

where $\mathbf{b}$ is a column vector of the magnetic field measurements and $\mathbf{g}$ is the column vector of the coefficients. The coefficients can be solved like

$$
\mathbf{A}^{-1}\mathbf{Ag}=\mathbf{A}^{-1}\mathbf{b} \\
$$

$$
\mathbf{g}=\mathbf{A}^{-1}\mathbf{b}
$$

A key note is that the matrix $\mathbf{A}$ is not square, so we cannot calculate its inverse directly. This is because there are are more magnetic field measurements (300) than coefficients to solve for (15). So instead, you need to perform the Moore-penrose pseudo-inverse which works for non-square matrices. This algorithm is provided in `scipy.linalg` and is simple as `pinv(A)`.


### Solving the Problem with Python

My code would need to do the following:

1. Read in the magnetic field data collected by the satellite
2. Create a matrix for the coefficients
3. Solve for the coefficients using the psuedo-inverse algorithm
4. Use the coefficients to calculate the magnetic field at any point in space

The data provided were magnetic field measurements from the NASA MAGSAT spacecraft that operated around 1980. Here is what the provided `field.dat` data looks like as a 3D vector plot

![](images/sph-harm-geomagnetic/dataset_vector_field.gif)


Each of these vectors is a magnetic field reading from the satellite. It's interesting to see that the data covers the entire surface of a sphere which shows how the satellite travels across the entire Earth as it orbits. Also, visualizing the field data like this was good confirmation that I processed the data correctly. You can see at the end of the gif the vectors seem to be pointing from left to right which indicates a magnetic field with no divergence (as required by Maxwell's equations). Since the orbit of the satellite is much greater than the radius of the Earth, this would be considered in the far-field and the dominant term of the magnetic field should be the dipole term, and this data clearly shows the magnetic field due to a dipole.

Now it was time to create the coefficient matrix. I decided to normalize the expressions in the radial direction which means setting $a/r = 1$ . This means that the interpolation sphere is the radius of the Earth. If I kept this in, I could calculate what the magnetic field is like at some altitute above the surface, which might be useful for a plane that needs to know the Earth's magnetic field when it's in the air.

Here's what the expressions look like in code for the X coefficients up to $N=3$:

```python
x_coe = np.array([
    -sin(theta),                                                              # g10
    +cos(phi)*cos(theta),                                                     # g11
    -3*cos(theta)*sin(theta),                                                 # g20
    +3*cos(2*theta)*cos(phi),                                                 # g21
    +3*cos(2*phi)*sin(2*theta),                                               # g22
    +(-(15/2)*cos(theta)**2*sin(theta)+(3/2)*sin(theta)),                     # g30
    +((sqrt(6)/4)*cos(phi)*(15*cos(theta)**3-11*cos(theta))),                 # g31
    +(sqrt(15)/2)*(cos(2*phi)*(2*sin(theta)*cos(theta)**2-sin(theta)**3)),    # g32
    +(cos(3*phi)*(3*sqrt(10)/4)*sin(theta)**2*cos(theta)),                    # g33
    +sin(phi)*cos(theta),                                                     # h11
    +3*sin(phi)*cos(2*theta),                                                 # h21
    +3*sin(2*phi)*sin(2*theta),                                               # h22
    +((sqrt(6)/4)*sin(phi)*(15*cos(theta)**3-11*cos(theta))),                 # h31
    +(sqrt(15)/2)*(sin(2*phi)*(2*sin(theta)*cos(theta)**2-sin(theta)**3)),    # h32
    +(sin(3*phi)*(3*sqrt(10)/4)*sin(theta)**2*cos(theta)),                    # h33
]).T
```
This was messy but manageable for the problem.

Here are the results of my work. All plots that follow show the magnetic field magnitude as filled in contours. My results are good considering that the magnitude of the plot is on the same order as 50 $\mu T$ which are the values I've measured for Earth's field in the lab with a magnetometer.

![](images/sph-harm-geomagnetic/final_mine.png)

Here is the error from my calculation and the official coefficients from the NOAA model for the year 1980:

| Coefficient | NOAA 1980 [nT] | My Program [nT]| Error [%]
| ----------- | -------------- | -------------- | ------- |
$g_{10}$ | -29992 | -24768 | 17.42
$g_{11}$ | -1956 | -1792 | 8.34
$g_{20}$ | -1997 | -1285 | 35.64
$g_{21}$ | 3027 | 1454 | 51.96
$g_{22}$ | 1663 | 363 | 78.15
$g_{30}$ | 1281 | 816 | 36.28
$g_{31}$ | -2180 | -1357 | 37.71
$g_{32}$ | 1251 | 964 | 22.90
$g_{33}$ | 833 | 65 | 92.19
$h_{11}$ | 5604 | 4419 | 21.14
$h_{21}$ | -2129 | -1023 | 51.94
$h_{22}$ | -200 | -57 | 71.32
$h_{31}$ | -336 | -154 | 54.13
$h_{32}$ | 271 | 218 | 19.38
$h_{33}$ | -252 | -163 | 35.07

Since the error calculations are close for some of the coefficients like $g_{11}$, I think I did the problem correctly. I think the reason that some of the coefficients have large error is because the official NOAA model uses much more magnetic field measurements than the 100 measurements in my dataset, and they use a higher order in the interpolation like $N=15$ or $N=20$.

Another feature that was missing from my plot was the South Atlantic Anomaly, which is an area near Argentina and Brazil with an unusually low magnetic field magnitude. Here is the plot from the white paper that I was referencing that shows what the plots should look like:

![](images/sph-harm-geomagnetic/field_from_white_paper.png)


So to verify if my plotting code was working correctly, I decieded to use the official coefficients from 1980 with the same contour levels I specified in my plot

![](images/sph-harm-geomagnetic/final_1980.png)

Now you can see the South Atlantic Anomaly which is good. Looking closer, it looks the the magnitude of the field is actually reaching a maximum in this area, but it should be reaching a lower magnitude of 20 - 30 $\mu T$ like it shows in the official contour plot. I have no idea why this is happening. It might be because I am not including the radial terms in my calculation, and the south atlantic anomaly might be at some altitude above the Earth.

It looked like my work had a different overall field magnitude near the South Atlantic Anomaly, so here is my plot with different contour levels that highlight it better:

![](images/sph-harm-geomagnetic/final_mine_better_contours.png)

This allows you to see a little more detail and specifically looking at the shape of the lines over North America looks quite similar to the official plot. 

So it looks like my approach is simply limited by the accuracy of the coefficients, which comes from having more magnetic field data and using more terms in the expansion.

I found this research paper https://doi.org/10.1186/s40623-020-01252-9 where they used $N=20$ and were able to see slight differences in the Earth's field from 2014 compared to 2020. To highlight the South Atlantic anomaly, they used contour levels (white lines) in steps of 500 nT from 22,500 nT to 27,000 nT, which is twice as detailed as my plots which had steps of 1 $\mu T$ in areas of interest.

![](images/sph-harm-geomagnetic/south_atlantic_anomaly.png)

From https://en.wikipedia.org/wiki/South_Atlantic_Anomaly.


Unfortunately, if I wanted to do $N=20$ in my calculation, this would correspond to 440 terms in the expansion. This is larger than the number of satellite magnetic field measurements provided by the dataset (300), so my matrix would be underdetermined and I wouldn't get an accurate answer. If I wanted to do this, I would need many more magnetic field measurements, and in this paper they likely used 10,000+ measurements because they had satellite data over a 21 year period.



### Future Work

- Find another data set with more magnetic field measurements
- Compare my code to the official C code from the National Oceanic Atmospheric Assoociation (NOAA), the World Magnetic Model: [world-magnetic-model](https://www.ncei.noaa.gov/products/world-magnetic-model)
- Figure out how to interpolate more correctly near the Equator, which has values that blow up to infinity because there are terms like $1/\sin(\phi)$ and the equator corresponds to $\phi = 0$ degrees so $1/\sin(0)=\infty$
- Go higher than $N=3$, which would require an automated algebra system to calculate the terms using something like [sympy](https://www.sympy.org/en/index.html)


