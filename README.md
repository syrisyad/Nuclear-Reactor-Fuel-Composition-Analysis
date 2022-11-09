# Nuclear-Reactor-Fuel-Composition-Analysis
This is a Python program that simulates Uranium-Plutonium conversion chain in fuel burnup as a function of time using semi-implicit method. In a nuclear reactor with Uranium cycle, the changes in neutron from every isotope can be expressed by differential equations. These equations can be solved using semi-implicit differential method. In this program LMFBR type reactor, which is composed by Uranium-234, Uranium-235, and Uranium-238, is used. The isotope conversion model in this program is based on Chen's model (2012).

The Uranium cycle is simulated using Python with following algorithm:
1. Input fraction, time limit (in years), time step (in years), neutron flux, capture cross-section, absorption cross-section, decay constant, and initial number of neutrons.
2. Create the isotope matrix as a function of time.
3. Input the initial number of every isotope.
4. Iterate the semi-implicit reaction equations to every isotope as function of time.
5. Make a semilog plot for every isotope against time.

It needs to be noted that this program is solely based on Chen's model, which doesn't take Americium composition into consideration. This model also assumes that the initial composition of the reactor only consists of Uranium-238. Meanwhile, in a more realistic reactor, other compositions can affect some changes to the burnup. 
