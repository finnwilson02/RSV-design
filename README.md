# RSV-design

Code used in the design of a satellite remote servicing vehicle (RSV). Modules should be able to easily be implemented for your own purposes as code is split into functions.

Link-budget includes code used to calculate the TCR (telemetry, control and ranging) requirements of the satellite, including gain, loss etc.

Optics looks at the requirements for an inspection camera to sufficiently survey client satellites.

Rad-intensity will provide information on solar radiation intensity (for solar panel power calcs) for a given orbit from classical orbital elements. To find your orbit or express it in classical orbital elements, you can refer to my orbital plots repo (https://github.com/finnwilson02/orbital-plots)

Orbital manoeuvres (credit Oscar Ansted 2021) was developed for the assignment to determine appropriate orbits, orbital manoeuvres and inspection flybys for the servicing mission.
