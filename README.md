# Multiple-Loads-Moving-on-Ice-Sheets
Code created for my master's thesis.
This code numerically computes the surface displacement a large ice sheet floating on water experiences when a load exerts force from above.
The code was written in MATLAB, and so needs to be run in an environment that supports .m files.

The file Constants.m contains parameter values experimentally measured at Lake Saroma in Japan, and may be modified to suit the user's needs.

The remaining files are seperated into two categories based on their output. The files with Anim in their name produce an animated plot, and saves a video of this. The files without Anim produce regular plots.

Besides the normally modular values for number of grid points, grid size, and time values, the main thing to vary is the path the load takes.
In LinearVelocity.m and AnimLinearVelocity.m the load must move with a constant linear velocity given by the constant vector "velocity".

In GeneralSolution.m, AnimGeneralSolution.m, GeneralSolutionMultiLoad.m, and AnimGeneralSolutionMultiLoad.m the load path is given in the function named "path". This needs to a vector parameterization of the path the load follows.
Due to the way MATLAB parses script all function definitions need to be at the end of the file, so to modify the path function you need to navigate down a bit, particularly in AnimGeneralSolutionMultiLoad.m.
If the path function is modified to depend on some values, those values need to be inputted into the function when it is called in the compute_eta_hat function located at the bottom of each file. In order for this function to accomplished, it also needs access to those values, so they need to be inputted into where it is called in the section of the code labeled % Computing the deflection.

In the files GeneralSolutionMultiLoad.m, and AnimGeneralSolutionMultiLoad.m there is also a variable for the amount of loads in the system, which may be modified.
By using if statements in the path function it is possible for different loads to follow different path.

All files except LinearVelocity.m can take a long time to compute, particularly if the number of grid points is large.
If the user desires the code to produce plots for multiple specific values of t, the sections % Computing the deflection and % Plotting can be wrapped in a for-loop. Reminder that in MATLAB -syntax a for-loop needs to be closed by an "end".

As demonstrated in AnimGeneralSolutionMultiLoad.m the view function can be used to specify the viewpoint of the plots produced.

The function named "load_pressure" may also be modified, but it needs to be a rotationally symmetric load distribution.
