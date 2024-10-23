Setup:
Download this repo.
Within IntelliJ, go to file->new->project from existing sources.
Select the downloaded wgd_abm folder. Accept all defaults to initialize projects.
Click on the drop-down menu in the top-right (near the green run button), and select "Edit Configurations...".
Click the "+" icon to add a new configuration and choose "Application".
Set the Main class to be ExampleOffLattice.
Add program arguments. The arguments format is as follows:
abm/parameters.csv output/ -v -i -o
The first argument is the path to a valid parameter file.
The second argument is the path to a folder where output will be saved. 
Arguments 3-5 are optional flags: -v enables visualisation, -i enables immune infiltration, -o enables output saving. 
You are now ready to run the model.

The parameter file values can be modified to control the behaviour of different celltypes in the model. 
The modifiable values are:
radius - cell radius
inhib_weight - larger values suppress division when cell is crowded
div_bias - higher values lead to faster dividing cells
r	g	b	- controls the colors for plotting each cell type
friction - values between 0 and 1. Controls how much cells move in response to being pushed by other cells.

