This folder contains the ParameterOptimization.py script to do parameter estimation. It uses the odeFunctions.py script as sourcefile for the functions.

The ParmOptim script should be used with system arguments for i) the path to the data file, ii) path to the ode file and iii) path to the output folder, and optionally iv) the model name.

Example: 
python ParameterOptimization.py "/path/to/input/CDDPdata4Model.csv" "/path/to/input/M019_20220311_ODE_p53model.txt" "/path/to/output" "M019"

# To use the script successfully, follow these steps:
1. Install anaconda (latest Python version that worked is 3.7.3)
2. Create a virtual environment and activate it
3. Run the following commands
	conda install numpy pandas scipy cython seaborn
	conda install -c conda-forge pydoe
	pip install sympy==1.1.1     

