# The purpose of this progrom is to install dependent libraries 
# Create by Steven Fuller
# Create date 4-4-2023
#

#importing libraries to install dependent libraries.
import subprocess
import sys

# Defining a function that will be used to install dependent libraries
# Function will take 1 argument that is the dependent library to install
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Calling function to install each dependent library
install("matplotlib")
install("numpy")
install("python-math")
install("tk")
