from distutils.core import setup

#This is a list of files to install, and where
#(relative to the 'root' dir, where setup.py is)
#You could be more specific.
files = ["color_scales/*", "pysolo_config.xml"]

setup(name = "pysolo",
    version = "0.4",
    description = "Display CfRadial format radar/lidar file",
    author = "Xin Pan and Joe VanAndel",
    author_email = "vanandel@ucar.edu",
    url = "whatever",
    #Name the folder where your packages live:
    #(If you have other packages (dirs) or modules (py files) then
    #put them into the package directory - they will be found 
    #recursively.)
    packages = ['PySolo'],
    #'package' package must contain files (see list above)
    #This dict maps the package name =to=> directories
    #It says, package *needs* these files.
    package_data = {'PySolo' : files },
    #'start_pysolo' is in the root.
    scripts = ["start_pysolo"],
    long_description = """Really long text here.""" 
    #
    #This next part it for the Cheese Shop, look a little down the page.
    #classifiers = []     
) 
