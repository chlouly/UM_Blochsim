from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "UM_Blochsim",                # Module name
        sources=["blochsim.py"],  # Cython and C source files
        include_dirs=[np.get_include()],        # Add include directories if necessary
    )
]

setup(
    name="UM_Blochsim",
    ext_modules=cythonize(extensions),
    packages=["UM_Blochsim"],
    package_dir={"UM_Blochsim": "."},
    version="1.0",
    author="Christopher Louly",
    author_email="clouly@umich.edu",
    include_package_data=True,
    package_data={
        "UM_Blochsim": ["UM_Blochsim.so"],
    }
)
