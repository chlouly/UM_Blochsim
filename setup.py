from setuptools import setup, Extension
import numpy as np

extensions = [
    Extension(
        "UM_Blochsim",                # Module name
        sources=["blochsim.py", "blochsim_core.c"],  # Cython and C source files
        include_dirs=[np.get_include()],        # Add include directories if necessary
    )
]

setup(
    name="UM_Blochsim",
    ext_modules=extensions
)
