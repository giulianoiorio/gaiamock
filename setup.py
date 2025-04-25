from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import subprocess
import os
import sys
from pathlib import Path

class CustomBuild(build_py):
    def run(self):
        # Check if pkg-config is installed
        try:
            subprocess.check_output(['pkg-config', '--exists', 'gsl'])
        except subprocess.CalledProcessError:
            print("ERROR: GSL development libraries not found. Please install libgsl-dev (Linux) or gsl (Mac with Homebrew).", file=sys.stderr)
            sys.exit(1)

        # Get GSL flags
        cflags = subprocess.check_output(['pkg-config', '--cflags', 'gsl']).decode().strip().split()
        libs = subprocess.check_output(['pkg-config', '--libs', 'gsl']).decode().strip().split()


        # Destination for compiled shared library
        build_lib = Path(self.build_lib) / "gaiamock"
        build_lib.mkdir(parents=True, exist_ok=True)
        libname = "kepler_solve_astrometry.so"
        out_path = build_lib / libname

        # Build the shared library
        subprocess.check_call(
            ["gcc", "-shared", "-fPIC", "gaiamock/src/kepler_solve_astrometry.c", "-o", str(out_path)] + cflags + libs
        )

        super().run()

# Setup script
setup(
    name="gaiamock",
    version="1.1.0",
    description="A package for simulatiing Gaia astrometry at epoch level",
    author="Kareem El-Badry",
    author_email="",
    packages=find_packages(),
    cmdclass={"build_py": CustomBuild},
    include_package_data=True,
    package_data={"gaiamock": ["kepler_solve_astrometry.so"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">3.8",
)