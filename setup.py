from pkg_resources import parse_requirements
import setuptools

setuptools.setup(
    name="PyEnzymeKinetics",
    version="1.1.6",
    url="https://github.com/haeussma/PyEnzymeKinetics",
    author="Haeussler, Max",
    author_email="max.haeussler@gmail.com",
    packages=setuptools.find_packages(),
    install_requires=["lmfit", "matplotlib", "scipy", "numpy", ]
)