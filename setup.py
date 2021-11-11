from setuptools import setup, find_packages
import os

with open(os.path.dirname(os.path.abspath(__file__)) + "/README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["ipython>=6", "nbformat>=4", "nbconvert>=5", "requests>=2"]

setup(
    name="notebookc",
    version="0.0.1",
    author="Anastasios Tzanidakis",
    author_email="atzanida@uw.edu",
    description="A pacakge on time series sampling",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/andytza/ermis/",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
