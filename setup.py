from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="chromatinfeaturesanalysis",
    version="0.1",
    description="A suite of scripts for analyzing DNA damage/repair/mutagenesis with respect to chromatin features like transcription factors and histone modifications",
    long_description_content_type="text/markdown",
    url='https://github.com/bmorledge-hampton19/benbiohelpers',
    author='Ben Morledge-Hampton',
    author_email='b.morledge-hampton@wsu.edu',
    license='MIT',
    python_requires='>=3.7',
    packages=["chromatinfeaturesanalysis"],
    package_dir={"chromatinfeaturesanalysis":"python_scripts"} # NOTE: This only works if the package is NOT installed as editable.
)