"""Setup.py file
"""
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CITE-seq-Count",
    version="2.0.0",
    author="Roelli Patrick",
    author_email="patrick.roelli@gmail.com",
    description="A python package to map reads from CITE-seq or hashing data for single cell experiments",
    url="https://github.com/Hoohm/CITE-seq-Count/",
    download_url="https://github.com/Hoohm/CITE-seq-Count/archive/2.0.0.tar.gz",
    packages=setuptools.find_packages(),
    entry_points={"console_scripts": ["CITE-seq-Count = cite_seq_count.__main__:main"]},
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    install_requires=[
        "scipy>=1.1.0",
        "umi_tools==1.1.4",
        "pytest>=6.0.0",
        "pytest-dependency==0.4.0",
        "cython>=0.29.17",
        "pyyaml==6.0",
        "pooch==1.6.0",
        "six==1.16.0",
        "polars==0.20.3-rc2",
    ],
    python_requires="==3.11.6",
    package_data={"report_template": ["templates/*.json"]},
)
