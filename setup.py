import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CITE-seq-Count",
    version="1.4.2",
    author="Roelli Patrick",
    author_email="patrick.roelli@gmail.com",
    description="A python package that counts CITE seq or hashing data for single cell experiments",
    url="https://github.com/Hoohm/CITE-seq-Count/",
    packages=setuptools.find_packages(),
    entry_points={
          'console_scripts': [
              'CITE-seq-Count = cite_seq_count.__main__:main'
          ]
      },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    install_requires=[
          'python-levenshtein>=0.12.0',
          'scipy>=1.1.0',
          'multiprocess>=0.70.6.1',
          'umi_tools==1.0.0',
          'pytest==4.1.0',
          'pytest-dependency==0.4.0',
          'pandas>=0.23.4',
          'pybktree==1.1'
      ],
      python_requires='>=3.6'
)
