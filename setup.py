import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CITE-seq-Count",
    version="1.3.2",
    author="Roelli Patrick",
    author_email="patrick.roelli@gmail.com",
    description="A small python package that deals with counting CITE seq or hashing data for single cell",
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
          'regex>=2018.07.11',
          'python-levenshtein>=0.12.0',
          'pandas>=0.23.3'
      ],
      python_requires='>3'
)