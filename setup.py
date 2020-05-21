import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GasThermo",
    version="0.0.3",
    author="Robert F. De Jaco",
    author_email="dejac001@umn.edu",
    description="Gas-Phase Thermodynamics in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dejac001/GasThermo",
    packages=['GasThermo'],
    package_data={'GasThermo': [
        'cp_ig.csv', 'critical_constants.csv'
    ]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    package_dir={'': 'src'},
    python_requires='>=3.6',
    install_requires=["matplotlib==3.2.1", "numpy==1.18.4", "scipy==1.4.1", 'chem-util==0.0.4']
)
