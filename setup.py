import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SciThermo",
    version="0.0.1",
    author="Robert F. De Jaco",
    author_email="dejac001@umn.edu",
    description="Scientific Thermodynamics in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dejac001/SciThermo",
    packages=['scithermo'],
    package_data={'scithermo': [
        'cp_ig.csv', 'critical_constants.csv'
    ]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["matplotlib==3.2.1", "numpy==1.18.2", "scipy==1.4.1"]
)
