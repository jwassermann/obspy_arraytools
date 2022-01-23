#!/usr/bin/env python

"""
The setup script.
"""

from setuptools import setup, find_packages

with open("README.md") as readme_file:
    readme = readme_file.read()

requirements = ["obspy"]

test_requirements = [
    "pytest>=3",
]

setup(
    author="Joachim Wassermann",
    author_email="j.wassermann@lmu.de",
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Lesser GNU General Public License v3 (LGPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="Array tools for ObsPy.",
    install_requires=requirements,
    license="Lesser GNU General Public License v3",
    long_description=readme,
    include_package_data=True,
    keywords="obspy_arraytools",
    name="obspy_arraytools",
    packages=find_packages(include=["obspy_arraytools", "obspy_arraytools.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/jwassermann/obspy_arraytools",
    version="0.1.0",
    zip_safe=False,
)
