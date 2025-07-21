#!/usr/bin/env python3

from setuptools import setup, find_packages

with open("README", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="parfit",
    version="1.0.0",
    author="Federico Zahariev",
    author_email="federico.zahariev@gmail.com",
    description="ParFit automates the process of fitting molecular-mechanics parameters to data obtained by ab-initio calculations.",
    long_description=long_description,
    long_description_content_type="text/plain",
    url="https://github.com/fzahari/ParFit",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.16.0",
        "scipy>=1.3.0",
        "deap>=1.3.0",
        "matplotlib>=3.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black>=21.0",
            "flake8>=3.8",
        ],
    },
    entry_points={
        "console_scripts": [
            "parfit=ParFit.ParFit:main",
        ],
    },
    include_package_data=True,
    package_data={
        "ParFit": ["*.db", "*.prm"],
        "Mengine": ["*.h", "*.c", "makefile"],
        "Data": ["**/*"],
        "Doc": ["*.pdf", "*.doc"],
    },
)
