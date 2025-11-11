"""
Setup script for pyJuTrack - Python wrapper for JuTrack.jl
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""

# Read version from pyJuTrack.py
def get_version():
    version = {}
    with open(os.path.join(os.path.dirname(__file__), 'pyJuTrack.py'), 'r') as f:
        for line in f:
            if line.startswith('__version__'):
                exec(line, version)
                return version['__version__']
    return "1.0.0"

setup(
    name='pyJuTrack',
    version=get_version(),
    author='JuTrack.jl Contributors',
    author_email='',
    description='Python wrapper for JuTrack.jl - A particle accelerator simulation and optimization library',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/MSU-Beam-Dynamics/JuTrack.jl',
    py_modules=['pyJuTrack'],
    python_requires='>=3.8',
    install_requires=[
        'juliacall>=0.9.0',
        'numpy>=1.20.0',
    ],
    extras_require={
        'plotting': ['matplotlib>=3.3.0'],
        'dev': [
            'pytest>=6.0',
            'matplotlib>=3.3.0',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    keywords='particle accelerator simulation optimization julia TPSA automatic-differentiation',
    project_urls={
        'Documentation': 'https://msu-beam-dynamics.github.io/JuTrack.jl/',
        'Source': 'https://github.com/MSU-Beam-Dynamics/JuTrack.jl',
        'Bug Reports': 'https://github.com/MSU-Beam-Dynamics/JuTrack.jl/issues',
    },
)
